package beast.gss;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.NSLogger;
import beast.core.Operator;
import beast.core.State;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.Parameter;
import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.NSLogAnalyser;
import beast.util.Randomizer;


/* From wikipedia:
Start with  N points θ1,...,θN sampled from prior.
 for i = 1 to j do        % The number of iterations j is chosen by guesswork.          
 Li := min(current likelihood values of the points)
 Xi := exp(-i/N);
 wi := X{i-1}-X{i}
 Z  := Z+Li*wi
     Save the point with least likelihood as a sample point with weight wi
     Update the point with least likelihood with some Markov chain Monte Carlo steps according to the prior, accepting only steps that
     keep the likelihood above L_i.
 end
 return Z    
*/

@Description("Nested sampling for phylogenetics")
@Citation("Patricio Maturana, Brendon J. Brewer, Steffen Klaere, Remco Bouckaert. Model selection and parameter inference in phylogenetics using Nested Sampling. arXiv:1703.05471v3")
public class NS extends MCMC {
	final static int SAMPLE_COUNT = 100;
	
	public Input<Integer> particleCountInput = new Input<>("particleCount", "number of particles (default 100)", 100);
	public Input<Integer> subChainLengthInput = new Input<>("subChainLength",
			"number of MCMC samples for each epoch (default 1000)", 1000);

	final public Input<Distribution> samplingDistributionInput = new Input<>("samplingDistribution",
			"probability distribution to sample from (e.g. a prior). "
					+ "If not specified, everything but the likelihood will be used as sampling distribution.");

	public Input<Double> epsilonInput = new Input<>("epsilon", "stopping criterion: smallest change in ML estimate to accept", 1e-6);
	public Input<Double> stopFactorInput = new Input<>("stopFactor", "stopping criterion: use at least stopfactor * Information * particleCount steps", 2.0);
	public Input<Integer> minStepsInput = new Input<>("minSteps", "minimal number of steps to take. Useful to guarantee chain does not get stuck in prior", 0);

	public Input<Boolean> autoSubChainLengthInput = new Input<>("autoSubChainLength", "automatically determines subchain length based on number of accepted steps (unless subChainLength * paramCount is reached)", true);
	public Input<Double> paramCountFactorInput = new Input<>("paramCountFactor", "determines length of subchain as multiplier of accepted steps before returning divided by number of parameters in the analysis"
			+ "ignored if autoSubChainLengt=false", 1.0);

	private static final boolean printDebugInfo = false;

	protected int particleCount;
	protected int subChainLength;

	protected String[] particleStates;
	protected double[] particleLikelihoods;
	protected double minLikelihood;

	protected Distribution likelihood, originalPrior = null;
	protected Distribution[] samplingDistribution;

	double newLogPrior, oldLogPrior;

	/** marginal likelihood estimate **/
	double Z;
	/** estimate of information **/
	double H;
	
	List<NSLogger> NSloggers;
	
	List<Double> likelihoods = new ArrayList<>();
	boolean finished = false;

	int paramCount = 0;
	double paramCountFactor = 1.0;
	boolean autoSubChainLength = true;
	
	public NS() {}
	
	public NS(int chainLength, int preBurnin, int particleCount, int subChainLength, State state, List<Operator> operators, CompoundDistribution distribution, Double epsilon) {
		initByName("chainLength", chainLength, 
				"preBurnin", preBurnin, 
				"particleCount", particleCount,
				"subChainLength", subChainLength,
				"state", state,
				"operator", operators,
				"distribution", distribution,
				"epsilon", epsilon);
	}
	
	@Override
	public void initAndValidate() {
		paramCountFactor = paramCountFactorInput.get();
		autoSubChainLength = autoSubChainLengthInput.get();

//		if (historyLengthInput.get() < 2) {
//			throw new IllegalArgumentException("history must be 2 or larger");
//		}
		// grab priors
		CompoundDistribution d = (CompoundDistribution) posteriorInput.get();
		List<Distribution> list = d.pDistributions.get();
		if (list.size() < 2) {
			throw new IllegalArgumentException("Expected one likelihood and at least one prior distribution.");
		}

		samplingDistribution = new Distribution[list.size() - 1];
		int nextPriorIndex = 0;
		for (Distribution pDist : list) {
			final String distID = pDist.getID().toLowerCase();
			if (distID.startsWith("likelihood")) {
				if (likelihood == null) {
					likelihood = pDist;
				} else {
					throw new IllegalArgumentException(
							"Expected only one likelihood distribution -- change the ID to something that does not start with likelihood and it will be deemed to be a prior.");
				}
			} else {
				samplingDistribution[nextPriorIndex] = pDist;
				nextPriorIndex++;
			}
		}

		if (samplingDistributionInput.get() != null) {
			// remove other priors
			List<Distribution> originalPriors = new ArrayList<>();
			for (int i =  list.size() - 1; i >= 0; i--) {
				final String distID = list.get(i).getID().toLowerCase();
				if (!distID.startsWith("likelihood")) {
					originalPriors.add(list.get(i));
					list.remove(i);
				}
			}
			if (originalPriors.size() == 1) {
				originalPrior = originalPriors.get(0);
			} else {
				originalPrior = new CompoundDistribution();
				originalPrior.initByName("distribution", originalPriors);
			}
			originalPrior.setID("originalPrior");
			list.add(originalPrior);
			// set sampling distribution from Input
			samplingDistribution = new Distribution[1];
			samplingDistribution[0] = samplingDistributionInput.get();
			// make sure samplingDistribution is calculated first
			if (list.contains(samplingDistribution[0])) {
				list.remove(samplingDistribution[0]);				
			}
			list.add(0, samplingDistribution[0]);
		}
	
		super.initAndValidate();

		particleCount = particleCountInput.get();
		if (particleCount <= 0) {
			throw new IllegalArgumentException("particle count must be a positive number");
		}

		subChainLength = subChainLengthInput.get();
		if (subChainLength <= 0) {
			throw new IllegalArgumentException("sub chain length must be a positive number");
		}

		particleStates = new String[particleCount];
		particleLikelihoods = new double[particleCount];


	
		Z = -Double.MAX_VALUE;
		H = 0;
		
		NSloggers = new ArrayList<>();
		
		
		paramCount = 0;
		for (StateNode node : state.stateNodeInput.get()) {
			if (node instanceof Parameter<?>) {
				Parameter<?> param = (Parameter<?>) node;
				paramCount += param.getDimension();
			} else if (node instanceof Tree) {
				Tree tree = (Tree) node;
				paramCount += tree.getNodeCount() * 2;
			}
		}
		Log.warning("Counting " + paramCount + " parameters");
	}

	@Override
	public void run() throws IOException, SAXException, ParserConfigurationException {
		// set up state (again). Other beastObjects may have manipulated the
		// StateNodes, e.g. set up bounds or dimensions
		state.initAndValidate();
		// also, initialise state with the file name to store and set-up whether
		// to resume from file
		state.setStateFileName(stateFileName);
		operatorSchedule.setStateFileName(stateFileName);

		burnIn = burnInInput.get();
		chainLength = chainLengthInput.get();
		int initialisationAttempts = 0;
		state.setEverythingDirty(true);
		posterior = posteriorInput.get();

		if (restoreFromFile) {
			state.restoreFromFile();
			operatorSchedule.restoreFromFile();
			burnIn = 0;
			oldLogPrior = state.robustlyCalcPosterior(posterior);
		} else {
			do {
				for (final StateNodeInitialiser initialiser : initialisersInput.get()) {
					initialiser.initStateNodes();
				}
				oldLogPrior = state.robustlyCalcPosterior(posterior);
				initialisationAttempts += 1;
			} while (Double.isInfinite(oldLogPrior) && initialisationAttempts < numInitializationAttempts.get());
		}
		final long startTime = System.currentTimeMillis();

		state.storeCalculationNodes();

		// do the sampling
		logAlpha = 0;
		debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

		// System.err.println("Start state:");
		// System.err.println(state.toString());

		Log.info.println("Start likelihood: " + oldLogPrior + " "
				+ (initialisationAttempts > 1 ? "after " + initialisationAttempts + " initialisation attempts" : ""));
		if (Double.isInfinite(oldLogPrior) || Double.isNaN(oldLogPrior)) {
			reportLogLikelihoods(posterior, "");
			throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.");
		}

		loggers = loggersInput.get();

		// put the loggers logging to stdout at the bottom of the logger list so
		// that screen output is tidier.
		Collections.sort(loggers, (o1, o2) -> {
			if (o1.isLoggingToStdout()) {
				return o2.isLoggingToStdout() ? 0 : 1;
			} else {
				return o2.isLoggingToStdout() ? -1 : 0;
			}
		});
		// warn if none of the loggers is to stdout, so no feedback is given on
		// screen
		boolean hasStdOutLogger = false;
		boolean hasScreenLog = false;
		for (Logger l : loggers) {
			if (l.isLoggingToStdout()) {
				hasStdOutLogger = true;
			}
			if (l.getID() != null && l.getID().equals("screenlog")) {
				hasScreenLog = true;
			}
		}
		if (!hasStdOutLogger) {
			Log.warning.println(
					"WARNING: If nothing seems to be happening on screen this is because none of the loggers give feedback to screen.");
			if (hasScreenLog) {
				Log.warning.println("WARNING: This happens when a filename  is specified for the 'screenlog' logger.");
				Log.warning.println("WARNING: To get feedback to screen, leave the filename for screenlog blank.");
				Log.warning.println("WARNING: Otherwise, the screenlog is saved into the specified file.");
			}
		}

		// replace trace loggers to file with NSloggers
		for (int i = loggers.size() - 1; i >= 0; i--) {
			Logger logger = loggers.get(i);
			if (logger instanceof NSLogger) {
				NSloggers.add((NSLogger) logger);
				loggers.remove(i);
			} else if (logger.mode == Logger.LOGMODE.compound && logger.fileNameInput.get() != null) {
				NSLogger nslogger = replaceLogger(logger);
				NSloggers.add(nslogger);
				loggers.remove(i);
			} else if (logger.mode == Logger.LOGMODE.autodetect && logger.fileNameInput.get() != null) {
				if (!(logger.loggersInput.get().size() == 1 && logger.loggersInput.get().get(0) instanceof TreeInterface)) {
					NSLogger nslogger = replaceLogger(logger);
					NSloggers.add(nslogger);
					loggers.remove(i);
				}
			} else {
				logger.initByName("logEvery", 1);
			}
		}

		// initialises log so that log file headers are written, etc.
		for (final Logger log : loggers) {
			log.init();
		}
		for (final Logger log : NSloggers) {
			log.init();
		}

		if (NSloggers.size() == 0) {
			Log.warning("");
			Log.warning("WARNING: no NSLoggers found -- no posterior sample will be generated");
			Log.warning("");
		}

		doLoop();

		Log.info.println();
			operatorSchedule.showOperatorRates(System.out);

		Log.info.println();
		final long endTime = System.currentTimeMillis();
		Log.info.println("Total calculation time: " + (endTime - startTime) / 1000.0 + " seconds");
		close();

		Log.warning.println("End likelihood: " + oldLogPrior);
		// System.err.println(state);
		state.storeToFile(chainLength);
		operatorSchedule.storeToFile();
		// Randomizer.storeToFile(stateFileName);
		
		for (final Logger log : loggers) {
			log.close();
		}
		for (final Logger log : NSloggers) {
			log.close();
		}

		// produce posterior samples
		Log.warning("Producing posterior samples");
		NSLogger nslogger0 = NSloggers.get(0);
		for (Logger logger : loggers) {
			if (logger.mode == Logger.LOGMODE.tree ||
				logger.mode == Logger.LOGMODE.autodetect && logger.fileNameInput.get() != null &&
				logger.loggersInput.get().get(0) instanceof TreeInterface) {
				
				NSLogAnalyser.main(new String[]{"-log", nslogger0.fileNameInput.get(),
						"-tree", logger.fileNameInput.get(),
						"-out", logger.fileNameInput.get() + ".posterior",
						"-N", particleCount+"",
						"-quiet"});
			}
		}
		for (NSLogger nslogger : NSloggers) {
			NSLogAnalyser.main(new String[]{"-log", nslogger.fileNameInput.get(),
					"-out", nslogger.fileNameInput.get() + ".posterior",
					"-N", particleCount+"",
					"-quiet"});
		}
		

	} // run;

	private NSLogger replaceLogger(Logger logger) {
		Log.warning("replacing logger " + logger.getID() + " with NSLogger");
		NSLogger nslogger = new NSLogger();
		nslogger.initByName("logEvery", logger.everyInput.get(),
				"fileName", logger.fileNameInput.get(),
				"model", logger.modelInput.get(),
				"sort", logger.sortModeInput.get(),
				"sanitiseHeaders", logger.sanitiseHeadersInput.get(),
				"log", logger.loggersInput.get()
				);
		nslogger.setID(logger.getID());
		return nslogger;
	}

	protected void initParticles() {
		minLikelihood = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < particleCount; i++) {
			for (int j = 0; j < subChainLength; j++) {
				composeProposal(1);
			}

			particleStates[i] = state.toXML(0);
			particleLikelihoods[i] = likelihood.getArrayValue();
			if (originalPrior != null) {
				particleLikelihoods[i] += originalPrior.getCurrentLogP() - samplingDistribution[0].getCurrentLogP(); 
			}
		}
		Log.warning(particleCount + " particles initialised");
	}

	/**
	 * main MCMC loop
	 * 
	 * @throws IOException
	 * *
	 */
	@Override
	protected void doLoop() throws IOException {
		double N = particleCount;	
		
		oldLogPrior = 0;
		for (Distribution d : samplingDistribution) {
			oldLogPrior += d.getArrayValue();
		}

		// initialise states
		initParticles();

		if (burnIn > 0) {
			Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
		}
		

		// sample from prior
		minLikelihood = Double.NEGATIVE_INFINITY; // guarantee that likelihood is ignored when sampling from prior
		for (int sampleNr = -burnIn; sampleNr < 0; sampleNr++) {
			// find particle with minimum likelihood
			int iMin = Randomizer.nextInt(particleCount);

			// init state
			state.fromXML(particleStates[Randomizer.nextInt(particleCount)]);

			// init calculation nodes & oldLogPrior
			robustlyCalcPosterior(posterior);
			oldLogPrior = 0;
			for (Distribution d : samplingDistribution) {
				oldLogPrior += d.getArrayValue();
			}
			
			for (int j = 0; j < subChainLength; j++) {
				composeProposal(j + subChainLength * sampleNr);
			}

			particleStates[iMin] = state.toXML(sampleNr);
			particleLikelihoods[iMin] = likelihood.getArrayValue();
			if (originalPrior != null) {
				particleLikelihoods[iMin] += originalPrior.getCurrentLogP() - samplingDistribution[0].getCurrentLogP(); 
			}

			callUserFunction(sampleNr);

			if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
				throw new RuntimeException(
						"Encountered a positive infinite posterior. This is a sign there may be improper priors or numeric instability in the model.");
			}
		}
		

		doInnerLoop(likelihoods);
		
		
 		
		double [] Zestimates = new double[SAMPLE_COUNT];
		for (int k = 0; k < SAMPLE_COUNT; k++) {
 	 		double logX = 0.0;
	 		double Z = 0;
			Z = -Double.MAX_VALUE;
	 		for (int i = 0; i < likelihoods.size(); i++) {
	 			double u = NSLogAnalyser.nextBeta(N, 1.0); 			
	 			double lw = logX + Math.log(1.0 - u);
	 			double L = lw  + likelihoods.get(i);
	 			Z = NS.logPlus(Z, L);
	 			logX += Math.log(u);
	 		}
	 		Zestimates[k] = Z;
 		}
		double sum = 0;
		for (double d : Zestimates) {
			sum += d;
		}
		Z = sum / SAMPLE_COUNT;
		double sum2 = 0;
		for (double d : Zestimates) {
			sum2 += (d - Z) * (d - Z);
		}
		double var = sum2/(SAMPLE_COUNT - 1.0);
		double stdev = Math.sqrt(var);
 		Log.info("Marginal likelihood: " + Z + "(" + stdev +")");

 		Log.info("Information: " + H);
 		Log.info("SD: " + Math.sqrt(H/N));
 	}

	/** The likelihoods List will be populated with minimum likelihoods for the set of particles 
	 * corrected for when using a sampling distribution. 
	 * @param likelihoods
	 * @throws IOException
	 */
	protected void doInnerLoop(List<Double> likelihoods) throws IOException {
		// run nested sampling
		double N = particleCount;	
		double lw;
		//double logW = Math.log(1.0 - Math.exp(-1.0/N));
		double logW = Math.log(1.0 - Math.exp(-2.0/N)) - Math.log(2.0) ;
		 
		double Zb = -Double.MAX_VALUE; // Zb is the marginal likelihood estimate from the previous iteration

		/** number of steps without change before stopping **/
		int HISTORY_LENGTH = 2;//historyLengthInput.get();
		/** size of required change for stopping criterion **/
		double EPSILON = epsilonInput.get();

		double stopFactor = stopFactorInput.get();
		int minSteps = minStepsInput.get();

		double [] mlHistory = new double[HISTORY_LENGTH];
		mlHistory[0] = -1.0; // to pass stop criterion when sampleNr = 0
		int sampleNr = 0;
		// continue while
		// o we have not reached a user specified upper bound of steps (through chainLength) AND
		//      o the number of samples is less than 2 * Information * #particles (stopFactor = 2 can be changed) OR
		//      o the relative gain in ML estimate is less than ESPILON
			while (sampleNr <= chainLength && (
					sampleNr < minSteps ||
					sampleNr < stopFactor * H * particleCount ||
					Math.abs(mlHistory[(sampleNr +HISTORY_LENGTH-1) % HISTORY_LENGTH] - mlHistory[sampleNr % HISTORY_LENGTH])/Math.abs(mlHistory[(sampleNr +HISTORY_LENGTH- 1) % HISTORY_LENGTH]) > EPSILON)) {

			// find particle with minimum likelihood
			int iMin = 0;
			minLikelihood = particleLikelihoods[0];
			for (int i = 1; i < particleCount; i++) {
				if (particleLikelihoods[i] < minLikelihood) {
					minLikelihood = particleLikelihoods[i];
					iMin = i;
				}
			}


			// Store Li + discarded points				
			if (particleCount > 1) {
				state.fromXML(particleStates[iMin]);
			}
			robustlyCalcPosterior(posterior);
			for (Logger logger: NSloggers) {
				if (logger instanceof NSLogger) {
					((NSLogger) logger).log(sampleNr, minLikelihood);
				}
			}
			log(sampleNr);

			// mean(theta) = \sum_i theta * Li*wi/Z
			// Z = \sum_i Li*wi
			
			updateParticle(sampleNr);
			
			if (sampleNr > 0) {
//				double Xi = Math.exp(-sampleNr/N);
//				double Xi_1 = Math.exp(-(sampleNr-1.0)/N);
//				double wi = Xi_1 - Xi;
				//delta = Li * wi; 
				//Z  += delta;

				lw = logW - (sampleNr - 1.0) / N;
				double Li = minLikelihood;
//				if (originalPrior != null) {
//					Li += originalPrior.getCurrentLogP() - samplingDistribution[0].getCurrentLogP(); 
//				}
				likelihoods.add(Li);
				
				double L = lw  + Li;
				
				Z = logPlus(Z, L);
				H = Math.exp(L - Z) * Li - Z + Math.exp(Zb - Z)*(H + Zb);
				Zb = Z;
				if (sampleNr % 1 == 0) {
			 		Log.info("ML: " + Z + " Information: " + H);
				}
				mlHistory[sampleNr % HISTORY_LENGTH] = Z;
			}
			
			double c = 0;
			if (originalPrior != null) {
				c = originalPrior.getCurrentLogP() - samplingDistribution[0].getCurrentLogP(); 
			}
			updateParticleState(iMin, state.toXML(sampleNr), likelihood.getArrayValue() + c);

			callUserFunction(sampleNr);

			// make sure we always save just before exiting
			if (storeEvery > 0 && (sampleNr + 1) % storeEvery == 0 || sampleNr == chainLength) {
				/* final double logLikelihood = */
				state.robustlyCalcNonStochasticPosterior(posterior);
				state.storeToFile(sampleNr);
				operatorSchedule.storeToFile();
			}

			if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
				throw new RuntimeException(
						"Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
			}
			sampleNr++;
		}
			
		Log.info(sampleNr+"<="+ chainLength +"&& ("+ sampleNr +"<"+ stopFactor +"*"+ H +"*"+ particleCount +"||"
					+"Math.abs("+mlHistory[(sampleNr +HISTORY_LENGTH-1) % HISTORY_LENGTH]+" - "+mlHistory[sampleNr % HISTORY_LENGTH]+")/Math.abs("+mlHistory[(sampleNr +HISTORY_LENGTH- 1) % HISTORY_LENGTH]+") > "+EPSILON);

			
		Log.info("Finished in " + sampleNr + " steps!");
		
		if (System.getProperty("beast.debug") != null) { 
			Log.warning(Arrays.toString(mlHistory));
			Log.warning("("+mlHistory[sampleNr % HISTORY_LENGTH] +"-"+ mlHistory[(sampleNr +1 )% HISTORY_LENGTH]+")/"+mlHistory[sampleNr % HISTORY_LENGTH] +">"+ EPSILON);
			
			double m = likelihoods.get(likelihoods.size()-1);
	 		double Z2 = 0;
	 		for (int i = 0; i < likelihoods.size(); i++) {
	 			double Xi = Math.exp(-i/N);
	 			double Xi_1 = Math.exp(-(i-1.0)/N);
	 			double wi = Xi_1 - Xi;
	 			Z2 += wi * Math.exp(likelihoods.get(i) - m);
	 		}
	 		Z2 = Math.log(Z2) + m;
	 		
	 		Log.info("Marginal likelihood: " + Z + " " + Z2);
		}
		finished = true;
	}

	protected void updateParticleState(int iMin, String state, double likelihood) {
		particleStates[iMin] = state;
		particleLikelihoods[iMin] = likelihood;
	}

	// propose a new starting state for particle such that likelihood > minLikelihood
	protected void updateParticle(int sampleNr) {
		if (particleCount > 1) {
			// init state
			state.fromXML(particleStates[Randomizer.nextInt(particleCount)]);

			// init calculation nodes & oldLogPrior
			robustlyCalcPosterior(posterior);

		}
		oldLogPrior = 0;
		for (Distribution d : samplingDistribution) {
			oldLogPrior += d.getArrayValue();
		}
		
		
		if (autoSubChainLength) {
			int acceptCount = 0;
			int j = 0;
			while (acceptCount < paramCount * paramCountFactor && j < paramCount * subChainLength) {
				boolean accept = composeProposal(j + subChainLength * sampleNr);
				if (accept) {
					acceptCount++;
				}
				j++;
			}
			System.err.println("autoSubChainLength = " + j);
		} else {
			for (int j = 0; j < subChainLength; j++) {
				composeProposal(j + subChainLength * sampleNr);
			}
		}
	}

	/** return marginal likelihood estimate **/
	public double getMarginalLikelihood() {
		return Z;
	}

	/** return standard deviation of marginal likelihood estimate **/
	public double getStandardDeviation() {
		return Math.sqrt(H/particleCount);
	}
	
	/** return estimate of information in the model **/
	public double getInformation() {
		return H;
	}
	
	
	public static double logPlus(double x, double y) {
		return x > y ? x + Math.log(1.0 + Math.exp(y-x)) :  y + Math.log(1.0 + Math.exp(x-y));
	}
	
	

	/**
	 * The MCMC algorithm to compose proposal distributions and the operators.
	 * 
	 * @param currState
	 *            the current state
	 * @return the selected {@link beast.core.Operator}
	 */
	protected boolean composeProposal(final int currState) {
		boolean accept = false;
		state.store(currState);

		final Operator operator = operatorSchedule.selectOperator();
		if (printDebugInfo)
			System.err.print("\n" + currState + " " + operator.getName() + ":");

		final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
		Evaluator evaluator = null;

		if (evaluatorDistribution != null) {
			evaluator = new Evaluator() {
				@Override
				public double evaluate() {
					double logP = 0.0;

					state.storeCalculationNodes();
					state.checkCalculationNodesDirtiness();

					try {
						logP = evaluatorDistribution.calculateLogP();
					} catch (Exception e) {
						e.printStackTrace();
						System.exit(1);
					}

					state.restore();
					state.store(currState);

					return logP;
				}
			};
		}
		final double logHastingsRatio = operator.proposal(evaluator);

		if (logHastingsRatio != Double.NEGATIVE_INFINITY) {

			if (operator.requiresStateInitialisation()) {
				state.storeCalculationNodes();
				state.checkCalculationNodesDirtiness();
			}

			newLogPrior = posterior.calculateLogP();

			newLogPrior = 0;
			for (Distribution d : samplingDistribution) {
				newLogPrior += d.getArrayValue();
			}

			logAlpha = newLogPrior - oldLogPrior + logHastingsRatio; // CHECK
																		// HASTINGS
			double c = 0;
			if (originalPrior != null) {
				c = originalPrior.getCurrentLogP() - newLogPrior; 
			}
			if (printDebugInfo)
				System.err.print(logAlpha + " " + newLogPrior + " " + oldLogPrior);

			if ((logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha))
					&& likelihood.getArrayValue() + c > minLikelihood) {
				// accept
				oldLogPrior = newLogPrior;
				state.acceptCalculationNodes();

				if (currState >= 0) {
					operator.accept();
				}
				if (printDebugInfo)
					System.err.print(" accept");
				accept = true;
			} else {
				// reject
				if (currState >= 0) {
					operator.reject(newLogPrior == Double.NEGATIVE_INFINITY ? -1 : 0);
				}
				state.restore();
				state.restoreCalculationNodes();
				if (printDebugInfo)
					System.err.print(" reject");
			}
			state.setEverythingDirty(false);
		} else {
			// operation failed
			if (currState >= 0) {
				operator.reject(-2);
			}
			state.restore();
			if (!operator.requiresStateInitialisation()) {
				state.setEverythingDirty(false);
				state.restoreCalculationNodes();
			}
			if (printDebugInfo)
				System.err.print(" direct reject");
		}
		return accept;
	}
}
