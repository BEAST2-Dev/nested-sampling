package beast.gss;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.StateNodeInitialiser;
import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.core.util.Log;
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
public class NS extends MCMC {
	public Input<Integer> particleCountInput = new Input<>("particleCount", "number of particles (default 100)", 100);
	public Input<Integer> subChainLengthInput = new Input<>("subChainLength",
			"number of MCMC samples for each epoch (default 1000)", 1000);

	final public Input<Distribution> samplingDistributionInput = new Input<>("samplingDistribution",
			"probability distribution to sample from (e.g. a prior). "
					+ "If not specified, everything but the likelihood will be used as sampling distribution.");

	private static final boolean printDebugInfo = false;

	protected int particleCount;
	protected int subChainLength;

	protected String[] particleStates;
	protected double[] particleLikelihoods;
	protected double minLikelihood;

	protected Distribution likelihood;
	protected Distribution[] samplingDistribution;

	@Override
	public void initAndValidate() {
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
			samplingDistribution = new Distribution[1];
			samplingDistribution[0] = samplingDistributionInput.get();
			boolean alreadyAdded = false;
			for (Distribution d2 : d.pDistributions.get()) {
				if (d2 == samplingDistribution[0]) {
					alreadyAdded = true;
					break;
				}
			}
			if (!alreadyAdded) {
				d.pDistributions.setValue(samplingDistribution[0], d);
			}
		}
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

		// initialises log so that log file headers are written, etc.
		for (final Logger log : loggers) {
			log.init();
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
	} // run;

	private void initParticles() {
		minLikelihood = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < particleCount; i++) {
			for (int j = 0; j < subChainLength; j++) {
				composeProposal(1);
			}

			particleStates[i] = state.toXML(0);
			particleLikelihoods[i] = likelihood.getArrayValue();
		}
		Log.warning("particles initialised");
	}

	/**
	 * main MCMC loop
	 * 
	 * @throws IOException
	 * *
	 */
	@Override
	protected void doLoop() throws IOException {
		double Z = -Double.MAX_VALUE;
		double delta = Double.POSITIVE_INFINITY; 
		List<Double> likelihoods = new ArrayList<>();
		double N = particleCount;	
		double lw;
		
		oldLogPrior = 0;
		for (Distribution d : samplingDistribution) {
			oldLogPrior += d.getArrayValue();
		}

		// initialise states
		initParticles();

		if (burnIn > 0) {
			Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
		}
		
		double logW = Math.log(1.0 - Math.exp(-1.0/N));

		double H = 0;
		double Zb = -Double.MAX_VALUE; // Zb is the marginal likelihood estimate from the previous iteration

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

			callUserFunction(sampleNr);

			if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
				throw new RuntimeException(
						"Encountered a positive infinite posterior. This is a sign there may be improper priors or numeric instability in the model.");
			}
		}
		
		// run nested sampling
		for (int sampleNr = 0; sampleNr <= chainLength; sampleNr++) {

			// find particle with minimum likelihood
			int iMin = 0;
			minLikelihood = particleLikelihoods[0];
			for (int i = 1; i < particleCount; i++) {
				if (particleLikelihoods[i] < minLikelihood) {
					minLikelihood = particleLikelihoods[i];
					iMin = i;
				}
			}

			likelihoods.add(minLikelihood);
			// init state
			state.fromXML(particleStates[iMin]);
			// RRB: use random selected state instead of the worst one?
			// if so, replace above line with the next line
			//state.fromXML(particleStates[Randomizer.nextInt(particleCount)]);

			// init calculation nodes & oldLogPrior
			robustlyCalcPosterior(posterior);
			oldLogPrior = 0;
			for (Distribution d : samplingDistribution) {
				oldLogPrior += d.getArrayValue();
			}
			
			for (int j = 0; j < subChainLength; j++) {
				composeProposal(j + subChainLength * sampleNr);
			}

			if (sampleNr > 0) {
//				double Xi = Math.exp(-sampleNr/N);
//				double Xi_1 = Math.exp(-(sampleNr-1.0)/N);
//				double wi = Xi_1 - Xi;
				//delta = Li * wi; 
				//Z  += delta;

				lw = logW - (sampleNr - 1.0) / N;
				double Li = minLikelihood;
				double L = lw  + Li;
				
				Z = logPlus(Z, L);
				H = Math.exp(L - Z) * Li - Z + Math.exp(Zb - Z)*(H + Zb);
				Zb = Z;
				if (sampleNr % 100 == 0) {
			 		Log.info("Marginal likelihood: " + Z);
			 		Log.info("Information: " + H);
				}
			}

			particleStates[iMin] = state.toXML(sampleNr);
			particleLikelihoods[iMin] = likelihood.getArrayValue();

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
		}
		
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
 		Log.info("Information: " + H);
 		Log.info("SD: " + Math.sqrt(H/N));
 	}

	double logPlus(double x, double y) {
		return x > y ? x + Math.log(1.0 + Math.exp(y-x)) :  y + Math.log(1.0 + Math.exp(x-y));
	}
	
	
	double newLogPrior, oldLogPrior;

	/**
	 * The MCMC algorithm to compose proposal distributions and the operators.
	 * 
	 * @param currState
	 *            the current state
	 * @return the selected {@link beast.core.Operator}
	 */
	protected Operator composeProposal(final int currState) {
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
			if (printDebugInfo)
				System.err.print(logAlpha + " " + newLogPrior + " " + oldLogPrior);

			if ((logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha))
					&& likelihood.getArrayValue() > minLikelihood) {
				// accept
				oldLogPrior = newLogPrior;
				state.acceptCalculationNodes();

				if (currState >= 0) {
					operator.accept();
				}
				if (printDebugInfo)
					System.err.print(" accept");
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
		log(currState);
		return operator;
	}
}
