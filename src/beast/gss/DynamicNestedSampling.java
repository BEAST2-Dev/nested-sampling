package beast.gss;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.NSLogger;
import beast.core.Operator;
import beast.core.State;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.tree.TreeInterface;
import beast.util.NSLogAnalyser;



@Description("Dynamic nested sampling for phylogenetics")
public class DynamicNestedSampling extends NS {
	enum Goal {evidence, posterior, mixture}
	public Input<Goal> goalInput = new Input<>("goal", "fraction of importance to be revisited", Goal.posterior, Goal.values());
	
	public Input<Double> gInput = new Input<>("G", "G-parameter from DNS paper. Ignored unless goal=\"mixture\". "
			+ "Determines whether more particles are drawn to get accurate "
			+ "evidence/marginal likelihood estimate (goal=0.0) or more accurate posterior sample (goal=1.0)", 1.0);
	public Input<Double> fractionInput = new Input<>("fraction", "fraction of importance to be revisited", 0.99);
	public Input<Double> targetSDInput = new Input<>("targetSD", "target standard deviation for the marginal likelihood/evidence", 3.0);
	
	
	private Goal goal;
	private double g, fraction;
	private double [] importance;
	private int [] n0;
	private double [] L0;
	private String [] states0;
	private double [] weights0;
	private int [] order0;
	private double evidenceSD, targetSD;


	public DynamicNestedSampling() {
	}
	
	public DynamicNestedSampling(int chainLength, int preBurnin, int particleCount, int subChainLength, State state, List<Operator> operators, CompoundDistribution distribution, Double epsilon) {
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
		super.initAndValidate();
		
		if (fractionInput.get() <= 0 || fractionInput.get() >= 1) {
			throw new IllegalArgumentException("fraction should be between 0 and 1, not " + fractionInput.get());
		}
		fraction = fractionInput.get();
		
		goal = goalInput.get();
		switch (goal) {
		case mixture:
			if (gInput.get() < 0 || gInput.get() > 1) {
				throw new IllegalArgumentException("goal should be between 0 and 1, not " + gInput.get());
			}
			g = gInput.get();
			break;
		case posterior:
			g = 1;
			break;
		case evidence:
			g = 0;
			break;			
		}
		
		states = new ArrayList<>();
		
		targetSD = targetSDInput.get();
		if (targetSD <= 0) {
			throw new IllegalArgumentException("targetSD should be positive, not " + targetSDInput.get());
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
		state.setEverythingDirty(true);
		posterior = posteriorInput.get();

		if (restoreFromFile) {
			throw new RuntimeException("Resume is not implemented yet for dynamic nested sampling");
//			state.restoreFromFile();
//			operatorSchedule.restoreFromFile();
//			burnIn = 0;
//			oldLogPrior = state.robustlyCalcPosterior(posterior);
		}
		final long startTime = System.currentTimeMillis();

		state.storeCalculationNodes();

		// do the sampling
		logAlpha = 0;
		debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

		// System.err.println("Start state:");
		// System.err.println(state.toString());

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


		
		
		doDynamicLoop();
		


		
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

 		// max nr of posterior samples
 		double ESS = 0;
 		for (int i = 0; i < weights0.length; i++) {
 			if (weights0[i] > 0) {
 				ESS -= weights0[i]* Math.log(weights0[i]);
 			}
 		}
 		ESS = Math.exp(ESS);
 		
 		Log.warning("Max ESS: " + ESS);

		// produce posterior samples
		Log.warning("Producing posterior samples");
		NSLogger nslogger0 = NSloggers.get(0);
		for (Logger logger : loggers) {
			if (logger.mode == Logger.LOGMODE.tree ||
				logger.mode == Logger.LOGMODE.autodetect && logger.fileNameInput.get() != null &&
				logger.loggersInput.get().get(0) instanceof TreeInterface) {
				
				NSLogAnalyser.resampleToFile(true, nslogger0.fileNameInput.get(), ESS, weights0, order0);
			}
		}
		for (NSLogger nslogger : NSloggers) {
			NSLogAnalyser.resampleToFile(false, nslogger0.fileNameInput.get(), ESS, weights0, order0);
		}
		

	} // run;




// Algorithm 1 from Dynamic Nested Sampling paper https://arxiv.org/pdf/1704.03459.pdf		
//	Generate a nested sampling run with a constant number of live points n_init;
//	while (dynamic termination condition not satisfied) do
//		recalculate importance I(G,i) of all points;
//		find first point j and last point k with importance of greater than some fraction f (we use f=0.9) of the largest importance;
//		generate a additional thread (or alternatively n batch	additional threads) starting at L_{j−1} 
//		and ending with the first sample taken with likelihood greater than L_{k+1};
//	end
	private void doDynamicLoop() throws IOException {
//		Generate a nested sampling run with a constant number of live points n_init;
		doLoop();
		updateLikelihoodsAndN();

		restoreFromFile = false;

		while (true) {
//			recalculate importance I(G,i) of all points;
			calcImportancePerPoint();

//			if (dynamic termination condition not satisfied) return
			if (evidenceSD < targetSD) {
				return;
			}
			
//			find first point j and last point k with importance of greater than some fraction f (we use f=0.9) of the largest importance;
			int [] startend = new int[2];
			findFirstLast(startend);
			int j = startend[0];
			int k = startend[1];
			
//			generate a additional thread (or alternatively n batch	additional threads) starting at L_{j−1}
			if (j == 0) {
				initParticles(null, Double.NEGATIVE_INFINITY);
			} else {
				initParticles(states0[j-1], L0[j-1]);
			}
			maxLikelihood = k+1 < L0.length ? L0[k + 1] : Double.POSITIVE_INFINITY;
			
//			and ending with the first sample taken with likelihood greater than L_{k+1};
			doInnerLoop();
			updateLikelihoodsAndN();
		}		
//		end
	}


	// merge in likelihoods & states and update n accordingly
	private void updateLikelihoodsAndN() {
		
		// initialise from start run
		if (n0 == null) {
			n0 = new int[likelihoods.size()];
			L0 = new double[likelihoods.size()];
			states0 = new String[likelihoods.size()];
			order0 = new int[likelihoods.size()];
			for (int i = 0; i < likelihoods.size(); i++) {
				n0[i] = particleCount; 
				L0[i] = likelihoods.get(i);
				states0[i] = states.get(i);
				order0[i] = i;
			}
			
			likelihoods.clear();
			states.clear();
			return;
		}
		int size = n0.length;

		int [] n2 = new int[size + likelihoods.size()];
		System.arraycopy(n0, 0, n2, 0, size);
		n0 = n2;
		double [] L2 = new double[size + likelihoods.size()];
		System.arraycopy(L0, 0, L2, 0, size);
		L0 = L2;
		String [] S2 = new String[size + likelihoods.size()];
		System.arraycopy(states0, 0, S2, 0, size);
		states0 = S2;		
		n2 = new int[size + likelihoods.size()];
		System.arraycopy(order0, 0, n2, 0, size);
		order0 = n2;

		
		// merge follow up runs
		int start = 0;
					
		int iStart = -1, iEnd = -1;
		for (int k = 0; k < likelihoods.size(); k++) {
			double likelihood = likelihoods.get(k);
			int i = Arrays.binarySearch(L0, 0, size, likelihood);
			if (i < 0) {
				// (-(insertion point) - 1
				i = 1 - i;
			}
			System.arraycopy(L0, i, L0, i+1, size - i);
			L0[i] = likelihood;
			System.arraycopy(states0, i, states0, i+1, size - i);
			states0[i] = states.get(k);
			System.arraycopy(n0, i, n0, i+1, size - i);
			System.arraycopy(order0, i, order0, i+1, size - i);
			order0[i] = size;
			size++;
			
			if (k == start) {
				iStart = i;
			} else {
				iEnd = i;
			}
		}
		for (int k = iStart; k < iEnd; k++) {
			n0[k] = n0[k] + particleCount;
		}
		
		likelihoods.clear();
		states.clear();
	}
	
	
	private void calcImportancePerPoint() {
		double Z = -Double.MAX_VALUE; 

		// evidence importance
		double [] IZ = new double[n0.length];
		// parameter importance
		weights0 = new double[n0.length];
		
//		double logX = 0;
//		double logX1 = - Math.log(n0[0]);
//		double logX2 = logX1 - Math.log(n0[1]);
//		double log2 = Math.log(2.0);
//				
//		for (int sampleNr = 0; sampleNr < n0.length-2; sampleNr++) {
//			
//			double logW = Math.log(1.0 - Math.exp(-2.0/n0[sampleNr])) - Math.log(2.0) ;
//			double lw = logW - (sampleNr - 1.0) / n0[sampleNr];
//
////			double lw = -log2 + logX + Math.log(1.0-Math.exp(-n[sampleNr]-n[sampleNr + 1]));
//			logX = logX1;
//			logX1 = logX2;
//			logX2 += 1.0/n0[sampleNr+2];
//			
//			double Li = L0[sampleNr];						
//			double L0 = lw  + Li;
//			
//			IP[sampleNr] = L0/n0[sampleNr];
//			Z = logPlus(Z, L0);
//			IZ[sampleNr] = Z;
//		}
		
		
		
		
		
 		double zMean = 0;
 		double v = 0;
 		double hMean = 0, H = 0;
 		final double RESAMPLE_COUNT = 100;
 		Arrays.fill(IZ,  0);
 		Arrays.fill(weights0,  0);
 		
 		for (int k = 0; k < RESAMPLE_COUNT; k++) {
 			double logX = 0;
			H = 0;
		 	double Zb = -Double.MAX_VALUE;
	 		Z = -Double.MAX_VALUE;
	 		for (int i = 0; i < n0.length; i++) {
	 			double u = NSLogAnalyser.nextBeta(n0[i], 1.0); 			
	 			double lw = logX + Math.log(1.0 - u);
	 			double L = lw  + L0[i];
	 			weights0[i] += L;
	 			Z = NS.logPlus(Z, L);
	 			IZ[i] += Z;
	 			// weights[i] += L; 			
	 			logX += Math.log(u);
	 			if (i > 0) {
	 				double Li = L0[i-1];
	 				H = Math.exp(L - Z) * Li - Z + Math.exp(Zb - Z)*(H + Zb);
	 				// Log.info("Math.exp("+L+" - "+Z+") * "+Li+ " -"+ Z +"+ Math.exp("+Zb+"- "+ Z+")*("+H+" + "+Zb+")");
	 				Zb = Z;
	 			}
	 		}
			
	 		zMean += Z;
			v += Z * Z;
	 		hMean += H;
		}
		for (int i = 0; i < n0.length; i++) {
			weights0[i] /= RESAMPLE_COUNT;
		}
		for (int i = 0; i < n0.length; i++) {
			IZ[i] /= RESAMPLE_COUNT;
		}
		
		Z = zMean / RESAMPLE_COUNT;
		v = Math.sqrt(v/RESAMPLE_COUNT-Z*Z);
		H = hMean / RESAMPLE_COUNT;
		Log.warning("\nMarginal likelihood: " + Z + " SD=(" + v + ") Information: " + H);
		evidenceSD = v;

		
		
		
		// IZ[n] now contains \sum_i<n log(Z_i)
		// we want IZ[n] = \sum_i>n Z
		IZ[n0.length-2] = IZ[n0.length-3]; 
		IZ[n0.length-1] = IZ[n0.length-2]; 
		double max = IZ[IZ.length - 1];
		for (int i = 0; i < n0.length; i++) {
			IZ[i] -= max;
			IZ[i] = 1.0 - Math.exp(IZ[i]);
		}

		// IP[n] now contains log(L_iw_i)
		max = weights0[0];
		for (double d : weights0) {
			if (d > max) {
				max = d;
			}
		}
		for (int i = 0; i < n0.length; i++) {
			weights0[i] = Math.exp(weights0[i] - max);
		}
		
		double sumZ = 0;
		for (double d : IZ) {
			sumZ += d;
		}

		double sumP = 0;
		for (double d : weights0) {
			sumP += d;
		}

		
		// importance per point
		importance = new double[n0.length];
		for (int i = 0; i < n0.length; i++) {
			importance[i] = (1-g) * IZ[i]/sumZ + g * weights0[i]/sumP;
		}		
	}

	// find first point j and last point k with importance of greater than 
	// some fraction f (we use f=0.9) of the largest importance;
	private void findFirstLast(int[] startend) {
		if (goal.equals(Goal.posterior)) {
			// determine HPD over `fraction` as start & end points
			int start = 0;
			int end = importance.length - 1;
			double sum = 0;
			while (sum < 1 - fraction) {
				if (importance[start]*n0[start] <= importance[end]*n0[end]) {
					sum += importance[start];
					start++;
				} else {
					sum += importance[end];
					end--;
				}
			}
			startend[0] = start;
			startend[1] = end;
			return;
		}
		
		double max = 0;
		for (double d : importance) {
			if (d > max) {
				max = d;
			}
		}
		
		double maxFraction = max * fraction;
		startend[0] = -1;
		startend[1] = -1;
		for (int i = 0; i < importance.length; i++) {
			if (importance[i] > maxFraction) {
				if (startend[0] == -1) {
					startend[0]  = i;
					startend[1] = i;
				} else {
					startend[1] = i;
				}
			}
		}
	}

}
