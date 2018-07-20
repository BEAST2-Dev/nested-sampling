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

	public Input<Double> goalInput = new Input<>("goal", "determines whetehr more particles are drawn to get accurate "
			+ "evidence/marginal likelihood estimate (goal=0.0) or more accurate posterior sample (goal=1.0)", 0.0);
	public Input<Double> fractionInput = new Input<>("fraction", "fraction of importance to be revisited", 0.9);
	
	private double goal, fraction;
	private double [] importance;
	
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
		if (goalInput.get() < 0 || goalInput.get() > 1) {
			throw new IllegalArgumentException("goal should be between 0 and 1, not " + goalInput.get());
		}
		goal = goalInput.get();
		
		states = new ArrayList<>();
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

//		while (dynamic termination condition not satisfied) do
		while (true) {
//			recalculate importance I(G,i) of all points;
			calcImportancePerPoint();
			
//			find first point j and last point k with importance of greater than some fraction f (we use f=0.9) of the largest importance;
			int [] startend = new int[2];
			findFirstLast(startend);
			int j = startend[0];
			int k = startend[1];
			
//			generate a additional thread (or alternatively n batch	additional threads) starting at L_{j−1}
			if (j == 0) {
				initParticles(null, Double.NEGATIVE_INFINITY);
			} else {
				initParticles(states0[j-1], L[j-1]);
			}
			maxLikelihood = k+1 < L.length ? L[k + 1] : Double.POSITIVE_INFINITY;
			
//			and ending with the first sample taken with likelihood greater than L_{k+1};
			doInnerLoop();
			updateLikelihoodsAndN();
		}		
//		end
	}


	int [] n;
	double [] L;
	String [] states0;
	
	// merge in likelihoods & states and update n accordingly
	private void updateLikelihoodsAndN() {
		
		// initialise from start run
		if (n == null) {
			n = new int[likelihoods.size()];
			L = new double[likelihoods.size()];
			states0 = new String[likelihoods.size()];
			for (int i = 0; i < likelihoods.size(); i++) {
				n[i] = particleCount; 
				L[i] = likelihoods.get(i);
				states0[i] = states.get(i);
			}
			
			likelihoods.clear();
			states.clear();
			return;
		}
		int size = n.length;

		int [] n2 = new int[n.length + size];
		System.arraycopy(n, 0, n2, 0, n.length);
		n = n2;
		double [] L2 = new double[L.length + size];
		System.arraycopy(L, 0, L2, 0, L.length);
		L = L2;
		String [] S2 = new String[states0.length + size];
		System.arraycopy(states0, 0, S2, 0, states0.length);
		states0 = S2;
		
		// merge follow up runs
		int start = 0;
					
		int iStart = -1, iEnd = -1;
		for (int k = 0; k < likelihoods.size(); k++) {
			double likelihood = likelihoods.get(k);
			int i = Arrays.binarySearch(L, 0, size + k, likelihood);
			if (i < 0) {
				// (-(insertion point) - 1
				i = 1 - i;
			}
			System.arraycopy(L, i, L, i+1, size - i);
			L[i] = likelihood;
			System.arraycopy(states0, i, states0, i+1, size - i);
			states0[i] = states.get(k);
			System.arraycopy(n, i, n, i+1, size - i);
			size++;
			
			if (k == start) {
				iStart = i;
			} else {
				iEnd = i;
			}
		}
		for (int k = iStart; k < iEnd; k++) {
			n[k] = n[k] + particleCount;
		}
		
		likelihoods.clear();
		states.clear();
	}
	
	
	private void calcImportancePerPoint() {
		double Z = -Double.MAX_VALUE; 

		// evidence importance
		double [] IZ = new double[n.length];
		// parameter importance
		double [] IP = new double[n.length];
		
		double logX = 0;
		double logX1 = - Math.log(n[0]);
		double logX2 = logX1 - Math.log(n[1]);
		double log2 = Math.log(2.0);
				
		for (int sampleNr = 0; sampleNr < n.length-2; sampleNr++) {
			
			double logW = Math.log(1.0 - Math.exp(-2.0/n[sampleNr])) - Math.log(2.0) ;
			double lw = logW - (sampleNr - 1.0) / n[sampleNr];

//			double lw = -log2 + logX + Math.log(1.0-Math.exp(-n[sampleNr]-n[sampleNr + 1]));
			logX = logX1;
			logX1 = logX2;
			logX2 += 1.0/n[sampleNr+2];
			
			double Li = L[sampleNr];						
			double L0 = lw  + Li;
			
			IP[sampleNr] = L0/n[sampleNr];
			Z = logPlus(Z, L0);
			IZ[sampleNr] = Z;
		}
		
		// IZ[n] now contains \sum_i<n log(Z_i)
		// we want IZ[n] = \sum_i>n Z
		IZ[n.length-2] = IZ[n.length-3]; 
		IZ[n.length-1] = IZ[n.length-2]; 
		double max = IZ[IZ.length - 1];
		for (int i = 0; i < n.length; i++) {
			IZ[i] -= max;
			IZ[i] = 1.0 - Math.exp(IZ[i]);
		}

		// IP[n] now contains log(L_iw_i)
		max = IP[0];
		for (double d : IP) {
			if (d > max) {
				max = d;
			}
		}
		for (int i = 0; i < n.length; i++) {
			IP[i] = Math.exp(IP[i] - max);
		}
		
		double sumZ = 0;
		for (double d : IZ) {
			sumZ += d;
		}

		double sumP = 0;
		for (double d : IP) {
			sumP += d;
		}

		
		// importance per point
		importance = new double[n.length];
		for (int i = 0; i < n.length; i++) {
			importance[i] = (1-goal) * IZ[i]/sumZ + goal * IP[i]/sumP;
		}		
	}

	// find first point j and last point k with importance of greater than 
	// some fraction f (we use f=0.9) of the largest importance;
	private void findFirstLast(int[] startend) {
		
		double max = 0;
		for (double d : importance) {
			if (d > max) {
				max = d;
			}
		}
		
		double maxFraction = max * fraction;
		
		for (int i = 0; i < importance.length; i++) {
			if (importance[i] > maxFraction) {
				if (startend[0] == -1) {
					startend[0]  = i;
				} else {
					startend[1] = i;
				}
			}
		}
	}

}
