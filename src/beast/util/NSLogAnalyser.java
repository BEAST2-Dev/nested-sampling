package beast.util;

import static beast.util.OutputUtils.format;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.app.BEASTVersion2;
import beast.app.util.Utils;
import beast.core.util.Log;
import beast.gss.NS;
import beast.util.LogAnalyser;

public class NSLogAnalyser extends LogAnalyser {
	public int particleCount = -1;
	
    public NSLogAnalyser(String absolutePath, int burninPercentage, boolean quiet, int particleCount) throws IOException {
		super(absolutePath, burninPercentage, quiet, false);
		this.particleCount = particleCount;
	}


	public NSLogAnalyser(NSLogAnalyser analyser, NSLogAnalyser analyser2) {
		// make sure the same labels are used
		m_sLabels = analyser.m_sLabels;
		m_types = analyser.m_types;
		String [] m_sLabels2 = analyser2.m_sLabels;
		if (m_sLabels.length != m_sLabels2.length) {
			throw new IllegalArgumentException("Incompatible log files: different number of columns");
		}
		for (int i = 0; i < m_sLabels.length; i++) {
			if (!m_sLabels[i].equals(m_sLabels2[i])) {
				throw new IllegalArgumentException("Incompatible log files: column " + m_sLabels[i] + " != " + m_sLabels2[i]);
			}
		}
		
		particleCount = analyser.particleCount + analyser2.particleCount;
		
		int n = m_sLabels.length;
		
		int m1 = analyser.m_fTraces[0].length;
		int m2 = analyser2.m_fTraces[0].length;
		int m = m1 + m2; 
		m_fTraces = new Double[n][m];

		int i = 0;
		int i1 = 0;
		int i2 = 0;
		while (i < m) {
			if (i1 != m1 && (i2 == m2 || analyser.m_fTraces[1][i1] < analyser2.m_fTraces[1][i2])) {
				for (int j = 0; j < n; j++) {
					m_fTraces[j][i] = analyser.m_fTraces[j][i1];
				}
				i1++;
			} else {
				for (int j = 0; j < n; j++) {
					m_fTraces[j][i] = analyser2.m_fTraces[j][i2];
				}
				i2++;
			}
			i++;
		}
		
	}

	

	/**
     * calculate statistics on the data, one per column.
     * First column (sample nr) is not set *
     */
    @Override
    public void calcStats() {
    	if (!m_sLabels[1].equals("NSlikelihood")) {
    		throw new IllegalArgumentException("This does not appear to be generated by NS sampling since the second column is not 'NSlikelihood'");
    	}
    	// calc marginal likelihood Z
    	Double [] NSLikelihoods = m_fTraces[1];
		// Z = \sum_i Li*wi
 		double Z = 0;
		Z = -Double.MAX_VALUE;
 		int N = particleCount;
 		double [] weights = new double[NSLikelihoods.length];
 		double logW = Math.log(1.0 - Math.exp(-1.0/N));
 		for (int i = 0; i < NSLikelihoods.length; i++) {
 			double lw = logW - (i - 1.0) / N;
 			double L = lw  + NSLikelihoods[i];
 			Z = NS.logPlus(Z, L);
 			weights[i] = L;
 		}
 		Log.warning("Marginal likelihood: " + Z);


 		double max = weights[0];
 		for (double d : weights) {
 			max = Math.max(d,  max);
 		}
 		for (int i = 0; i < weights.length; i++) {
 			weights[i] = Math.exp(weights[i] - max); 			
 		} 		
 		// normalise weights so they sum to 1
 		double w = Randomizer.getTotal(weights);
 		for (int i = 0; i < NSLikelihoods.length; i++) {
 			weights[i] /= w;
 		}
 		
 		// max nr of posterior samples
 		double ESS = 0;
 		for (int i = 0; i < NSLikelihoods.length; i++) {
 			if (weights[i] > 0) {
 				ESS -= weights[i]* Math.log(weights[i]);
 			}
 		}
 		ESS = Math.exp(ESS);
 		
 		Log.warning("Max ESS: " + ESS);
        
        logln("\nCalculating statistics\n\n" + BAR);
        int stars = 0;
        int items = m_sLabels.length;
        m_fMean = newDouble(items);
        m_fStdError = newDouble(items);
        m_fStdDev = newDouble(items);
        m_fMedian = newDouble(items);
        m_f95HPDlow = newDouble(items);
        m_f95HPDup = newDouble(items);
        m_fESS = newDouble(items);
        m_fACT = newDouble(items);
        m_fGeometricMean = newDouble(items);
//        int sampleInterval = (int) (m_fTraces[0][1] - m_fTraces[0][0]);
        for (int i = 2; i < items; i++) {
            // calc mean and standard deviation
            Double[] trace = m_fTraces[i];
            double sum = 0, sum2 = 0;
            for (int k = 0; k < trace.length; k++) {
                double f = trace[k];
     			double wi = weights[k];
                sum += f * wi;
                sum2 += f * f * wi;
            }
            if (m_types[i] != type.NOMINAL) {
                m_fMean[i] = sum;// / trace.length;
                m_fStdDev[i] = Math.sqrt(sum2 /* trace.length */ - m_fMean[i] * m_fMean[i]);
            } else {
                m_fMean[i] = Double.NaN;
                m_fStdDev[i] = Double.NaN;
            }

//            if (m_types[i] == type.REAL || m_types[i] == type.INTEGER) {
//                // calc median, and 95% HPD interval
//                Double[] sorted = trace.clone();
//                Arrays.sort(sorted);
//                m_fMedian[i] = sorted[trace.length / 2];
//                // n instances cover 95% of the trace, reduced down by 1 to match Tracer
//                int n = (int) ((sorted.length - 1) * 95.0 / 100.0);
//                double minRange = Double.MAX_VALUE;
//                int hpdIndex = 0;
//                for (int k = 0; k < sorted.length - n; k++) {
//                    double range = sorted[k + n] - sorted[k];
//                    if (range < minRange) {
//                        minRange = range;
//                        hpdIndex = k;
//                    }
//                }
//                m_f95HPDlow[i] = sorted[hpdIndex];
//                m_f95HPDup[i] = sorted[hpdIndex + n];
//
//                // calc effective sample size
//                m_fACT[i] = ESS.ACT(m_fTraces[i], sampleInterval);
//                m_fStdError[i] = ESS.stdErrorOfMean(trace, sampleInterval);
//                m_fESS[i] = trace.length / (m_fACT[i] / sampleInterval);
//
//                // calc geometric mean
//                if (sorted[0] > 0) {
//                    // geometric mean is only defined when all elements are positive
//                    double gm = 0;
//                    for (double f : trace)
//                        gm += Math.log(f);
//                    m_fGeometricMean[i] = Math.exp(gm / trace.length);
//                } else
//                    m_fGeometricMean[i] = Double.NaN;
//            } else {
                m_fMedian[i] = Double.NaN;
                m_f95HPDlow[i] = Double.NaN;
                m_f95HPDup[i] = Double.NaN;
                m_fACT[i] = Double.NaN;
                m_fESS[i] = Double.NaN;
                m_fGeometricMean[i] = Double.NaN;
//            }
            while (stars < 80 * (i + 1) / items) {
                log("*");
                stars++;
            }
        }
        logln("\n");
    } // calcStats

	private Double[] newDouble(int items) {
		Double [] array = new Double[items];
		Arrays.fill(array, Double.NaN);
		return array;
	}

    public void print(PrintStream out) {
    	out.println("#Particles = " + particleCount);
    	
    	// set up header for prefix, if any is specified
    	String prefix = System.getProperty("prefix");
    	String prefixHead = (prefix == null ? "" : "prefix ");
    	if (prefix != null) {
	    	String [] p = prefix.trim().split("\\s+");
	    	if (p.length > 1) {
	    		prefixHead = "";
	    		for (int i = 0; i < p.length; i++) {
	    			prefixHead += "prefix" + i + " ";
	    		}
	    	}
    	}
    	
        try {
            // delay so that stars can be flushed from stderr
            Thread.sleep(100);
        } catch (Exception e) {
        }
        int max = 0;
        for (int i = 1; i < m_sLabels.length; i++)
            max = Math.max(m_sLabels[i].length(), max);
        String space = "";
        for (int i = 0; i < max; i++)
            space += " ";

        out.println("item" + space.substring(4) + " " + prefixHead +
        		format("mean")  + format("stddev"));
        for (int i = 2; i < m_sLabels.length; i++) {
            out.println(m_sLabels[i] + space.substring(m_sLabels[i].length()) + SPACE + (prefix == null ? "" : prefix + SPACE) +
                    format(m_fMean[i]) + SPACE + format(m_fStdDev[i]));
        }
    }

	
    static void printUsageAndExi() {
        	System.out.println("NSLogAnalyser [options] -N n1 ... n2 -log [file1] ... [filen]");
        	System.out.println("-N space separated list of particle sizes in order of files provided\n"
        			+ "If only 1 argument is provided it is assumed to apply to all files\n"
        			+ "If no -N item is provided it is assumed N=1");
            System.out.println("-quiet Quiet mode.  Avoid printing status updates to stderr.");
        	System.out.println("-help");
        	System.out.println("--help");
        	System.out.println("-h print this message");
        	System.out.println("[fileX] log file to analyse. Multiple files are allowed, each is analysed separately");
        	System.exit(0);
	}
    
	public static void main(String[] args) {
    try {
        NSLogAnalyser analyser;
        	// process args
        	int burninPercentage = 0;
            boolean quiet = false;
        	List<String> files = new ArrayList<>();
        	List<Integer> N = new ArrayList<>();
        	int i = 0;
        	while (i < args.length) {
        		String arg = args[i];
                switch (arg) {

                case "-quiet":
                    quiet = true;
                    i += 1;
                    break;
                case "-N":
                	do {
                    	i++;
                		int n = Integer.parseInt(args[i]);
                		if (n < 0) {
                			throw new IllegalArgumentException("Number of particles must be specified with the -N argument");
                		}
                		N.add(n);
                	} while (i+2 < args.length && !args[i+1].startsWith("-"));
                    i++;
                    break;
                case "-log":
                	do {
                    	i++;
                		files.add(args[i]);
                	} while (i+2 < args.length && !args[i+1].startsWith("-"));
                    i++;
                    break;

        		case "-h":
        		case "-help":
        		case "--help":
        			printUsageAndExit();
        			break;
        		default:
        			if (arg.startsWith("-")) {
        				Log.warning.println("unrecognised command " + arg);
        				printUsageAndExit();
        			}
        			files.add(arg);
        			i++;
        		}
        	}
        	if (files.size() == 0) {
        		// no file specified, open file dialog to select one
                BEASTVersion2 version = new BEASTVersion2();
                File file = Utils.getLoadFile("LogAnalyser " + version.getVersionString() + " - Select log file to analyse",
                        null, "BEAST log (*.log) Files", "log", "txt");
                if (file == null) {
                    return;
                }
                analyser = new NSLogAnalyser(file.getAbsolutePath(), burninPercentage, quiet, 1);
                analyser.calcStats();
                analyser.print(System.out);
        	} else {
        		// process files
        		if (N.size() == 0) {
        			N.add(1);
        		}
                analyser = new NSLogAnalyser(files.get(0), burninPercentage, quiet, N.get(0));
                for (int j = 1; j < files.size(); j++) {
                	String file = files.get(j);
                	NSLogAnalyser analyser2 = new NSLogAnalyser(file, burninPercentage, quiet, N.get(j % N.size()));
                	analyser = new NSLogAnalyser(analyser, analyser2);
                }
                analyser.calcStats();
                analyser.print(System.out);
                
        	}
    } catch (Exception e) {
        e.printStackTrace();
    }
}
}
