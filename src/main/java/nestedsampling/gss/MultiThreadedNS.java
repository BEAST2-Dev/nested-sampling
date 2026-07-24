package nestedsampling.gss;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beast.base.parser.XMLProducer;

@Description("Multi threaded nested sampling")
public class MultiThreadedNS extends NS {
	final public Input<Integer> threadsInput = new Input<>("threads","number of threads to use", 2);
	
	/** set of particles for all threads **/
	Map<String, Double> particlePool;

	/** deals with running and syncing threads **/
	int threadCount;
	ExecutorService exec;
	CountDownLatch countDown;
	CoreRunnable [] runnable;

	/** the classes doing the actual work **/
	NSThread [] NS;
	
	
	@Override
	public void initAndValidate() {
		particlePool = new LinkedHashMap<>();
		threadCount = threadsInput.get();
	    exec = Executors.newFixedThreadPool(threadCount);
	    
	    
	    XMLProducer xmlProducer = new XMLProducer();
	    String xml = xmlProducer.toRawXML(this);
	    xml = "<beast version='2.5'>\n" +
	    		xml.replaceAll("spec='" + this.getClass().getCanonicalName() + "'",
	    		"spec='" + NSThread.class.getCanonicalName() + "'").
	    		replaceAll("threads='" + threadCount+"'", "")
	    	+ "\n</beast>";
	    
	    NS = new NSThread[threadCount];
	    runnable = new CoreRunnable[threadCount];
	    
	    for (int i = 0; i < threadCount; i++) {
	    	try {
	    	    XMLParser xmlParser = new XMLParser();
	    		String xml2 = xml.replaceAll("fileName='", "fileName='" + i);
	    		Object o = xmlParser.parseFragment(xml2, true);
				NS[i] = (NSThread) o;
				runnable[i] = new CoreRunnable(NS[i]);
			} catch (XMLParserException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	    
	}
	
	@Override
	public void run() throws IOException, SAXException, ParserConfigurationException {
		
		// start up threads, perform initialisations, halt on updating of particles
        countDown = new CountDownLatch(threadCount);
        for (int i = 0; i < threadCount; i++) {
        	NS[i].countDown = countDown;
        	NS[i].particlePool = particlePool;
        }
        
        for (CoreRunnable coreRunnable : runnable) {
            exec.execute(coreRunnable);
        }
        
        // wait till all threads are ready to update particles
        try {
			countDown.await();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
//        
//        // carry on till all threads finished.
//        // threads are synchronised every time a particle was updated
//        // TODO: perhaps syncing is not necessary, making things more efficient?
//        int particlesFinished = 0;
//        while (particlesFinished < threadCount) {
//            countDown = new CountDownLatch(threadCount - particlesFinished);
//            for (int i = 0; i < threadCount; i++) {
//            	if (!NS[i].isFinished()) { 
//                	NS[i].countDown = countDown;
//            		NS[i].nsCountDown.countDown();
//            	}
//            }
//        	
//            try {
//    			countDown.await();
//    		} catch (InterruptedException e) {
//    			e.printStackTrace();
//    		}
//        	particlesFinished = 0;
//            for (int i = 0; i < threadCount; i++) {
//            	if (NS[i].isFinished()) { 
//            		particlesFinished++;
//            	}
//            }
//        }
    }
	

   class CoreRunnable implements Runnable {
	   NSThread NS;

        CoreRunnable(NSThread NS) {
            this.NS = NS;
        }

        @Override
		public void run() {
            try {
            	NS.run();
            } catch (Exception e) {
                e.printStackTrace();
                Log.err("Something went wrong running thread: " + e.getMessage());
                NS.countDown.countDown();
            }
        }

    } // CoreRunnable
}
