package beast.util;

import java.io.FileWriter;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.OperatorSchedule;
import beast.core.Runnable;
import beast.core.Input.Validate;
import beast.core.Logger;
import beast.core.util.Log;
import beast.gss.NS;

@Description("Convert MCMC analysis to nested sampling analysis")
public class MCMC2NS extends Runnable {
	
	public Input<XMLFile> model1Input = new Input<>("xml",
			"file name of BEAST XML file containing the model for which to create a GSS XML file for",
			new XMLFile("examples/normalTest-1XXX.xml"), Validate.REQUIRED);
	public Input<OutFile> outputInput = new Input<>("output", "where to save the file", new OutFile("beast.xml"));
	
	public Input<Integer> particleCountInput = new Input<>("particleCount", "number of particles (default 10)", 10);
	public Input<Integer> subChainLengthInput = new Input<>("subChainLength",
			"number of MCMC samples for each epoch (default 10000)", 10000);
	public Input<Double> epsilonInput = new Input<>("epsilon", "stopping criterion: smallest change in ML estimate to accept", 1e-6);
	
	
	@Override
	public void initAndValidate() {
	}

	
	@Override
	public void run() {
		XMLParser parser = new XMLParser();
		MCMC mcmc;
		try {
			mcmc = (MCMC) parser.parseFile(model1Input.get());
			
			NS ns = new NS();
			// nested sampling options
			ns.particleCountInput.setValue(particleCountInput.get(), ns);
			ns.subChainLengthInput.setValue(subChainLengthInput.get(), ns);
			ns.epsilonInput.setValue(epsilonInput.get(), ns);
			
			// mcmc options
			ns.posteriorInput.setValue(mcmc.posteriorInput.get(), ns);
			ns.startStateInput.setValue(mcmc.startStateInput.get(), ns);
			ns.operatorsInput.setValue(mcmc.operatorsInput.get(), ns);
			ns.chainLengthInput.setValue(mcmc.chainLengthInput.get(), ns);
			ns.initialisersInput.setValue(mcmc.initialisersInput.get(), ns);			

			// loggers
			ns.loggersInput.setValue(mcmc.loggersInput.get(), ns);
			for (Logger logger : ns.loggersInput.get()) {
				String fn = logger.fileNameInput.get();
				if (fn != null && fn.length() > 0 && fn.contains(".")) {
					fn = fn.trim();
					int i = fn.lastIndexOf('.');
					fn = fn.substring(0, i) + "-NS." + fn.substring(i + 1);
					logger.fileNameInput.setValue(fn, logger);
				}
			}
			
			// operator schedule
			OperatorSchedule os = mcmc.operatorScheduleInput.get();
			ns.operatorScheduleInput.setValue(os, ns);
			os.autoOptimiseInput.setValue(false, os);

		
	        Log.warning("Writing to file " + outputInput.get().getPath());
			XMLProducer producer = new XMLProducer();
			String xml = producer.toXML(ns);
	        FileWriter outfile = new FileWriter(outputInput.get());
	        outfile.write(xml);
	        outfile.close();
	        
			
		} catch (SAXException | IOException | ParserConfigurationException | XMLParserException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		
        Log.warning("Done");
	}

	public static void main(String[] args) throws Exception {
		new Application(new MCMC2NS(), "MCMC2NS", args);
	}

}
