package test.beast.gss.integration;


import java.io.File;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.junit.Test;
import org.xml.sax.SAXException;

import beast.core.Logger;
import beast.util.Randomizer;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import junit.framework.TestCase;
import nestedsampling.gss.NS;

public class NS4TaxaTest extends TestCase {
	
	@Test
	public void testNS4Taxa() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
        Logger.FILE_MODE = Logger.LogFileMode.overwrite;

        int seed = 127;
        String fileName = "examples/NS_4taxa_NormalBirthRate.xml";
        Randomizer.setSeed(seed);
        System.out.println("Processing " + fileName);
        XMLParser parser = new XMLParser();
        beast.core.Runnable runable = parser.parseFile(new File(fileName));
        if (runable instanceof NS) {
            NS ns = (NS) runable;
            ns.run();
            double Z = ns.getEvidence();
            assertEquals(Z, -2349.5333124902536, 1.0);
        }
        System.out.println("Done " + fileName);
	}

}
