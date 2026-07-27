package test.beast.gss.integration;


import java.io.File;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;

import org.junit.jupiter.api.Test;
import org.xml.sax.SAXException;

import beast.base.inference.Logger;
import beast.base.util.Randomizer;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import nestedsampling.gss.NS;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class NS4TaxaTest {

	@Test
	public void testNS4Taxa() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
        Logger.FILE_MODE = Logger.LogFileMode.overwrite;

        int seed = 127;
        String fileName = "nestedsampling/examples/NS_4taxa_NormalBirthRate.xml";
        Randomizer.setSeed(seed);
        System.out.println("Processing " + fileName);
        XMLParser parser = new XMLParser();
        beast.base.inference.Runnable runable = parser.parseFile(new File(fileName));
        if (runable instanceof NS ns) {
            ns.run();
            double Z = ns.getEvidence();
            assertEquals(-2349.5333124902536, Z, 1.0);
        }
        System.out.println("Done " + fileName);
	}

}
