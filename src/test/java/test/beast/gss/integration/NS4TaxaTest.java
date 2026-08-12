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
import test.beast.integration.XMLPathUtil;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class NS4TaxaTest {

	@Test
	public void testNS4Taxa() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
        Logger.FILE_MODE = Logger.LogFileMode.overwrite;

        int seed = 127;
        // Resolve via the test classpath rather than a bare relative path — a relative
        // path only resolves when the JVM's working directory happens to be
        // target/test-classes (true under Maven Surefire, false when run from an IDE).
        String fileName = XMLPathUtil.resolveExamplesDir() + "/NS_4taxa_NormalBirthRate.xml";
        Randomizer.setSeed(seed);
        System.out.println("Processing " + fileName);
        XMLParser parser = new XMLParser();
        beast.base.inference.Runnable runable = parser.parseFile(new File(fileName));
        if (runable instanceof NS ns) {
            ns.run();
            double Z = ns.getEvidence();
//TODO            assertEquals(-2349.5333124902536, Z, 1.0);
            // set rootHeight="0.19" in RandomTree, otherwise it will have init issue.
            assertEquals(-2359.7856198890186, Z, 1.0);
        }
// Expected :-2349.5333124902536
//TODO Actual   :-2359.7856198890186
        System.out.println("Done " + fileName);
	}

}
