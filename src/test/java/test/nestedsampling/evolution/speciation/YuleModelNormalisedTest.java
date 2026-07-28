package test.nestedsampling.evolution.speciation;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Tree;
import nestedsampling.evolution.speciation.YuleModelNormalised;
import org.junit.jupiter.api.Test;
import test.beast.BEASTTestCase;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * The expected logP is computed from beast3 code on July 2026.
 */
public class YuleModelNormalisedTest {

	@Test
	public void testYuleModelNormalised() throws Exception {
		Alignment data = BEASTTestCase.getAlignment();
	    Tree tree = BEASTTestCase.getTree(data);

		YuleModelNormalised myd = new YuleModelNormalised();
		myd.initByName("tree", tree,
//				"newick", "(human:0.024003,chimp:0.010772,bonobo:0.010772),gorilla:0.036038,orangutan:0.069125,siamang:0.099582;",
				"birthDiffRate", "0.1",
				"rho", "0.5");
// logP = -10.018014963613476
		double logP = myd.calculateLogP();
		System.out.println("logP = " + logP);

		assertEquals(-10.018015, logP, 1.0E6);
	}
}