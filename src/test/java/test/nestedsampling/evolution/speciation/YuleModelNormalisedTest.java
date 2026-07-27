//package test.nestedsampling.evolution.speciation;
//
//import beast.base.evolution.alignment.Alignment;
//import beast.base.evolution.tree.Tree;
//import nestedsampling.evolution.speciation.YuleModelNormalised;
//import org.junit.jupiter.api.Test;
//import test.beast.BEASTTestCase;
//
//public class YuleModelNormalisedTest {
//
//	@Test
//	public void testYuleModelNormalised() throws Exception {
//		Alignment data = BEASTTestCase.getAlignment();
//	    Tree tree = BEASTTestCase.getTree(data);
//
//		YuleModelNormalised myd = new YuleModelNormalised();
//		myd.initByName("tree", tree,
//				"newick", "(human:0.024003,chimp:0.010772,bonobo:0.010772),gorilla:0.036038,orangutan:0.069125,siamang:0.099582;",
//				"birthDiffRate", "0.1",
//				"gamma", "0.5");
//
//		System.err.println("logP = " + myd.calculateLogP());
//	}
//}
// TODO: required beast-base-2.8.0-tests.jar,
//       see https://github.com/CompEvol/beast3/issues/127