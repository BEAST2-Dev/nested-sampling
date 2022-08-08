package nestedsampling.evolution.speciation;

import java.util.*;

import org.apache.commons.math3.util.FastMath;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;
import test.beast.BEASTTestCase;


@Description("Yule model with normalisation constant so density integrates to 1")
public class YuleModelNormalised extends TreeDistribution {
	public Input<Double> rhoInput = new Input<>("rho", "Extant sampling proportion, default 1", 1.0);
	public Input<RealParameter> lambdaInput = new Input<>("birthDiffRate", "birth rate, one for each clade", Validate.REQUIRED);

	protected double rho;
	
	RealParameter lambda;
	
	double [] oldLength;
	double [] oldRate;
	double [] oldLength1;
	double [] oldRate1;

	Tree tree;

	public void initAndValidate() {
		tree = (Tree) treeInput.get();
		rho = rhoInput.get();
		lambda = lambdaInput.get();

		if (lambda.getDimension() != 1) {
			throw new IllegalArgumentException("birth rate input should have dimension 1");
		}
	}

	@Override
	public double calculateLogP() {
		
		if (oldLength == null) {
			initialise();
		}

		logP = 0;
		
		final double gamma = 0;

		for (Node node : tree.getNodesAsArray()) {
			final double lambda = this.lambda.getValue();
			final double t = node.getHeight();
			if (node.isRoot()) {
				if( rho != 1 ) {
					logP += -2 * FastMath.log1p(-p0(lambda, gamma, t));
				}
			} else {
				final double tp = node.getParent().getHeight();
				logP += logftip(lambda, gamma, tp, t);
				if( ! node.isLeaf() ) {
					logP += FastMath.log(lambda);
				}
			}
		}
	
		for (Node node : tree.getNodesAsArray()) {
			if (node.getLength() < 1e-6 && !node.isRoot() && node.getLength() > 0) {
				logP -= 1.0/node.getLength(); 
			}
		}
		return logP;
	}
	
	void initialise() {
		oldLength = new double[tree.getNodeCount()];
		oldRate = new double[tree.getNodeCount()];
		oldLength1 = new double[tree.getNodeCount()];
		oldRate1 = new double[tree.getNodeCount()];
	}

	protected double p0(double lambda, /* double mu, */double gamma, double time) {
		//double c = Math.sqrt(Math.pow(gamma + 0 + lambda, 2) - 4 * 0 * lambda * (1 - 0));
		//double c = gamma + lambda;
		//double x1 = -(gamma + 0 + lambda + c) / 2;
		double x1 = -(gamma + lambda);
		//double x2 = -(gamma + 0 + lambda - c) / 2; == 0
		//double A1 = lambda * (1 - rho) + x1;
		double A1 = -lambda * rho - gamma;
		double A2 = lambda * (1 - rho);
		//double res = (A1 * 0 - A2 * x1 * Math.exp(-c * time)) / (lambda * (A2 * Math.exp(-c * time) - A1));
		double res = (A2 * x1 * FastMath.exp(x1 * time)) / (lambda * (A2 * FastMath.exp(x1 * time) - A1));
		return res;
	}

	protected double logftip(double lambda, /* double mu, */double gamma, double time) {
		//double c = Math.sqrt(Math.pow(gamma + 0 + lambda, 2) - 4 * 0 * lambda * (1 - 0));
		//double c = gamma + lambda;
		double x1 = -(gamma + lambda);
		//double x2 = -(gamma + 0 + lambda - c) / 2; == 0
		//double A1 = lambda * (1 - rho) + x1;
		double A1 = -lambda * rho - gamma;
		double A2 = lambda * (1 - rho);
		double y = A2 * FastMath.exp(x1 * time) - A1;
		double res = x1 * x1 / (FastMath.exp(-x1 * time) * y * y);
		return FastMath.log(res);
	}

	protected double logftip(double lambda, /* double mu, */double gamma, double time1, double time2) {
		double v;
		final double x1 = -(gamma + lambda);
		if( rho == 1 ) {
			v = x1 * (time1 - time2);
		} else {
			final double mlamrho = -lambda * rho;
			final double A1 = mlamrho - gamma;
			final double A2 = lambda + mlamrho;

			final double x1t1 = x1 * time1;
			final double y1 = A2 * FastMath.exp(x1t1) - A1;

			final double x1t2 = x1 * time2;
			final double y2 = A2 * FastMath.exp(x1t2) - A1;

			v = (x1t1 - x1t2) + 2 * FastMath.log(y2 / y1);
		}
		return v;
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}
	
	
	@Override
	protected boolean requiresRecalculation() {
		return true;
	}

	
	public static void main(String[] args) throws Exception {
		Alignment data = BEASTTestCase.getAlignment();
	    Tree tree = BEASTTestCase.getTree(data);
	    
		YuleModelNormalised myd = new YuleModelNormalised();
		myd.initByName("tree", tree, 
				"newick", "(human:0.024003,chimp:0.010772,bonobo:0.010772),gorilla:0.036038,orangutan:0.069125,siamang:0.099582;",
				"birthDiffRate", "0.1",
				"gamma", "0.5");
		
		System.err.println("logP = " + myd.calculateLogP());
	}
}
