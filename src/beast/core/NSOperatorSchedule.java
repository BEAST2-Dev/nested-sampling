package beast.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.StateNode;
import beast.core.parameter.Parameter;
import beast.core.util.Log;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

@Description("Operator scedule tuned for NS")
public class NSOperatorSchedule extends OperatorSchedule {
	enum WeightSchedule {STATIC, CYCLIC, AUTO};
	public Input<WeightSchedule> cycleInput = new Input<>("weightSchedule", "adapt operator weights to cycle over all state nodes till every one is changed", WeightSchedule.AUTO, WeightSchedule.values());
	
	public NSOperatorSchedule() {
		// set default autoOptimise to false
		autoOptimiseInput.setValue(false, this);
	}
	
	List<StateNode> stateNodes;
	double [] stateNodeWeights;
	double [] stateNodeAcceptCounts;
	double [] stateNodeDimensions;
	double totalDimension;
	double totalAccepted;
	List<Operator> [] stateNodeOperators;
	double [][] operatorWeights;
	int consecutiveRejectCount;
	
	Operator prevOperator = null;
	int prevAcceptedCount, prevStateNodeIndex;
	
	WeightSchedule weightSchedule;
	
	void initialise() {
		// the following triggers the private reweightOperators, which sets up the operators
		super.selectOperator();
//		showOperatorRates(System.out);
		weightSchedule = cycleInput.get();
		Log.warning("weight schedule = " + weightSchedule);
		
    	// never optimise for NS
    	autoOptimise = autoOptimiseInput.get();
    	Log.warning("autoOptimise = " + autoOptimise);

    	// collect stateNodes
		stateNodes = new ArrayList<>();		
        for (final Operator op : operators) {
        	List<StateNode> nodes = op.listStateNodes();
        	for (StateNode sn : nodes) {
        		if (!stateNodes.contains(sn)) {
        			stateNodes.add(sn);
        		}
        	}
        }
        
        // initialise dimensions
        stateNodeDimensions = new double[stateNodes.size()];
        totalDimension = 0;
        for (int i = 0; i < stateNodeDimensions.length; i++) {
        	StateNode sn = stateNodes.get(i);
        	if (sn instanceof Parameter) {
        		stateNodeDimensions[i] = ((Parameter) sn).getDimension();
        	} else {
        		stateNodeDimensions[i] = 2 * ((TreeInterface) sn).getNodeCount();
        	}
        	totalDimension += stateNodeDimensions[i];
        }
        
        stateNodeAcceptCounts = new double[stateNodes.size()];
        Arrays.fill(stateNodeAcceptCounts, 1.0);
        Log.warning("#parameters = " + totalDimension);
        
        // initialise stateNodeWeights
        stateNodeWeights = new double[stateNodes.size()];
        System.arraycopy(stateNodeDimensions, 0, stateNodeWeights, 0, stateNodeWeights.length);
        
        // initialise stateNodeOperators
        stateNodeOperators = new List[stateNodes.size()];
        for (int i = 0; i < stateNodeOperators.length; i++) {
        	stateNodeOperators[i] = new ArrayList<>();
        }
        for (final Operator op : operators) {
        	List<StateNode> nodes = op.listStateNodes();
        	for (StateNode sn : nodes) {
            	int i = stateNodes.indexOf(sn);
            	stateNodeOperators[i].add(op);
        	}
        }
        
        // determine total sum of dimensions of parameters operators operate on
        int [] operatorDimensions = new int[operators.size()];
        for (int i = 0; i < operators.size(); i++) {
        	for (StateNode sn : operators.get(i).listStateNodes()) {
            	if (sn instanceof Parameter) {
            		operatorDimensions[i] += ((Parameter) sn).getDimension();
            	} else {
            		operatorDimensions[i] += 2 * ((TreeInterface) sn).getNodeCount();
            	}        		 
        	}
        }
        
        // initialise operatorsWeights
        operatorWeights = new double[stateNodes.size()][];
        for (int i = 0; i < stateNodes.size(); i++) {
        	List<Operator> sOperators = stateNodeOperators[i];
        	operatorWeights[i] = new double[sOperators.size()];
        	double dim = stateNodeDimensions[i];
        	for (int j = 0; j < sOperators.size(); j++) {
        		Operator operator = sOperators.get(j);
        		int k = operators.indexOf(operator);
        		operatorWeights[i][j] = operator.getWeight() * dim / operatorDimensions[k];
        	}
        	operatorWeights[i] = Randomizer.getNormalized(operatorWeights[i]);
        	for (int j = 1; j < operatorWeights[i].length; j++) {
        		operatorWeights[i][j] += operatorWeights[i][j-1]; 
        	}
        
        }

        // calculate normalised weights
        normalizedWeights = new double[operators.size()];
        for (int i = 0; i < operatorWeights.length; i++) {
        	for (int j = 0; j < operatorWeights[i].length; j++) {
        		int k = operators.indexOf(stateNodeOperators[i].get(j));
        		double w = (j == 0 ? operatorWeights[i][0] : operatorWeights[i][j] - operatorWeights[i][j-1]); 
        		normalizedWeights[k] += stateNodeDimensions[i] * w;
        	}
        }
//        System.out.println(Arrays.toString(normalizedWeights));
        normalizedWeights = Randomizer.getNormalized(normalizedWeights);
//        System.out.println(Arrays.toString(normalizedWeights));
        
        
        // calculate cumulative operator probs
        cumulativeProbs = new double[normalizedWeights.length];
        cumulativeProbs[0] = normalizedWeights[0];
        for (int i = 1; i < normalizedWeights.length; i++) {
        	cumulativeProbs[i] = cumulativeProbs[i-1] + normalizedWeights[i];
        }
//        System.out.println(Arrays.toString(cumulativeProbs));

        prevOperator = null;
    	prevAcceptedCount = 0;
    	totalAccepted = 0;
    	consecutiveRejectCount = 0;

//    	showOperatorRates(System.out);
    	
    }

	
    public Operator selectOperator() {
    	if (stateNodes == null) {
    		initialise();
    	}
    	
    	switch (weightSchedule) {
    	case CYCLIC :
	    	{
		    	if (prevOperator != null) {
		    		if (prevAcceptedCount != prevOperator.m_nNrAccepted) {
		    			// the prevOperator was accepted, so remove 1 from stateNodeWeights
		    			stateNodeWeights[prevStateNodeIndex]--;
		    			totalAccepted++;
		    			if (totalAccepted == totalDimension) {
		    				// reset
		    				totalAccepted = 0;
		        			consecutiveRejectCount = 0;
		    		        System.arraycopy(stateNodeDimensions, 0, stateNodeWeights, 0, stateNodeWeights.length);
		    			}
		    		} else {
		    			consecutiveRejectCount++;
		    			if (consecutiveRejectCount > totalDimension) {
		    				// reset
		    				totalAccepted = 0;
		        			consecutiveRejectCount = 0;
		    		        System.arraycopy(stateNodeDimensions, 0, stateNodeWeights, 0, stateNodeWeights.length);
		    			}
		    		}
		    	}
		    	
		        final int stateNodeIndex = Randomizer.randomChoicePDF(stateNodeWeights);
		        double [] ow = operatorWeights[stateNodeIndex];
		        final int operatorIndex = ow.length == 1 ? 0 : Randomizer.randomChoice(ow);
		        prevOperator = stateNodeOperators[stateNodeIndex].get(operatorIndex);
		        prevAcceptedCount = prevOperator.m_nNrAccepted;
		        prevStateNodeIndex = stateNodeIndex;
		        return prevOperator;
	    	}
	    	case AUTO:
	    	{
		    	if (prevOperator != null) {
		    		if (prevAcceptedCount != prevOperator.m_nNrAccepted) {
		    			// the prevOperator was accepted, so add 1 from stateNodeAcceptCounts
		    			stateNodeAcceptCounts[prevStateNodeIndex]++;
		    		}
		    	}
		    	
		    	double [] weights = new double[stateNodeWeights.length];
		    	for (int i = 0; i < weights.length; i++) {
		    		weights[i] = stateNodeWeights[i] / stateNodeAcceptCounts[i];
		    	}
		    	
		        final int stateNodeIndex = Randomizer.randomChoicePDF(weights);
		        double [] ow = operatorWeights[stateNodeIndex];
		        final int operatorIndex = ow.length == 1 ? 0 : Randomizer.randomChoice(ow);
		        prevOperator = stateNodeOperators[stateNodeIndex].get(operatorIndex);
		        prevAcceptedCount = prevOperator.m_nNrAccepted;
		        prevStateNodeIndex = stateNodeIndex;
		        return prevOperator;
	    	}

    	case STATIC:
    	default:
    		// don't cycle, just pick an operator with static weight distribution
    		{
    			final int operatorIndex = Randomizer.randomChoice(cumulativeProbs);
    			return operators.get(operatorIndex);
    		}
    	}
    }
    
}
