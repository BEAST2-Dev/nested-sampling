package nestedsampling.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.Parameter;
import beast.base.core.Log;
import beast.base.evolution.operator.TipDatesRandomWalker;
import beast.base.evolution.operator.TipDatesScaler;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

@Description("Operator scedule that automatically determines operator weights")
public class NSOperatorSchedule extends OperatorSchedule {
	enum WeightSchedule {STATIC, CYCLIC, AUTO};
	public Input<WeightSchedule> weighteScheduleInput = new Input<>("weightSchedule", "adapt operator weights to cycle over all state nodes till every one is changed", null, WeightSchedule.values());
	
	enum TreeWeightSchedule {x1, x1_5, x1_5_plus, x2};
	public Input<TreeWeightSchedule> treeWeightScheduleInput = new Input<>("treeWeightSchedule", "determine weight factor for nodes in the tree", TreeWeightSchedule.x1_5_plus, TreeWeightSchedule.values());
	
	public Input<Integer> cycleRefreshThresholdInput = new Input<>("cycleRefreshThreshold", "number of dimensions below which the cyclic weight schedule resets", 0);
	public Input<String> weightInput = new Input<>("weights", "space separated list of doubles. If specified will be used as state node dimensions");
	
	
	public NSOperatorSchedule() {
		// set default autoOptimise to false
		autoOptimiseInput.setValue(false, this);
		autoOptimiseInput.setRule(Validate.REQUIRED);
		weighteScheduleInput.setValue(null, this);
		weighteScheduleInput.setRule(Validate.REQUIRED);
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
	int cycleRefreshThreshold;
	
	void initialise() {
		// the following triggers the private reweightOperators, which sets up the operators
		super.selectOperator();
//		showOperatorRates(System.out);
		weightSchedule = weighteScheduleInput.get();
		cycleRefreshThreshold = cycleRefreshThresholdInput.get();
		Log.warning("weight schedule = " + weightSchedule);
		
    	// never optimise for NS
    	setAutoOptimise(autoOptimiseInput.get());
    	Log.warning("autoOptimise = " + autoOptimiseInput.get());

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
        if (weightInput.get() != null) {
        	String [] strs = weightInput.get().trim().split("\\s+");
        	if (strs.length != stateNodes.size()) {
        		throw new IllegalArgumentException("expected " + stateNodes.size() + " weights, but got " + strs.length);
        	}
            totalDimension = 0;
        	for (int i = 0; i < strs.length; i++) {
        		stateNodeDimensions[i] = Double.parseDouble(strs[i]);
        		Log.warning(stateNodes.get(i).getID() + ": " + stateNodeDimensions[i]);
            	totalDimension += stateNodeDimensions[i];
        	}
        } else {
            totalDimension = 0;
            for (int i = 0; i < stateNodeDimensions.length; i++) {
            	StateNode sn = stateNodes.get(i);
            	if (sn instanceof Parameter) {
            		stateNodeDimensions[i] = ((Parameter) sn).getDimension();
            	} else {
            		switch (treeWeightScheduleInput.get()) {
            		case x1:
                		stateNodeDimensions[i] = 1 * ((TreeInterface) sn).getNodeCount();
            			break;
            		case x1_5:
                		stateNodeDimensions[i] = 1.5 * ((TreeInterface) sn).getNodeCount();
            			break;
            		case x1_5_plus:
                		stateNodeDimensions[i] = 1.5 * ((TreeInterface) sn).getNodeCount();
                		for (BEASTInterface o : ((BEASTInterface)sn).getOutputs()) {
                			if (o instanceof TipDatesRandomWalker) {
                				TipDatesRandomWalker op = (TipDatesRandomWalker) o;
                				stateNodeDimensions[i] += op.m_taxonsetInput.get().getTaxonCount();
                			} else if (o instanceof TipDatesScaler) {
                				TipDatesScaler op = (TipDatesScaler) o;
                				stateNodeDimensions[i] += op.taxonsetInput.get().getTaxonCount();
                				
                			}
                		}
            			break;
            		case x2:
                		stateNodeDimensions[i] = 2 * ((TreeInterface) sn).getNodeCount();
            			break;
            		}
            	}
            	totalDimension += stateNodeDimensions[i];
            }
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
        double [] normalizedWeights = new double[operators.size()];
        for (int i = 0; i < operatorWeights.length; i++) {
        	for (int j = 0; j < operatorWeights[i].length; j++) {
        		int k = operators.indexOf(stateNodeOperators[i].get(j));
        		double w = (j == 0 ? operatorWeights[i][0] : operatorWeights[i][j] - operatorWeights[i][j-1]); 
        		normalizedWeights[k] += stateNodeDimensions[i] * w;
        	}
        }
//        System.out.println(Arrays.toString(normalizedWeights));
        normalizedWeights = Randomizer.getNormalized(normalizedWeights);
        setNormalizedWeights(normalizedWeights);
//        System.out.println(Arrays.toString(normalizedWeights));
        
        
        // calculate cumulative operator probs
        double [] cumulativeProbs = new double[getNormalizedWeights().length];
        cumulativeProbs[0] = getNormalizedWeights()[0];
        for (int i = 1; i < getNormalizedWeights().length; i++) {
        	cumulativeProbs[i] = cumulativeProbs[i-1] + getNormalizedWeights()[i];
        }
        setCumulativeProbs(cumulativeProbs);
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
		    		if (prevAcceptedCount != prevOperator.get_m_nNrAccepted()) {
		    			// the prevOperator was accepted, so remove 1 from stateNodeWeights
		    			stateNodeWeights[prevStateNodeIndex]--;
		    			totalAccepted++;
		    			if (totalAccepted >= totalDimension - cycleRefreshThreshold) {
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
		        prevAcceptedCount = prevOperator.get_m_nNrAccepted();
		        prevStateNodeIndex = stateNodeIndex;
		        return prevOperator;
	    	}
	    	case AUTO:
	    	{
		    	if (prevOperator != null) {
		    		if (prevAcceptedCount != prevOperator.get_m_nNrAccepted()) {
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
		        prevAcceptedCount = prevOperator.get_m_nNrAccepted();
		        prevStateNodeIndex = stateNodeIndex;
		        return prevOperator;
	    	}

    	case STATIC:
    	default:
    		// don't cycle, just pick an operator with static weight distribution
    		{
    			final int operatorIndex = Randomizer.randomChoice(getCumulativeProbs());
    			return operators.get(operatorIndex);
    		}
    	}
    }
    
}
