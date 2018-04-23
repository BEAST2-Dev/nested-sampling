package beast.core;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.StateNode;
import beast.core.parameter.Parameter;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

@Description("Operator scedule tuned for NS")
public class NSOperatorSchedule extends OperatorSchedule {
	
	
	List<StateNode> stateNodes;
	double [] stateNodeWeights;
	double [] stateNodeDimensions;
	double totalDimension;
	double totalAccepted;
	List<Operator> [] stateNodeOperators;
	double [][] operatorWeights;
	int consecutiveRejectCount;
	
	Operator prevOperator = null;
	int prevAcceptedCount, prevStateNodeIndex;
	
	void initialise() {		
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
        
        // initialise operatorsWeights
        operatorWeights = new double[stateNodes.size()][];
        for (int i = 0; i < stateNodes.size(); i++) {
        	List<Operator> operators = stateNodeOperators[i];
        	operatorWeights[i] = new double[operators.size()];
        	for (int j = 0; j < operators.size(); j++) {
        		operatorWeights[i][j] = operators.get(j).getWeight();
        	}
        	operatorWeights[i] = Randomizer.getNormalized(operatorWeights[i]);
        	for (int j = 1; j < operatorWeights[i].length; j++) {
        		operatorWeights[i][j] += operatorWeights[i][j-1]; 
        	}
        
        }

        
    	prevOperator = null;
    	prevAcceptedCount = 0;
    	totalAccepted = 0;
    	consecutiveRejectCount = 0;
	}

	
    public Operator selectOperator() {
    	if (stateNodes == null) {
    		initialise();
    	}
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


    @Override
    public void showOperatorRates(PrintStream out) {
    	out.println("NSOperatorSchedule::showOperatorRates(); not implemented yet");
    }
}
