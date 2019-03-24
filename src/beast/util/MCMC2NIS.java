package beast.util;

import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Runnable;
import beast.gss.MCMC2IS;
import beast.gss.NS;

@Description("Convert MCMC analysis to nested importance sampling analysis and optionally run the analysis")
public class MCMC2NIS extends MCMC2IS {
	public Input<Integer> particleCountInput = new Input<>("particleCount", "number of particles (default 100)", 100);
	public Input<Integer> subChainLengthInput = new Input<>("subChainLength",
			"number of MCMC samples for each epoch (default 1000)", 1000);
	public Input<Double> epsilonInput = new Input<>("epsilon", "stopping criterion: smallest change in ML estimate to accept", 1e-6);
	
	
	@Override
	protected Runnable newInstance(MCMC mcmc) {		
		NS ns = new NS();
		ns.particleCountInput.setValue(particleCountInput.get(), ns);
		ns.subChainLengthInput.setValue(subChainLengthInput.get(), ns);
		ns.epsilonInput.setValue(epsilonInput.get(), ns);
		
		for (Input<?> input : mcmc.listInputs()) {
			ns.setInputValue(input.getName(), mcmc.getInput(input.getName()).get());
		}
		
		return ns;
	}

	public static void main(String[] args) throws Exception {
		new Application(new MCMC2NIS(), "MCMC2NIS", args);
	}
}
