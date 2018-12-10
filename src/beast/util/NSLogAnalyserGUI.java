package beast.util;

import java.util.ArrayList;
import java.util.List;

import beast.app.util.LogFile;
import beast.app.util.OutFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.Input.Validate;

@Description("GUI wrapper for NSLogAnalyser")
public class NSLogAnalyserGUI extends Runnable {
	public Input<LogFile> logFileInput = new Input<>("log","input log file produced with nested sampling analysis", Validate.REQUIRED);
	public Input<Integer> NInput = new Input<>("N","particle size used  in generating the log file", 1);
	public Input<OutFile> outputInput = new Input<>("out","output file. Print to stdout if not specified", Validate.REQUIRED);
	public Input<Boolean> quietInput = new Input<>("quiet","Quiet mode.  Avoid printing status updates to stderr.", false);
	public Input<Boolean> noposteriorInput = new Input<>("noposterior","do not produce posterior", false);
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		List<String> args = new ArrayList<>();
		args.add("-log");
		args.add(logFileInput.get().getAbsolutePath());
		
		args.add("-N");
		args.add(NInput.get() + "");
		
		if (quietInput.get()) {
			args.add("-quiet");
		}
		
		if (noposteriorInput.get()) {
			args.add("-noposterior");
		}
		
		args.add("-out");
		args.add(outputInput.get().getAbsolutePath());
		
		NSLogAnalyser.main(args.toArray(new String[]{}));
	}
	
}