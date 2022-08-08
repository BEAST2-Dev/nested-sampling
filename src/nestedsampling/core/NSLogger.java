package nestedsampling.core;


import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;

import beast.base.core.Description;
import beast.base.core.Loggable;
import beast.base.inference.Logger;
import beast.base.parser.XMLProducer;

@Description("Trace logger that logs NSLikelihood as second column (used by NSLogAnalyser to prodoce posterior sample)")
public class NSLogger extends Logger {
	
	@Override
	public void init() throws IOException {
        final boolean needsHeader = openLogFile();
        if (needsHeader) {
            if (modelInput.get() != null) {
                // print model at top of log
                String xml = new XMLProducer().modelToXML(modelInput.get());
                xml = "#" + xml.replaceAll("\\n", "\n#");
                m_out.println("#\n#model:\n#");
                m_out.println(xml);
                m_out.println("#");
            }
            ByteArrayOutputStream baos = null;
            if (m_out == System.out) {
                baos = new ByteArrayOutputStream();
                m_out = new PrintStream(baos);
            }
            final ByteArrayOutputStream rawbaos = new ByteArrayOutputStream();
            final PrintStream out = new PrintStream(rawbaos);
            if (mode == LOGMODE.compound) {
                out.print("Sample\t");
                out.print("NSlikelihood\t");
            }
            for (final Object m_logger : loggersInput.get()) {
                ((Loggable)m_logger).init(out);
            }

            // Remove trailing tab from header
            String header = rawbaos.toString().trim();

            if (sanitiseHeadersInput.get()) {
            	m_out.print(sanitiseHeader(header));
            } else {
            	m_out.print(header);
            }
            
            m_out.println();
        }
	}
	
	public void log(int sampleNr, double NSlikelihood) {
        if (sampleNr < 0) {
            return;
        }
        if (sampleOffset >= 0) {
            if (sampleNr == 0) {
                // don't need to duplicate the last line in the log
                return;
            }
            sampleNr += sampleOffset;
        }

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream out = new PrintStream(baos);

        if (mode == LOGMODE.compound) {
            out.print((sampleNr) + "\t");
            out.print(NSlikelihood + "\t");
        }

        for (final Object m_logger : loggersInput.get()) {
        	((Loggable)m_logger).log(sampleNr, out);
        }

        // Acquire log string and trim excess tab
        String logContent;
        try {
            logContent = baos.toString("ASCII").trim();
        } catch (UnsupportedEncodingException e) {
            throw new RuntimeException("ASCII string encoding not supported: required for logging!");
        }

        m_out.println(logContent);
	}

}
