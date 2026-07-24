open module nested.sampling {
    requires beast.pkgmgmt;
    requires beast.base;
    requires java.xml;

    requires org.apache.commons.statistics.distribution;
    requires org.apache.commons.math4.core;

    requires model.selection;

    // GUI (optional at runtime)
    requires static beast.fx;
    requires static javafx.controls;

    exports nestedsampling.core;
    exports nestedsampling.evolution.speciation;
    exports nestedsampling.gss;
    exports nestedsampling.util;

    provides beast.base.core.BEASTInterface with
        nestedsampling.core.NSLogger,
        nestedsampling.core.NSOperatorSchedule,
        nestedsampling.evolution.speciation.YuleModelNormalised,
        nestedsampling.gss.DynamicNestedSampling,
        nestedsampling.gss.MultiThreadedNS,
        nestedsampling.gss.NIS,
        nestedsampling.gss.NS,
        nestedsampling.gss.NSThread,
        nestedsampling.util.MCMC2NIS,
        nestedsampling.util.MCMC2NS,
        nestedsampling.util.NSLogAnalyserGUI;
}
