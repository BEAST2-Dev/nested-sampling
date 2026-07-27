package test.beast.integration;

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;

/**
 * Shared test infrastructure for beast-base integration tests.
 *
 * <p>Two distinct concerns are kept as separate methods on purpose:
 * <ul>
 *   <li>{@link #resolveExamplesDir()} — pure function; finds where BEAST <em>reads</em>
 *       XML/JSON input examples from the test classpath.</li>
 *   <li>{@link #setUpOutputDir()} — side-effectful; creates the {@code ./test/} directory
 *       and sets {@code file.name.prefix} so BEAST <em>writes</em> log/tree output there.
 *       Call from {@code @BeforeEach}.</li>
 * </ul>
 * Merging them would couple a pure query to a mutating side effect, forcing every
 * caller of {@code resolveExamplesDir()} to trigger directory creation implicitly.
 *
 * <p>Individual test classes are responsible for naming their own XML/JSON files.
 */
public class XMLPathUtil {

    private static final String EXAMPLES_CLASSPATH = "beast.base/examples";

    /**
     * Returns the absolute path to the beast.base examples directory.
     * Resolves via the test classpath (works on any machine or CI runner),
     * falling back to {@code user.dir} if the resource is not found.
     */
    public static String resolveExamplesDir() {
        URL url = XMLPathUtil.class.getClassLoader().getResource(EXAMPLES_CLASSPATH);
        if (url != null) {
            try {
                return new File(url.toURI()).getAbsolutePath();
            } catch (URISyntaxException e) {
                // fall through to user.dir fallback
            }
        }
        return System.getProperty("user.dir") + "/" + EXAMPLES_CLASSPATH;
    }

    /**
     * Creates the {@code ./test/} output directory if absent and sets
     * {@code file.name.prefix=test/} so BEAST logger output is written there.
     * Call from {@code @BeforeEach} in each integration test class.
     */
    public static void setUpOutputDir() {
        setUpOutputDir("");
    }

    /**
     * Creates {@code ./test/<subdir>/} and sets {@code file.name.prefix} to that
     * path, isolating log and tree files for one specific test from those of others.
     * Use when multiple tests in the same class share log-file names (e.g. when
     * XMLs all write to {@code test.$(seed).log}).
     *
     * <p>Note: {@code file.name.prefix} is a JVM-wide system property. Callers that
     * mutate it must declare {@code @ResourceLock("beast.logger.globals")} so JUnit 5's
     * parallel scheduler serializes them; see {@code junit-platform.properties} for the
     * full list of affected classes.
     */
    public static void setUpOutputDir(String subdir) {
        String path = subdir == null || subdir.isEmpty() ? "test/" : "test/" + subdir + "/";
        File dir = new File("./" + path);
        if (!dir.exists())
            dir.mkdirs();
        System.setProperty("file.name.prefix", path);
    }
}
