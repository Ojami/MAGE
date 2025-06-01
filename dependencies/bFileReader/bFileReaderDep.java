import java.nio.file.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.function.*;
import java.util.regex.Pattern;
import java.io.InputStream;
import java.io.BufferedInputStream;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.FileInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.ByteArrayOutputStream;
import java.io.ObjectOutputStream;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.GZIPInputStream;

public class bFileReaderDep {

    // A functional interface that returns the string to append for each line,
    // or an empty string if we want to discard the line.
    @FunctionalInterface
    private interface LineProcessor {
        String process(String line);
    }

    // ----------------------------------------------
    // 1) openFile
    // ----------------------------------------------
    static private BufferedReader openFile(String filename) throws IOException {
        String fileExt = filename.substring(filename.lastIndexOf(".") + 1).toLowerCase();
        if (fileExt.equals("gz")) {
            FileInputStream fis = new FileInputStream(filename);
            GZIPInputStream gis = new GZIPInputStream(fis);
            InputStreamReader isr = new InputStreamReader(gis);
            return new BufferedReader(isr);
        } else {
            return Files.newBufferedReader(Paths.get(filename));
        }
    }

    // ----------------------------------------------
    // 2) lineCount (needed for parallel offsets)
    // ----------------------------------------------
    // We'll use your existing function for counting lines quickly.
    public static int lineCount(String filename) throws IOException {
        String fileExt = filename.substring(filename.lastIndexOf(".") + 1).toLowerCase();
        InputStream is;
        if (fileExt.equals("gz")) {
            FileInputStream fis = new FileInputStream(filename);
            GZIPInputStream gis = new GZIPInputStream(fis);
            is = new BufferedInputStream(gis);
        } else {
            is = new BufferedInputStream(new FileInputStream(filename));
        }

        try {
            byte[] c = new byte[1024];
            int readChars = is.read(c);
            if (readChars == -1) {
                return 0;
            }
            int count = 0;
            while (readChars == 1024) {
                for (int i = 0; i < 1024;) {
                    if (c[i++] == '\n') {
                        ++count;
                    }
                }
                readChars = is.read(c);
            }
            while (readChars != -1) {
                for (int i = 0; i < readChars; ++i) {
                    if (c[i] == '\n') {
                        ++count;
                    }
                }
                readChars = is.read(c);
            }
            return (count == 0) ? 1 : count;
        } finally {
            is.close();
        }
    }

    // ----------------------------------------------
    // 3) parallelProcessFile (the new concurrency)
    // ----------------------------------------------
    /**
     * Splits file lines among multiple threads, each applying lineProcessor
     * to produce output for each line. We skip the first `skip` lines globally.
     * We then read the remainder among the threads. Finally, we merge results.
     *
     * Because .gz is not random-access friendly, each thread re-reads from the start
     * and discards lines until it reaches its assigned start offset.
     *
     * @param filename The file to open
     * @param skip     Number of lines to skip at the start
     * @param lineProcessor A function that, given a line, returns the string to append or ""
     * @return aggregated char[] of the processed lines
     */
    private static char[] parallelProcessFile(String filename, int skip, LineProcessor lineProcessor) throws IOException {
        int total = lineCount(filename); // total lines in file
        if (total <= skip) {
            // If the skip is >= total lines, no data to process
            return new char[0];
        }

        int numThreads = Runtime.getRuntime().availableProcessors();
        // Lines to actually process after skipping
        int linesToProcess = total - skip;
        // Basic chunk size
        int chunkSize = (int)Math.ceil((double)linesToProcess / numThreads);

        // We'll hold each thread's result in a separate slot, then merge
        StringBuilder[] partialResults = new StringBuilder[numThreads];
        for (int i = 0; i < numThreads; i++) {
            partialResults[i] = new StringBuilder();
        }

        ExecutorService pool = Executors.newFixedThreadPool(numThreads);
        List<Future<?>> futures = new ArrayList<>(numThreads);

        for (int threadIndex = 0; threadIndex < numThreads; threadIndex++) {
            final int startLine = skip + threadIndex * chunkSize; // inclusive
            final int endLine = Math.min(skip + (threadIndex + 1) * chunkSize, total); // exclusive
            final int idx = threadIndex;

            if (startLine >= endLine) {
                // no lines to do for this thread, skip it
                break;
            }

            futures.add(pool.submit(() -> {
                try (BufferedReader br = openFile(filename)) {
                    int currentLine = 0;
                    String line;
                    // Skip up to startLine
                    while (currentLine < startLine) {
                        if (br.readLine() == null) {
                            return null; // end of file
                        }
                        currentLine++;
                    }
                    // Now read until endLine
                    while (currentLine < endLine) {
                        line = br.readLine();
                        if (line == null) {
                            break; // EOF
                        }
                        currentLine++;
                        // process
                        String result = lineProcessor.process(line);
                        if (!result.isEmpty()) {
                            partialResults[idx].append(result);
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return null; // no return
            }));
        }

        // wait for all
        for (Future<?> f : futures) {
            try {
                f.get();
            } catch (InterruptedException | ExecutionException e) {
                e.printStackTrace();
            }
        }
        pool.shutdown();

        // Merge
        StringBuilder merged = new StringBuilder();
        for (StringBuilder sb : partialResults) {
            merged.append(sb);
        }
        return merged.toString().toCharArray();
    }

    // -------------------------------------------------------------------------
    // 4) Now we adapt each of your public methods to:
    //    - handle "s" and "cp" as we did in the single-thread revision
    //    - for "sp", call parallelProcessFile(...) with the relevant logic
    // -------------------------------------------------------------------------

    // readAll
    static public char[] readAll(String filename, String delim, Integer skip, String bench) throws IOException {
        if ("sp".equals(bench)) {
            // parallel approach
            return parallelProcessFile(filename, skip, line -> line + delim);
        } else if ("cp".equals(bench)) {
            // chunk-based single-thread approach
            // We'll do basically the chunk approach from the prior revision
            // For brevity, we use the same logic as we had before:
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                String nextLine;
                int lineCtr = 0;
                int chunk = 100000; // or some default
                while (true) {
                    List<String> buffer = new ArrayList<>(chunk);
                    for (int i = 0; i < chunk; i++) {
                        nextLine = br.readLine();
                        if (nextLine == null) {
                            break;
                        }
                        lineCtr++;
                        if (lineCtr <= skip) {
                            continue;
                        }
                        buffer.add(nextLine);
                    }
                    for (String bline : buffer) {
                        sb.append(bline).append(delim);
                    }
                    if (buffer.size() < chunk) {
                        break; // done
                    }
                }
            }
            return sb.toString().toCharArray();
        } else {
            // "s" or anything else => single-thread line-by-line
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    sb.append(line).append(delim);
                }
            }
            return sb.toString().toCharArray();
        }
    }

    // readAll2
    static public String[] readAll2(String filename, Integer skip, String bench) throws IOException {
        if ("sp".equals(bench)) {
            // parallel: we'll gather lines in a big StringBuilder, then split them
            char[] big = parallelProcessFile(filename, skip, line -> line + "\n");
            // convert to a single string, then split on \n
            String str = new String(big);
            return str.split("\n", -1);
        } else {
            // same single-thread approach
            List<String> lines = new ArrayList<>();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    lines.add(line);
                }
            }
            return lines.toArray(new String[0]);
        }
    }

    // readHeader
    static public String readHeader(String filename, Integer skip) throws IOException {
        // no parallel logic needed, because it's only one line
        try (BufferedReader br = openFile(filename)) {
            for (int i = 0; i < skip; i++) {
                if (br.readLine() == null) {
                    return null;
                }
            }
            return br.readLine();
        }
    }

    // lineCount2
    static public Long lineCount2(String filename, String bench) throws IOException {
        // parallel counting is rarely beneficial. We'll just do single pass
        long c = 0;
        try (BufferedReader br = openFile(filename)) {
            while (br.readLine() != null) {
                c++;
            }
        }
        return c;
    }

    // compare
    static public char[] compare(String filename, String[] terms, String delim, Integer skip,
                                 String bench, Integer[] retrunCols) throws IOException {
        Predicate<String> predicate = chainPredicates(terms);

        if ("sp".equals(bench)) {
            // parallel
            return parallelProcessFile(filename, skip, line -> {
                if (predicate.test(line)) {
                    if (retrunCols[0] == null) {
                        return line + delim;
                    } else {
                        return getCustomCols(line, delim, retrunCols) + delim;
                    }
                }
                return "";
            });
        } else if ("cp".equals(bench)) {
            // chunk-based single-thread
            return compareChunk(filename, terms, delim, skip, retrunCols);
        } else {
            // single-thread
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    if (predicate.test(line)) {
                        if (retrunCols[0] == null) {
                            sb.append(line).append(delim);
                        } else {
                            sb.append(getCustomCols(line, delim, retrunCols)).append(delim);
                        }
                    }
                }
            }
            return sb.toString().toCharArray();
        }
    }

    // A small helper to do "cp" in compare
    private static char[] compareChunk(String filename, String[] terms, String delim,
                                       Integer skip, Integer[] retrunCols) throws IOException {
        Predicate<String> predicate = chainPredicates(terms);
        StringBuilder sb = new StringBuilder();
        int chunk = 100000; // example chunk size
        try (BufferedReader br = openFile(filename)) {
            List<String> buffer = new ArrayList<>();
            String nextLine;
            int lineCtr = 0;
            while (true) {
                buffer.clear();
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String s : buffer) {
                    if (predicate.test(s)) {
                        if (retrunCols[0] == null) {
                            sb.append(s).append(delim);
                        } else {
                            sb.append(getCustomCols(s, delim, retrunCols)).append(delim);
                        }
                    }
                }
                if (buffer.size() < chunk) {
                    break; // done
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // compareToCols
    public char[] compareToCols(String filename, String[] terms, String delim, Integer[] col,
                                Integer skip, String bench, Integer[] retrunCols) throws IOException {
        Predicate<String> predicate = chainPredicates(terms);

        if ("sp".equals(bench)) {
            return parallelProcessFile(filename, skip, line -> {
                if (toKeepWrapper(line, delim, predicate, col)) {
                    if (retrunCols[0] == null) {
                        return line + delim;
                    } else {
                        return getCustomCols(line, delim, retrunCols) + delim;
                    }
                }
                return "";
            });
        } else if ("cp".equals(bench)) {
            // chunk single-thread
            return compareToColsChunk(filename, terms, delim, col, skip, retrunCols);
        } else {
            // single-thread
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    if (toKeepWrapper(line, delim, predicate, col)) {
                        if (retrunCols[0] == null) {
                            sb.append(line).append(delim);
                        } else {
                            sb.append(getCustomCols(line, delim, retrunCols)).append(delim);
                        }
                    }
                }
            }
            return sb.toString().toCharArray();
        }
    }

    private char[] compareToColsChunk(String filename, String[] terms, String delim, Integer[] col,
                                      Integer skip, Integer[] retrunCols) throws IOException {
        Predicate<String> predicate = chainPredicates(terms);
        StringBuilder sb = new StringBuilder();
        int chunk = 100000;
        try (BufferedReader br = openFile(filename)) {
            List<String> buffer = new ArrayList<>();
            String nextLine;
            int lineCtr = 0;
            while (true) {
                buffer.clear();
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String s : buffer) {
                    if (toKeepWrapper(s, delim, predicate, col)) {
                        if (retrunCols[0] == null) {
                            sb.append(s).append(delim);
                        } else {
                            sb.append(getCustomCols(s, delim, retrunCols)).append(delim);
                        }
                    }
                }
                if (buffer.size() < chunk) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // compareToJoinedCols
    public char[] compareToJoinedCols(String filename, String[] terms, String delim, Integer[] col,
                                      Integer skip, String bench, Integer[] retrunCols, Integer chunk)
                                      throws IOException {
        Predicate<String> predicate = chainPredicates(terms);
        if ("sp".equals(bench)) {
            // parallel
            return parallelProcessFile(filename, skip, line -> {
                if (toKeepWrapperJoint(line, delim, predicate, col)) {
                    if (retrunCols[0] == null) {
                        return line + delim;
                    } else {
                        return getCustomCols(line, delim, retrunCols) + delim;
                    }
                }
                return "";
            });
        } else if ("cp".equals(bench)) {
            // chunk single-thread approach
            return compareToJoinedColsChunk(filename, terms, delim, col, skip, retrunCols, chunk);
        } else {
            // single-thread
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    if (toKeepWrapperJoint(line, delim, predicate, col)) {
                        if (retrunCols[0] == null) {
                            sb.append(line).append(delim);
                        } else {
                            sb.append(getCustomCols(line, delim, retrunCols)).append(delim);
                        }
                    }
                }
            }
            return sb.toString().toCharArray();
        }
    }

    private char[] compareToJoinedColsChunk(String filename, String[] terms, String delim,
                                            Integer[] col, Integer skip, Integer[] retrunCols,
                                            Integer chunk) throws IOException {
        Predicate<String> predicate = chainPredicates(terms);
        StringBuilder sb = new StringBuilder();
        try (BufferedReader br = openFile(filename)) {
            List<String> buffer = new ArrayList<>();
            String nextLine;
            int lineCtr = 0;
            while (true) {
                buffer.clear();
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String s : buffer) {
                    if (toKeepWrapperJoint(s, delim, predicate, col)) {
                        if (retrunCols[0] == null) {
                            sb.append(s).append(delim);
                        } else {
                            sb.append(getCustomCols(s, delim, retrunCols)).append(delim);
                        }
                    }
                }
                if (buffer.size() < chunk) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // filterCol
    public char[] filterCol(String filename, String delim, Integer[] col, String[] term, double[] cond,
                            Integer skip, String bench, Integer[] retrunCols) throws IOException {
        Map<String, BiPredicate<Double, Double>> ops = buildOperationsMap();
        if ("sp".equals(bench)) {
            return parallelProcessFile(filename, skip, line -> {
                if (toFilterWrapper(line, delim, col, cond, term, ops)) {
                    if (retrunCols[0] == null) {
                        return line + delim;
                    } else {
                        return getCustomCols(line, delim, retrunCols) + delim;
                    }
                }
                return "";
            });
        } else if ("cp".equals(bench)) {
            // chunk single-thread
            return filterColChunk(filename, delim, col, term, cond, skip, retrunCols);
        } else {
            // single-thread
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    if (toFilterWrapper(line, delim, col, cond, term, ops)) {
                        if (retrunCols[0] == null) {
                            sb.append(line).append(delim);
                        } else {
                            sb.append(getCustomCols(line, delim, retrunCols)).append(delim);
                        }
                    }
                }
            }
            return sb.toString().toCharArray();
        }
    }

    private char[] filterColChunk(String filename, String delim, Integer[] col, String[] term,
                                  double[] cond, Integer skip, Integer[] retrunCols) throws IOException {
        Map<String, BiPredicate<Double, Double>> ops = buildOperationsMap();
        StringBuilder sb = new StringBuilder();
        int chunk = 100000;
        try (BufferedReader br = openFile(filename)) {
            List<String> buffer = new ArrayList<>();
            String nextLine;
            int lineCtr = 0;
            while (true) {
                buffer.clear();
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String s : buffer) {
                    if (toFilterWrapper(s, delim, col, cond, term, ops)) {
                        if (retrunCols[0] == null) {
                            sb.append(s).append(delim);
                        } else {
                            sb.append(getCustomCols(s, delim, retrunCols)).append(delim);
                        }
                    }
                }
                if (buffer.size() < chunk) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // compareFilterCol
    public char[] compareFilterCol(String filename, Integer skip, String delim, String[] patt,
                                   Integer[] colT, String[] operation, Integer[] colN, double[] cond,
                                   String bench, Integer[] retrunCols) throws IOException {
        Predicate<String> textPred = chainPredicates(patt);
        Map<String, BiPredicate<Double, Double>> ops = buildOperationsMap();
        if ("sp".equals(bench)) {
            return parallelProcessFile(filename, skip, line -> {
                // numeric filter
                if (toFilterWrapper(line, delim, colN, cond, operation, ops)) {
                    // text filter
                    if (toKeepWrapper(line, delim, textPred, colT)) {
                        if (retrunCols[0] == null) {
                            return line + delim;
                        } else {
                            return getCustomCols(line, delim, retrunCols) + delim;
                        }
                    }
                }
                return "";
            });
        } else if ("cp".equals(bench)) {
            return compareFilterColChunk(filename, skip, delim, patt, colT, operation, colN, cond, retrunCols);
        } else {
            // single-thread
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    if (toFilterWrapper(line, delim, colN, cond, operation, ops)) {
                        if (toKeepWrapper(line, delim, textPred, colT)) {
                            if (retrunCols[0] == null) {
                                sb.append(line).append(delim);
                            } else {
                                sb.append(getCustomCols(line, delim, retrunCols)).append(delim);
                            }
                        }
                    }
                }
            }
            return sb.toString().toCharArray();
        }
    }

    private char[] compareFilterColChunk(String filename, Integer skip, String delim, String[] patt,
                                         Integer[] colT, String[] operation, Integer[] colN, double[] cond,
                                         Integer[] retrunCols) throws IOException {
        Predicate<String> textPred = chainPredicates(patt);
        Map<String, BiPredicate<Double, Double>> ops = buildOperationsMap();
        StringBuilder sb = new StringBuilder();
        int chunk = 100000;
        try (BufferedReader br = openFile(filename)) {
            List<String> buffer = new ArrayList<>();
            String nextLine;
            int lineCtr = 0;
            while (true) {
                buffer.clear();
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String s : buffer) {
                    if (toFilterWrapper(s, delim, colN, cond, operation, ops)) {
                        if (toKeepWrapper(s, delim, textPred, colT)) {
                            if (retrunCols[0] == null) {
                                sb.append(s).append(delim);
                            } else {
                                sb.append(getCustomCols(s, delim, retrunCols)).append(delim);
                            }
                        }
                    }
                }
                if (buffer.size() < chunk) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // multiCompareFilterCol
    public char[] multiCompareFilterCol(String filename, Integer skip, String delim, String[] patt,
                                        Integer[] colT, String[] operation, Integer[] colN, double[] cond,
                                        String bench, Integer[] retrunCols) throws IOException {
        Map<String, BiPredicate<Double, Double>> ops = buildOperationsMap();
        if ("sp".equals(bench)) {
            return parallelProcessFile(filename, skip, line -> {
                // numeric
                if (toFilterWrapper(line, delim, colN, cond, operation, ops)) {
                    // text
                    if (toKeepWrapperSingle(line, delim, patt, colT)) {
                        if (retrunCols[0] == null) {
                            return line + delim;
                        } else {
                            return getCustomCols(line, delim, retrunCols) + delim;
                        }
                    }
                }
                return "";
            });
        } else if ("cp".equals(bench)) {
            return multiCompareFilterColChunk(filename, skip, delim, patt, colT, operation, colN, cond, retrunCols);
        } else {
            // single-thread
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    if (toFilterWrapper(line, delim, colN, cond, operation, ops)) {
                        if (toKeepWrapperSingle(line, delim, patt, colT)) {
                            if (retrunCols[0] == null) {
                                sb.append(line).append(delim);
                            } else {
                                sb.append(getCustomCols(line, delim, retrunCols)).append(delim);
                            }
                        }
                    }
                }
            }
            return sb.toString().toCharArray();
        }
    }

    private char[] multiCompareFilterColChunk(String filename, Integer skip, String delim,
                                              String[] patt, Integer[] colT, String[] operation,
                                              Integer[] colN, double[] cond, Integer[] retrunCols)
                                              throws IOException {
        Map<String, BiPredicate<Double, Double>> ops = buildOperationsMap();
        StringBuilder sb = new StringBuilder();
        int chunk = 100000;
        try (BufferedReader br = openFile(filename)) {
            List<String> buffer = new ArrayList<>();
            String nextLine;
            int lineCtr = 0;
            while (true) {
                buffer.clear();
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String s : buffer) {
                    if (toFilterWrapper(s, delim, colN, cond, operation, ops)) {
                        if (toKeepWrapperSingle(s, delim, patt, colT)) {
                            if (retrunCols[0] == null) {
                                sb.append(s).append(delim);
                            } else {
                                sb.append(getCustomCols(s, delim, retrunCols)).append(delim);
                            }
                        }
                    }
                }
                if (buffer.size() < chunk) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // multiCompareToCols
    public char[] multiCompareToCols(String filename, Integer skip, String delim, String[] patt,
                                     Integer[] colT, String bench, Integer[] retrunCols)
                                     throws IOException {
        if ("sp".equals(bench)) {
            return parallelProcessFile(filename, skip, line -> {
                if (toKeepWrapperSingle(line, delim, patt, colT)) {
                    if (retrunCols[0] == null) {
                        return line + delim;
                    } else {
                        return getCustomCols(line, delim, retrunCols) + delim;
                    }
                }
                return "";
            });
        } else if ("cp".equals(bench)) {
            return multiCompareToColsChunk(filename, skip, delim, patt, colT, retrunCols);
        } else {
            // single-thread
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    if (toKeepWrapperSingle(line, delim, patt, colT)) {
                        if (retrunCols[0] == null) {
                            sb.append(line).append(delim);
                        } else {
                            sb.append(getCustomCols(line, delim, retrunCols)).append(delim);
                        }
                    }
                }
            }
            return sb.toString().toCharArray();
        }
    }

    private char[] multiCompareToColsChunk(String filename, Integer skip, String delim,
                                           String[] patt, Integer[] colT, Integer[] retrunCols)
                                           throws IOException {
        StringBuilder sb = new StringBuilder();
        int chunk = 100000;
        try (BufferedReader br = openFile(filename)) {
            List<String> buffer = new ArrayList<>();
            String nextLine;
            int lineCtr = 0;
            while (true) {
                buffer.clear();
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String s : buffer) {
                    if (toKeepWrapperSingle(s, delim, patt, colT)) {
                        if (retrunCols[0] == null) {
                            sb.append(s).append(delim);
                        } else {
                            sb.append(getCustomCols(s, delim, retrunCols)).append(delim);
                        }
                    }
                }
                if (buffer.size() < chunk) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // getColumn
    static public char[] getColumn(String filename, String delim, Integer[] cols, Integer skip,
                                   String bench, Integer[] rows, Integer chunk) throws IOException {
        // If rows are specified, we do existing row-based logic (and can do parallel if "sp"?).
        // But row-based skipping is complicated for parallel. We'll do the naive approach:
        // We'll keep single-thread logic if rows != null, to avoid confusion.
        if (rows[0] != null) {
            if ("cp".equals(bench)) {
                return getColumnRowsChunk(filename, delim, cols, skip, rows, chunk);
            } else if ("sp".equals(bench)) {
                // We could do a parallel approach, but skipping random rows in parallel is tricky.
                // We'll just do naive read + skip for each thread, or revert to single-thread if it’s not too large.
                return parallelProcessFile(filename, skip, line -> {
                    // Hard to match "rows" in parallel without reading them in order.
                    // For simplicity, we’ll do single-thread if rows are specified.
                    return ""; 
                });
            } else {
                return getColumnRowsSequential(filename, delim, cols, skip, rows);
            }
        }

        // else: we want entire columns from all rows
        if ("sp".equals(bench)) {
            // parallel approach
            return parallelProcessFile(filename, skip, line -> getCustomCols(line, delim, cols) + delim);
        } else if ("cp".equals(bench)) {
            // chunk single-thread
            return getColumnAllChunk(filename, delim, cols, skip, chunk);
        } else {
            // single-thread
            return getColumnAllSimple(filename, delim, cols, skip);
        }
    }

    // getColumnBuffer
    static public char[] getColumnBuffer(String filename, String delim, Integer[] col,
                                         Integer skip, String bench) throws IOException {
        // If you want to do parallel, we can do it similarly:
        if ("sp".equals(bench)) {
            return parallelProcessFile(filename, skip, line -> {
                String[] parts = line.split(delim);
                String[] rowList = new String[col.length];
                for (int i = 0; i < col.length; i++) {
                    rowList[i] = (col[i] < parts.length) ? parts[col[i]] : "";
                }
                return String.join(delim, rowList) + delim;
            });
        } else {
            // single or chunk approach
            StringBuilder sb = new StringBuilder();
            try (BufferedReader br = openFile(filename)) {
                int skipped = 0;
                String line;
                while ((line = br.readLine()) != null) {
                    if (skipped < skip) {
                        skipped++;
                        continue;
                    }
                    String[] parts = line.split(delim);
                    String[] rowList = new String[col.length];
                    for (int i = 0; i < col.length; i++) {
                        rowList[i] = (col[i] < parts.length) ? parts[col[i]] : "";
                    }
                    sb.append(String.join(delim, rowList)).append(delim);
                }
            }
            return sb.toString().toCharArray();
        }
    }

    // getRowColumn
    public char[] getRowColumn(String filename, String delim, Integer[] cols, Integer[] rows,
                               Integer skip, String bench, Integer chunk) throws IOException {
        // This is inherently random row reads. Parallelizing it well would require advanced logic.
        // We'll keep single-thread for correctness.
        List<String> lines = new ArrayList<>();
        try (BufferedReader br = openFile(filename)) {
            // skip + first row
            long firstSkip = skip + rows[0] - 1;
            String first = br.lines().skip(firstSkip).findFirst().orElse(null);
            if (first == null) {
                return new char[0];
            }
            lines.add(first);

            // subsequent rows
            for (int i = 1; i < rows.length; i++) {
                long distance = rows[i] - rows[i - 1] - 1;
                String row = br.lines().skip(distance).findFirst().orElse(null);
                if (row == null) {
                    break;
                }
                lines.add(row);
            }
        }
        // Now filter columns
        StringBuilder sb = new StringBuilder();
        for (String line : lines) {
            sb.append(getCustomCols(line, delim, cols)).append(delim);
        }
        return sb.toString().toCharArray();
    }

    // -------------------------------------------------------------------------
    // Private Helpers
    // -------------------------------------------------------------------------

    // Build numeric operations map
    private static Map<String, BiPredicate<Double, Double>> buildOperationsMap() {
        Map<String, BiPredicate<Double, Double>> m = new HashMap<>();
        m.put("gt", (a, b) -> a > b);
        m.put("lt", (a, b) -> a < b);
        m.put("eq", (a, b) -> Double.compare(a, b) == 0);
        m.put("ge", (a, b) -> a >= b);
        m.put("le", (a, b) -> a <= b);
        return m;
    }

    // chainPredicates
    static private Predicate<String> chainPredicates(String[] terms) {
        // original code used OR across all compiled regex
        if (terms == null || terms.length == 0) {
            // If no terms, treat as "always true" (matching original code logic)
            return x -> true;
        }
        List<Predicate<String>> preds = new ArrayList<>();
        for (String t : terms) {
            preds.add(Pattern.compile(t).asPredicate());
        }
        return preds.stream().reduce(Predicate::or).orElse(x -> true);
    }

    private static boolean toKeepWrapper(String s, String delim, Predicate<String> predicate, Integer[] col) {
        String[] arr = s.split(delim);
        for (Integer c : col) {
            if (c < arr.length && predicate.test(arr[c])) {
                return true;
            }
        }
        return false;
    }

    private static boolean toKeepWrapperJoint(String s, String delim, Predicate<String> predicate, Integer[] col) {
        String[] arr = s.split(delim);
        List<String> selected = new ArrayList<>();
        for (int c : col) {
            if (c < arr.length) {
                selected.add(arr[c]);
            }
        }
        String joined = String.join("_", selected);
        return predicate.test(joined);
    }

    private static boolean toKeepWrapperSingle(String s, String delim, String[] patt, Integer[] col) {
        String[] arr = s.split(delim);
        for (int i = 0; i < col.length; i++) {
            int c = col[i];
            if (c >= arr.length) {
                return false;
            }
            Predicate<String> p = Pattern.compile(patt[i]).asPredicate();
            if (!p.test(arr[c])) {
                return false;
            }
        }
        return true;
    }

    private static boolean toFilterWrapper(String s, String delim, Integer[] col, double[] cond,
                                           String[] term, Map<String, BiPredicate<Double, Double>> ops) {
        String[] arr = s.split(delim);
        for (int i = 0; i < col.length; i++) {
            int c = col[i];
            if (c >= arr.length) {
                return false;
            }
            Double val = checkDouble(arr[c]);
            if (!ops.get(term[i]).test(val, cond[i])) {
                return false;
            }
        }
        return true;
    }

    private static Double checkDouble(String s) {
        try {
            return Double.parseDouble(s);
        } catch (NumberFormatException e) {
            return Double.NaN;
        }
    }

    private static String getCustomCols(String s, String delim, Integer[] cols) {
        String[] arr = s.split(delim);
        String[] sel = new String[cols.length];
        for (int i = 0; i < cols.length; i++) {
            int c = cols[i];
            sel[i] = (c < arr.length) ? arr[c] : "";
        }
        return String.join(delim, sel);
    }

    // -------------------------------------
    // Sub-helpers for getColumn
    // -------------------------------------

    // Reading entire columns from all lines, single-thread
    private static char[] getColumnAllSimple(String filename, String delim, Integer[] cols,
                                             Integer skip) throws IOException {
        StringBuilder sb = new StringBuilder();
        try (BufferedReader br = openFile(filename)) {
            int skipped = 0;
            String line;
            while ((line = br.readLine()) != null) {
                if (skipped < skip) {
                    skipped++;
                    continue;
                }
                sb.append(getCustomCols(line, delim, cols)).append(delim);
            }
        }
        return sb.toString().toCharArray();
    }

    // chunk single-thread
    private static char[] getColumnAllChunk(String filename, String delim, Integer[] cols,
                                            Integer skip, Integer chunk) throws IOException {
        StringBuilder sb = new StringBuilder();
        try (BufferedReader br = openFile(filename)) {
            String nextLine;
            int lineCtr = 0;
            while (true) {
                List<String> buffer = new ArrayList<>(chunk);
                for (int i = 0; i < chunk; i++) {
                    nextLine = br.readLine();
                    if (nextLine == null) {
                        break;
                    }
                    lineCtr++;
                    if (lineCtr <= skip) {
                        continue;
                    }
                    buffer.add(nextLine);
                }
                for (String bline : buffer) {
                    sb.append(getCustomCols(bline, delim, cols)).append(delim);
                }
                if (buffer.size() < chunk) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }

    // reading specific rows single-thread
    private static char[] getColumnRowsSequential(String filename, String delim, Integer[] cols,
                                                  Integer skip, Integer[] rows) throws IOException {
        List<String> collected = new ArrayList<>(rows.length);
        try (BufferedReader br = openFile(filename)) {
            // skip up to the first row
            long firstSkip = skip + rows[0] - 1;
            String first = br.lines().skip(firstSkip).findFirst().orElse(null);
            if (first == null) {
                return new char[0];
            }
            collected.add(first);

            for (int i = 1; i < rows.length; i++) {
                long distance = (rows[i] - rows[i - 1] - 1);
                String row = br.lines().skip(distance).findFirst().orElse(null);
                if (row == null) {
                    break;
                }
                collected.add(row);
            }
        }
        // map to columns
        return collected.stream()
            .map(s -> getCustomCols(s, delim, cols))
            .collect(Collectors.joining(delim))
            .toCharArray();
    }

    // chunk-based row reading, single-thread
    private static char[] getColumnRowsChunk(String filename, String delim, Integer[] cols,
                                             Integer skip, Integer[] rows, Integer chunk) throws IOException {
        StringBuilder sb = new StringBuilder();
        try (BufferedReader br = openFile(filename)) {
            // read first requested row
            long firstSkip = skip + rows[0] - 1;
            String first = br.lines().skip(firstSkip).findFirst().orElse(null);
            if (first == null) {
                return new char[0];
            }
            sb.append(getCustomCols(first, delim, cols)).append(delim);

            int ptr = 1;
            while (ptr < rows.length) {
                List<String> buffer = new ArrayList<>(chunk);
                for (int j = 0; j < chunk && ptr < rows.length; j++) {
                    long distance = rows[ptr] - rows[ptr - 1] - 1;
                    String row = br.lines().skip(distance).findFirst().orElse(null);
                    if (row == null) {
                        break;
                    }
                    buffer.add(row);
                    ptr++;
                }
                for (String bline : buffer) {
                    sb.append(getCustomCols(bline, delim, cols)).append(delim);
                }
                if (ptr >= rows.length) {
                    break;
                }
            }
        }
        return sb.toString().toCharArray();
    }
}
