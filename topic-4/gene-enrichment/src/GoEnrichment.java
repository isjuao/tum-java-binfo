import org.apache.commons.math3.distribution.HypergeometricDistribution;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class GoEnrichment {
    public static File oboFile = null;
    public static String mapFileStr = null;
    public static File enrichFile = null;
    public static boolean mapTypeGo;
    public static char root;
    public static GoTerm rootNode;
    public static int N;
    public static int K;
    public static int minSize;
    public static int maxSize;

    public static HashMap<String, GoTerm> dag = new HashMap<>();
    public static HashMap<String, Gene> allGenes = new HashMap<>();
    public static HashSet<GoTerm> trues = new HashSet<>();

    /**
     * Main procedure.
     * @param outputTsvString
     * @param overlapOutTsvString
     * @throws IOException
     */
    private static void performEnrichment(String outputTsvString, String overlapOutTsvString)
        throws IOException {
        // 1) read the obo file into a HashMap
        collectTerms();

        // 2) build the DAG by adding pointers by using the HashMap
        buildDag();

        // 3) assign genes to associated GO nodes
        geneMapping();

        // 4) add enrichment info to genes and annotate true nodes
        addEnrichment();
        int x = 6;  // breakpoint

        // 5) calculate values for statistics
        int compN = 0;
        int compK = 0;
        for(Gene gene : allGenes.values()) {
            if(gene.significant != null) {
                compN++;
                if(gene.significant) {
                    compK++;
                }
            }
        }
        N = compN;
        K = compK;

        int discardedNodes = calculateStatistics();
        mtCorrection(discardedNodes);

        if(!overlapOutTsvString.equals("")) {
            // 6) overlap output
            HashSet<String> overlapping = calcOverlapPairs();
            writeOverlapOutput(overlapping, overlapOutTsvString);
        }

        // 7) info output
        writeOutput(outputTsvString, minSize, maxSize);
    }

    /**
     * Writes the whole output (not only one line!) for the enrichment information file.
     * @param outputTsvString
     * @param minSize
     * @param maxSize
     * @throws IOException
     */
    private static void writeOutput(String outputTsvString, int minSize, int maxSize) throws IOException {
        File out = new File(outputTsvString);
        out.delete();
        FileOutputStream outStream = new FileOutputStream(out);
        BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(outStream));
        String outLine;
        // write header
        bwOut.write("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\t" +
                        "ks_fdr\tshortest_path_to_a_true");
        bwOut.newLine();
        // write file content
        for(GoTerm node : dag.values()) {
            HashSet<Gene> genesOfNode = node.geneSet;
            int size = genesOfNode.size();
            if(size < minSize || size > maxSize) continue;
            int[] geneSizes = node.getMeasuredSize();
            int measuredSize = geneSizes[0];
            int nOverlap = geneSizes[1];
            String pathString = shortestPathToTrue(node);

            // write
            outLine = "" + node.id + "\t" + node.name + "\t" + measuredSize + "\t" + node.simulTrue + "\t" + nOverlap
                    +"\t" + node.hgP + "\t" + node.hgPCorr + "\t" + node.fejP + "\t" + node.fejPCorr + "\t" + node.kStat
                    +"\t" + node.kP + "\t" + node.kPCorr + "\t" + pathString;
            bwOut.write(outLine);
            bwOut.newLine();
        }
        bwOut.close();
    }

    /**
     * Calculates the exact overlap of given overlapping term pairs and writes the output (one line per pair).
     * @param termPairs
     * @param overlapOutTsvString
     * @throws IOException
     */
    private static void writeOverlapOutput(HashSet<String> termPairs, String overlapOutTsvString) throws IOException {
        File overlapOut = new File(overlapOutTsvString);
        overlapOut.delete();
        FileOutputStream outStream = new FileOutputStream(overlapOut);
        BufferedWriter bwOOut = new BufferedWriter(new OutputStreamWriter(outStream));
        String outLine;
        // write header
        bwOOut.write("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent");
        bwOOut.newLine();

        String firstID;
        String secondID;
        GoTerm first;
        GoTerm second;
        HashSet<String> genesFromBoth;
        for(String pair : termPairs) {
            firstID = pair.substring(0, 10);
            secondID = pair.substring(11);
            if(firstID.equals("GO:0006487") && secondID.equals("GO:0009100")) {
                int z = 2;  // breakpoint
            }

            first = dag.get(firstID);
            second = dag.get(secondID);
            // check if second over first
            /*boolean related = first.isDescendantOf(second);
            if(!related) {
                // check if first over second
                related = second.isDescendantOf(first);
            }*/
            // check if second over first
            boolean related = true;
            int pathLength = first.getPathLengthToAscendant(second, 0);
            if(pathLength == -1) {
                // check if first over second
                pathLength = second.getPathLengthToAscendant(first, 0);
            }
            if(pathLength == -1) {
                related = false;
                // have to find shortest path over lca
                GoTerm lca = getLCA(first, second);
                int firstPathLength = first.getPathLengthToAscendant(lca, 0);
                int secondPathLength = second.getPathLengthToAscendant(lca, 0);
                pathLength = firstPathLength + secondPathLength;
            }

            // calculate overlap
            int numOverlapping = 0;
            for(Gene g : second.geneSet) {
                if(first.geneSet.contains(g)) {
                    numOverlapping++;
                }
            }

            /*genesFromBoth = new HashSet<>();
            for(Gene g : first.geneSet) {
                genesFromBoth.add(g.id);
            }
            for(Gene g : second.geneSet) {
                id = g.id;
                if(!genesFromBoth.contains(id)) {
                    // non-overlapping
                    genesFromBoth.add(id);
                } else {
                    // overlapping
                    numOverlapping++;
                }
            }
            float maxOvPercent = (float)numOverlapping/(float)(numOverlapping + genesFromBoth.size());*/
            int firstSize = first.geneSet.size();
            int secondSize = second.geneSet.size();
            float maxOvPercent = 100*((float)numOverlapping/(float)Math.min(firstSize, secondSize));

            outLine = firstID + "\t" + secondID + "\t" + related + "\t" + pathLength + "\t" +
                    numOverlapping + "\t" + String.format("%.02f", maxOvPercent);
            bwOOut.write(outLine);
            bwOOut.newLine();
        }
        bwOOut.close();
    }

    /**
     * Finds the overlapping GoTerms in the dag by going through all genes.
     * @return a HashSet of all pairs in the format of "*smallerID*.*biggerID*"
     */
    private static HashSet<String> calcOverlapPairs() {
        HashSet<String> termPairs = new HashSet<>();
        HashSet<GoTerm> mappedTerms;
        String firstID;
        String secondID;
        for(Gene gene : allGenes.values()) {
            mappedTerms = gene.belongsTo;

            term: for(GoTerm mappedTerm1 : mappedTerms) {
                firstID = mappedTerm1.id;
                int size1 = mappedTerm1.geneSet.size();
                if(size1 < minSize || size1 > maxSize) continue term;

                partner: for(GoTerm mappedTerm2 : mappedTerms) {
                    secondID = mappedTerm2.id;
                    int size2 = mappedTerm2.geneSet.size();
                    if(size2 < minSize || size2 > maxSize) continue partner;

                    int order = firstID.compareTo(secondID);
                    if(order == 0) continue partner;    // skip equal case

                    String toAdd;
                    if(order < 0) {
                        toAdd = firstID + "." + secondID;
                    } else {
                        toAdd = secondID + "." + firstID;
                    }
                    if(!termPairs.contains(toAdd)) {
                        termPairs.add(toAdd);
                        //System.out.println(toAdd);
                    } else {
                        int z = 5;
                    }
                }
            }
        }
        return termPairs;
    }

    /**
     * Calculates the LCA for 2 given GoTerms.
     * @param first
     * @param second
     * @return the LCA GoTerm
     */
    public static GoTerm getLCA(GoTerm first, GoTerm second) {
        HashMap<String, Integer> height1 = new HashMap<>();
        HashMap<String, Integer> height2 = new HashMap<>();
        HashMap<String, Integer> possibleLCAs = new HashMap<>();
        HashSet<String> ancestors1 = new HashSet<>();
        HashSet<String> ancestors2 = new HashSet<>();
        HashSet<GoTerm> lastParents1 = new HashSet<>();
        HashSet<GoTerm> lastParents2 = new HashSet<>();
        HashSet<GoTerm> nextLastParents1 = new HashSet<>();
        HashSet<GoTerm> nextLastParents2 = new HashSet<>();
        String parentID = "";
        String lcaID = "";

        lastParents1.add(first);
        lastParents2.add(second);
        int currentHeight = 0;
        outer: while(lastParents1.size() != 0 || lastParents2.size() != 0) {
            currentHeight++;
            // node
            for(GoTerm current1 : lastParents1) {
                for(GoTerm parent : current1.parents) {
                    parentID = parent.id;
                    /*if(ancestors2.contains(parentID)) {
                        // found LCA???
                        lcaID = parentID;
                        break outer;
                    } else {
                        ancestors1.add(parent.id);
                        nextLastParents1.add(parent);
                    }*/
                    if(height2.containsKey(parentID)) {
                        // found POSSIBLE LCA
                        possibleLCAs.put(parentID, height2.get(parentID));
                    } else {
                        if(!height1.containsKey(parentID)) {
                            height1.put(parentID, currentHeight);
                        }
                        nextLastParents1.add(parent);
                    }
                }
            }
            lastParents1.clear();
            lastParents1.addAll(nextLastParents1);
            nextLastParents1.clear();

            // node 2
            for(GoTerm current2 : lastParents2) {
                for(GoTerm parent : current2.parents) {
                    parentID = parent.id;
                    /*if(ancestors1.contains(parentID)) {
                        // found LCA???
                        lcaID = parentID;
                        break outer;
                    } else {
                        ancestors2.add(parent.id);
                        nextLastParents2.add(parent);
                    }*/
                    if(height1.containsKey(parentID)) {
                        // found POSSIBLE LCA
                        possibleLCAs.put(parentID, height1.get(parentID));
                    } else {
                        if(!height2.containsKey(parentID)) {
                            height2.put(parentID, currentHeight);
                        }
                        nextLastParents2.add(parent);
                    }
                }
            }
            lastParents2.clear();
            lastParents2.addAll(nextLastParents2);
            nextLastParents2.clear();

            // look for LCAs
            if(!possibleLCAs.isEmpty()) {
                lcaID = "";
                int minHeight = Integer.MAX_VALUE;
                for(String possLcaID : possibleLCAs.keySet()) {
                    int height = possibleLCAs.get(possLcaID);
                    if(height < minHeight) {
                        lcaID = possLcaID;
                        minHeight = height;
                    }
                }
                break outer;
            }

            // update ancestors1
            /*for(GoTerm newAncestor : nextLastParents1) {
                ancestors1.add(newAncestor.id);
            }
            nextLastParents1.clear();*/
        }
        GoTerm lca;
        if(!lcaID.equals("")) {
            // found lca early
            lca = dag.get(lcaID);
        } else {
            // lca is root
            lca = rootNode;
        }
        return lca;
    }
    /*
    public static String getLCAID(GoTerm first, GoTerm second) {
        HashSet<String> ancestors1 = new HashSet<>();
        HashSet<String> ancestors2 = new HashSet<>();
        HashSet<GoTerm> lastParents1 = new HashSet<>();
        HashSet<GoTerm> lastParents2 = new HashSet<>();
        HashSet<GoTerm> nextLastParents1 = new HashSet<>();
        HashSet<GoTerm> nextLastParents2 = new HashSet<>();
        String parentID = "";
        String lcaID = "";

        lastParents1.add(first);
        lastParents2.add(second);
        outer: while(lastParents1.size() != 0 || lastParents2.size() != 0) {
            // node
            for(GoTerm current1 : lastParents1) {
                for(GoTerm parent : current1.parents) {
                    parentID = parent.id;
                    if(ancestors2.contains(parentID)) {
                        // found LCA???
                        lcaID = parentID;
                        break outer;
                    } else {
                        ancestors1.add(parent.id);
                        nextLastParents1.add(parent);
                    }
                }
            }
            lastParents1.clear();
            lastParents1.addAll(nextLastParents1);
            nextLastParents1.clear();
            // node 2
            for(GoTerm current2 : lastParents2) {
                for(GoTerm parent : current2.parents) {
                    parentID = parent.id;
                    if(ancestors1.contains(parentID)) {
                        // found LCA???
                        lcaID = parentID;
                        break outer;
                    } else {
                        ancestors2.add(parent.id);
                        nextLastParents2.add(parent);
                    }
                }
            }
            lastParents2.clear();
            lastParents2.addAll(nextLastParents2);
            nextLastParents2.clear();
        }
        return lcaID;
    }*/

    /**
     * Calculates the shortest path between 2 given nodes.
     * @param first
     * @param second
     * @return the direct path, if directly related, or the path over the least common ancestor.
     */
    private static String shortestPathBetweenTwo(GoTerm first, GoTerm second) {
        String path = first.getPathToAscendant3(second, "");
        if(path != null) {
            // directly related 1
            path = idToName(path);
            if(path.startsWith("|")) {
                path = path.substring(1);
            }
            path += " *";   // ends with LCA
        } else {
            // check if first over second
             path = second.getPathToAscendant3(first, "");
            if(path != null) {
                // directly related 2
                path = idToName(path);
                path += " * ";
                path = reversePath(path, 1);
            }
        }
        if(path == null) {
            // have to find shortest path over lca
            GoTerm lca = getLCA(first, second);
            String firstPath = first.getPathToAscendant3(lca, "");
            String secondPath = second.getPathToAscendant3(lca, "");

            firstPath = idToName(firstPath);    // last one is lca
            if(firstPath.startsWith("|")) {
                firstPath = firstPath.substring(1);
            }
            secondPath = idToName(secondPath);  // last one is lca
            String reverseSecondPath = reversePath(secondPath, 2);

            path = firstPath + " * |" + reverseSecondPath;
        }
        return path;
    }

    /**
     * Calculates the shortest path to a simulated-true entry in the DAG.
     * @param query
     * @return path in String format, starting at query and ending in closest simulated-true node
     */
    private static String shortestPathToTrue(GoTerm query) {
        if(query.id.equals("GO:0051224")) {
            int x = 4;  // breakpoint
        }

        String path = "";
        long lengthOld = Integer.MAX_VALUE;
        if(query.simulTrue || trues.isEmpty()) {
            return "";
        }
        for(GoTerm trueTerm : trues) {
            String newPath = shortestPathBetweenTwo(query, trueTerm);
            // calculate length of new path
            long lengthNew = newPath.chars().filter(counter -> counter == '|').count();
            if(path.equals("") || lengthNew < lengthOld) {
                path = newPath;
                lengthOld = lengthNew;
            }
        }
        return path;
    }

    /**
     * Calculates the reverse of a given path, ending at the LCA.
     * @param secondPath
     * @return the reversed (now correct) path, starting with the node after the LCA
     */
    private static String reversePath(String secondPath, int ifCutLCA2Else1) {
        StringBuilder reverseSecondPathB = new StringBuilder();
        String[] names = secondPath.split("\\|");
        for(int i = names.length - ifCutLCA2Else1; i > 0; i--) {  // i > 0 bc last one should be empty (starts with "|")
            reverseSecondPathB.append(names[i]).append("|");
        }
        String reverseSecondPath = reverseSecondPathB.toString();
        if(reverseSecondPath.endsWith("|")) {
            reverseSecondPath = reverseSecondPath.substring(0, reverseSecondPath.length() - 1);
        }
        return reverseSecondPath;
    }

    /**
     * Collects the trivial names of a path of GoTerm ids.
     * @param path
     * @return the path containing names instead of ids.
     */
    private static String idToName(String path) {
        StringBuilder result = new StringBuilder();
        for(String s : path.split("\\|")) {
            if(!s.isEmpty()) {
                GoTerm term = dag.get(s);
                result.append("|").append(term.name);
            }
        }
        return result.toString();
    }

    /**
     * Performs multiple testing correction for all p-values of all genes.
     * @param discardedNodes number of nodes that were discarded due to empty gene sets / min-/maxSize issues
     */
    private static void mtCorrection(int discardedNodes) {
        int numOfPs = dag.size() - discardedNodes;
        /*double[] allHgPs = new double[numOfPs];
        double[] allFejPs = new double[numOfPs];
        double[] allKsPs = new double[numOfPs];*/
        int i = 0;
        ArrayList<GoTerm> notDiscarded = new ArrayList<>();
        for(GoTerm term : dag.values()) {
            if(!term.discardThis) {
                if(term.id.equals("GO:0006638")) {
                    int x = 3;
                }
                notDiscarded.add(term);
                /*allHgPs[i] = term.hgP;
                allFejPs[i] = term.fejP;
                allKsPs[i] = term.kP;
                i++;*/
            }
        }
        benjaminiHochberg2(notDiscarded, 'h');
        benjaminiHochberg2(notDiscarded, 'f');
        benjaminiHochberg2(notDiscarded, 'k');
        /*i = 0;
        for(GoTerm term : dag.values()) {
            if(!term.discardThis) {
                term.hgPCorr = adjHgPs[i];
                term.fejPCorr = adjFejPs[i];
                term.kPCorr = adjKsPs[i];
                i++;
            }
        }*/
    }

    /*
     * Performs Benjamini-Hochberg correction.
     * @param pValues
     * @return FDR
     */
    /*
    private static double[] benjaminiHochberg(double[] pValues) {
        // double[] unsorted = pValues; Array copy?
        java.util.Arrays.sort(pValues);
        int n = pValues.length;
        double[] fdr = new double[n];
        double prevAdjP = -1.0d;
        double p;
        for(int i = n - 1; i >= 0; i--) {
            if(prevAdjP == -1.0d) {
                p = pValues[i];
            } else {
                p = Math.min( prevAdjP, (pValues[i]*((double)n/(double)(i+1))) );
            }
            fdr[i] = p;
            prevAdjP = p;
        }
        return fdr;
    }*/

    /**
     * Performs Benjamini-Hochberg correction.
     * @param notDiscarded all GoTerms containing relevant p-values
     * @param mode char indicating which p-values to be used
     */
    private static void benjaminiHochberg2(ArrayList<GoTerm> notDiscarded, char mode) {
        int n = notDiscarded.size();
        double p;
        int rank = 0;
        GoTerm prevTerm = notDiscarded.get(0);
        if(mode == 'h') {
            // hg
            notDiscarded.sort(Comparator.comparing(term -> term.getHgP()));
            double prevP = prevTerm.hgP;
            for(int i = 0; i < n; i++) {
                GoTerm currTerm = notDiscarded.get(i);
                if(prevP != currTerm.hgP) {
                    rank = i;
                }
                p = currTerm.hgP * ((double)n/(double)rank);
                if(Double.isNaN(p) || p < 0) {
                    p = 0.0d;
                } else if(p > 1) {
                    p = 1.0d;
                }
                currTerm.hgPCorr = p;
                prevP = currTerm.hgP;
            }
        } else if(mode == 'f') {
            // fej
            notDiscarded.sort(Comparator.comparing(term -> term.getFejP()));
            double prevP = prevTerm.fejP;
            for(int i = 0; i < n; i++) {
                GoTerm currTerm = notDiscarded.get(i);
                if(prevP != currTerm.fejP) {
                    rank = i;
                }
                p = currTerm.fejP * ((double)n/(double)rank);
                if(Double.isNaN(p) || p < 0) {
                    p = 0.0d;
                } else if(p > 1) {
                    p = 1.0d;
                }
                currTerm.fejPCorr = p;
                prevP = currTerm.fejP;
            }
        } else {
            // ks
            notDiscarded.sort(Comparator.comparing(term -> term.getkP()));
            double prevP = prevTerm.kP;
            for(int i = 0; i < n; i++) {
                GoTerm currTerm = notDiscarded.get(i);
                if(prevP != currTerm.kP) {
                    rank = i;
                }
                p = currTerm.kP * ((double)n/(double)rank);
                if(Double.isNaN(p) || p < 0) {
                    p = 0.0d;
                } else if(p > 1) {
                    p = 1.0d;
                }
                currTerm.kPCorr = p;
                prevP = currTerm.kP;
            }
        }
    }

    /**
     * Calculates and sets the needed statistics in all nodes, for later multiple testing correction.
     */
    private static int calculateStatistics() {
        collectAllFoldChanges();
        int discardedNodes = 0;

        /** example */
        /*GoTerm example = dag.get("GO:0002683");
        int[] eGeneSizes = example.getMeasuredSize();
        int eNMeasured = eGeneSizes[0];
        int eNOverlap = eGeneSizes[1];
        // ORA
        example.calcOraPValues(N, K, eNMeasured, eNOverlap);
        // comparing distributions
        // example.calcDistPValues(allGenes, N, eNMeasured);
        int x = 4;  // breakpoint*/
        /** example end */

        for(GoTerm node : dag.values()) {
            if(node.geneSet.isEmpty()) {
                discardedNodes++;
                node.discardThis = true;
                continue;    // skip nodes without genes
            }
            int size = node.geneSet.size();
            if(size >= minSize && size <= maxSize) {
                int[] geneSizes = node.getMeasuredSize();
                int nMeasured = geneSizes[0];
                int nOverlap = geneSizes[1];
                // ORA
                node.calcOraPValues(N, K, nMeasured, nOverlap);
                // comparing distributions
                node.calcDistPValues(allGenes, N, nMeasured);
            } else {
                discardedNodes++;
                node.discardThis = true;
            }
        }
        return discardedNodes;
    }

    /**
     * Makes every node collect its fold changes.
     */
    private static void collectAllFoldChanges() {
        for(GoTerm node : dag.values()) {
            node.collectFoldChanges();
        }
    }

    /**
     * Collects the enrichment information of genes to all genes.
     * @throws IOException
     */
    private static void addEnrichment() throws IOException {
        BufferedReader brEn = new BufferedReader(new FileReader(enrichFile));
        String line;
        String[] splitLine;
        String id = "";
        String termID;
        GoTerm trueTerm;
        Gene currentGene;
        while((line = brEn.readLine()) != null) {
            if(line.startsWith("#")) {
                // annotate truly enriched GO terms
                termID = line.strip().substring(1);
                trueTerm = dag.get(termID);
                trueTerm.simulTrue = true;
                trues.add(trueTerm);
            } else {
                // add enrichment info to all genes
                if(line.startsWith("id\t")) continue;   // skip header
                splitLine = line.split("\t");
                id = splitLine[0];
                if(id.equals("CCHCR1")) {
                    int x = 3;  // breakpoint
                }
                currentGene = allGenes.get(id);
                if(currentGene != null) {
                    currentGene.foldChange = Float.parseFloat(splitLine[1]);
                    currentGene.significant = Boolean.parseBoolean(splitLine[2]);
                }
            }
        }
        brEn.close();
    }

    /**
     * Adds the associated genes to the DAG nodes.
     * @throws IOException
     */
    private static void geneMapping() throws IOException {
        BufferedReader brMap;
        if(mapTypeGo) {
            GZIPInputStream inputStream = new GZIPInputStream(new FileInputStream(mapFileStr));
            brMap = new BufferedReader(new InputStreamReader(inputStream));
            readGo(brMap);
        } else {
            brMap = new BufferedReader(new FileReader(new File(mapFileStr)));
            readEnsembl(brMap);
        }
        brMap.close();
    }

    /**
     * Performs the mapping association routine for the Ensembl mapping type.
     * @param br
     * @throws IOException
     */
    private static void readEnsembl(BufferedReader br) throws IOException {
        String line = br.readLine();    // skip header
        String[] splitLine;
        String id = "";
        String[] associatedTerms;
        Gene currentGene;
        while((line = br.readLine()) != null) {
            splitLine = line.split("\t");
            id = splitLine[1];
            associatedTerms = splitLine[2].split("\\|");

            if(id.equals("")) continue;     // skip if gene id empty
            if(associatedTerms.length == 0) continue;   // skip if no go entries

            boolean isNewGene = false;
            if(!allGenes.containsKey(id)) {
                // create new gene
                currentGene = new Gene(id);
                isNewGene = true;
            } else {
                // get existing gene
                currentGene = allGenes.get(id);
            }
            boolean addedToAtLeastOne = false;
            for(String termID : associatedTerms) {
                // recursive adding to DAG nodes
                GoTerm term = dag.get(termID);
                if(term != null) {
                    // else: term was discarded (obsolete, other namespace, ...)
                    term.addGene(currentGene);
                    addedToAtLeastOne = true;
                }
            }
            if(isNewGene && addedToAtLeastOne) {
                allGenes.put(id, currentGene);
            }
        }
    }

    /**
     * Performs the mapping association routine for the GO mapping type.
     * @param br
     * @throws IOException
     */
    private static void readGo(BufferedReader br) throws IOException {
        String line;
        String[] splitLine;
        String geneID = "";
        String termID = "";
        Gene currentGene;
        while((line = br.readLine()) != null) {
            if(line.startsWith("!")) continue;      // skip comment lines
            splitLine = line.split("\t");
            if(!splitLine[3].equals("")) continue;   // no association qualifier modifier allowed
            termID = splitLine[4];
            geneID = splitLine[2];
            if(termID.equals("") || geneID.equals("")) continue;    // skip if ids are empty

            boolean isNewGene = false;
            if(!allGenes.containsKey(geneID)) {
                // create new gene
                currentGene = new Gene(geneID);
                isNewGene = true;
            } else {
                // get existing gene
                currentGene = allGenes.get(geneID);
            }
            // recursive adding to DAG nodes
            GoTerm term = dag.get(termID);
            if(term != null) {
                term.addGene(currentGene);
                if(isNewGene) {
                    allGenes.put(geneID, currentGene);
                }
            }
        }
    }

    /**
     * Builds the DAG from the collected entries in the dag-HashMap.
     */
    private static void buildDag() {
        for(String id : dag.keySet()) {
            GoTerm node = dag.get(id);
            // for each node in the DAG

            String[] parentIDs = node.parentIDs;
            GoTerm[] addedParents = new GoTerm[parentIDs.length];
            for(int i = 0; i < parentIDs.length; i++) {
                // add pointers to the node's parents-array
                String pID = parentIDs[i];
                addedParents[i] = dag.get(pID);
            }
            node.parents = addedParents;
            if(node.parents.length == 0) {
                // must be root
                rootNode = node;
            }
        }
    }

    /**
     * Collects GO entries from the given .obo-file regarding the given namespace (root).
     * @throws IOException
     */
    private static void collectTerms() throws IOException {
        BufferedReader brObo = new BufferedReader(new FileReader(oboFile));
        String line;
        String id = "";
        String name = "";
        ArrayList<String> parentsList = new ArrayList();
        boolean discard = true;
        while((line = brObo.readLine()) != null) {
            if(line.strip().equals("[Term]")) {
                // finished GoTerm: store
                if(!discard) {
                    String[] parentsArray = new String[0];
                    if(!parentsList.isEmpty()) {
                        parentsArray = parentsList.toArray(parentsArray);
                    } else {
                        parentsArray = new String[0];
                    }
                    GoTerm term = new GoTerm(id, name, parentsArray);
                    dag.put(id, term);
                }
                // new GoTerm: load
                discard = false;
                id = "";
                name = "";
                parentsList = new ArrayList<>();
            } else if(line.startsWith("id: ")) {
                id = line.strip().substring(4);
            } else if(line.startsWith("namespace: ")) {
                if(line.charAt(11) != root) {
                    discard = true;
                }
            } else if(line.strip().equals("is_obsolete: true")) {
                discard = true;
            } else if(line.startsWith("name: ")) {
                name = line.strip().substring(6);
            } else if(line.startsWith("is_a: ")) {
                parentsList.add(line.strip().substring(6, 16));
            } else if(line.strip().equals("[Typedef]")) break;
        }
        // finished GoTerm: store
        if(!discard) {
            String[] parentsArray = new String[0];
            if(!parentsList.isEmpty()) {
                parentsArray = parentsList.toArray(parentsArray);
            } else {
                parentsArray = new String[0];
            }
            GoTerm term = new GoTerm(id, name, parentsArray);
            dag.put(id, term);
        }
        brObo.close();
    }

    /**
     * Prints usage when program gets invoked with false parameters.
     *
     */
    private static void printUsage() {
        // parameters
        System.out.println("Following parameters (in any order) are necessary:\n" +
                "-obo (file name): the DAG structure\n" +
                "-mapping (file name): defines gene <-> GO-class association\n" +
                "-enrich (file name): gene info for enrichment analysis" +
                "-o (file name): output tsv\n" +
                "-root (string): GO namespace specifying which DAG to use" +
                "-mappingtype (string): format how mapping is provided, either ensembl or go" +
                "-minsize (int): only consider GO entries with at least minsize associated genes" +
                "-maxsize (int): only consider GO entries with at most maxsize associated genes" +
                "*OPTIONAL* -overlapout (file name): additional output tsv"
        );
        // what will be done
        System.out.println("If invoked correctly, the program will annotate enrichment analysis on simulated data. " +
                "Optionally, the program will analyze GO-DAG overlap properties. " +
                "The output contains one line per GO entry with 13 columns. " +
                "If overlapout is not null, this output file will have one line per overlapping pair of DAG entries.");
    }

    public static void main(String[] args) {
        //java.awt.Toolkit.getDefaultToolkit().beep();
        // argparse
        String outputTsvString = null;
        String overlapOutTsvString = "";
        try {
            if(args.length < 6) {
                printUsage();
            } else if(args.length == 16 || args.length == 18){
                for(int i = 0; i < args.length - 1; i += 2) {
                    if(args[i].equals("-obo")) {
                        oboFile = new File(args[i + 1]);
                    } else if(args[i].equals("-mapping")) {
                        mapFileStr = args[i + 1];
                    } else if(args[i].equals("-enrich")) {
                        enrichFile = new File(args[i + 1]);
                    } else if(args[i].equals("-o")) {
                        outputTsvString = args[i + 1];
                    } else if(args[i].equals("-overlapout")) {
                        overlapOutTsvString = args[i + 1];
                    } else if(args[i].equals("-root")) {
                        String rootStr = args[i + 1];
                        if(rootStr.startsWith("bio")) {
                            // namespace = biological process
                            root = 'b';
                        } else if(rootStr.startsWith("mol")) {
                            // namespace = molecular function
                            root = 'm';
                        } else if(rootStr.startsWith("cel")) {
                            // namespace = cellular component
                            root = 'c';
                        }
                    } else if(args[i].equals("-mappingtype")) {
                        if(args[i + 1].equals("go")) {
                            mapTypeGo = true;
                        } else {
                            mapTypeGo = false;
                        }
                    } else if(args[i].equals("-minsize")) {
                        minSize = Integer.parseInt(args[i + 1]);
                    } else if(args[i].equals("-maxsize")) {
                        maxSize = Integer.parseInt(args[i + 1]);
                    } else {
                        printUsage();
                    }
                }
                performEnrichment(outputTsvString, overlapOutTsvString);
            } else {
                printUsage();
            }
        } catch (IOException e) {
            System.err.println("Error: Something went wrong.");
        }
    }
}
