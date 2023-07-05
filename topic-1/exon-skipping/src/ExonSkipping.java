import java.io.*;
import java.util.*;


public class ExonSkipping {
    private static HashMap<String, Gene> genes = new HashMap<>();
    private static String outputFilePath = "";

    /**
     * Reads and analyzes gtf files.
     *
     * @param fileName
     * @return
     * @throws IOException
     */
    public static void readFile(String fileName) throws IOException {
        try {
            if(fileName.endsWith(".gtf") == false) {
                throw new Exception("Please only use GTF-files!");
            }
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = "";
        String attributesLine = "";

        // help attributes for analyzing reasons
        /*HashMap<String, Integer[]> startAndEndOfExons = new HashMap<>();
        HashMap<String, Region> allCDSs = new HashMap<>();
        HashMap<String, Region> allExons = new HashMap<>();
        HashMap<String, String> allTranscriptsToGene = new HashMap<>();
        */

        while((line = br.readLine()) != null) {
            // get all fields
            String[] blocks = line.split("\t");

            // check if it's a valid entry
            if (blocks.length < 9) {
                continue;
            }

            attributesLine = blocks[blocks.length - 1];
            String[] attributes = cutEmptySpace(attributesLine.split(";"));

            // help methods here!
            /*analyzeCDSandExons(blocks, attributes, startAndEndOfExons);
            analyzeGenesAndStrands(blocks, attributes);
            analyzeExonFakeIDs(blocks, attributes, allExons, allCDSs);
            analyzeCDSbeforeTranscript(blocks, attributes, allTranscriptsToGene);*/

            // collect cdss for genes/transcripts
            if(blocks[2].equals("CDS") || blocks[2].equals("exon")) {
                String transcriptID = findArgument("transcript_id", attributes);
                String geneID = findArgument("gene_id", attributes);
                String idFake = geneID + "_" + transcriptID + "_" + findArgument("exon_number", attributes);

                Gene currentGene;
                // check if gene exists
                if(!genes.containsKey(geneID)) {
                    // create new gene and add
                    currentGene = new Gene();
                    currentGene.id = geneID;
                    currentGene.chromosome = blocks[0];
                    currentGene.strand = blocks[6].charAt(0);
                    currentGene.name = findArgument("gene_name", attributes);
                    genes.put(geneID, currentGene);
                } else {
                    // get existing gene
                    currentGene = genes.get(geneID);
                }

                Transcript currentTranscript;
                // check if transcript exists
                if (!currentGene.transcriptsAll.containsKey(transcriptID)) {
                    // create new transcript and add
                    currentTranscript = new Transcript(transcriptID);
                    currentTranscript.proteinID = findArgument("protein_id", attributes);
                    currentGene.transcriptsAll.put(transcriptID, currentTranscript);
                } else {
                    // get existing transcript
                    currentTranscript = currentGene.transcriptsAll.get(transcriptID);
                }

                int x1 = Integer.parseInt(blocks[3]);
                int x2 = Integer.parseInt(blocks[4]);
                Region tempRegion = new Region(x1, x2, transcriptID); // dont forget belongsTo

                if(blocks[2].equals("CDS")) {
                    // check if CDS is already in list of its transcript
                    if(!currentTranscript.cdss.regions.contains(tempRegion)) {
                        currentTranscript.cdss.regions.add(tempRegion);
                        currentTranscript.cdss.ids.add(idFake);
                        if(!currentTranscript.containsCDS) {
                            currentGene.numberOfCDStanscr++;
                            currentTranscript.containsCDS = true;
                        }
                    }
                    if(currentTranscript.proteinID.equals("")) {
                        // put protein ID if missing
                        currentTranscript.proteinID = findArgument("protein_id", attributes);
                        if(currentTranscript.proteinID.equals("ENSP00000381283")) {
                            System.out.println(transcriptID);
                        }
                    }
                } else {
                    // check if exon is already in list of its transcript
                    if(!currentTranscript.exons.regions.contains(tempRegion)) {
                        currentTranscript.exons.regions.add(tempRegion);
                        currentTranscript.exons.ids.add(idFake);
                    }
                }
            }
        }
    }

    /**
     * Calculates all introns, once the file is read and the data is collected.
     */
    private static void calculateAllIntrons() {
        int sumIntrons = 0;
        for(String currentGeneID : genes.keySet()) {
            Gene currentGene = genes.get(currentGeneID);
            currentGene.calcIntrons();
            if(currentGeneID.equals("ENSG00000131018")) {
                // set breakpoint here to check
                for(Transcript t : currentGene.transcriptsAll.values()) {
                    if(t.containsCDS) {
                        sumIntrons += t.introns.size();
                    }
                }
                // System.out.println(sumIntrons);      // debugging reasons
            }
            currentGene.collectIntrons();
            currentGene = currentGene;
        }
    }

    /**
     * Extracts all exon skipping events of all Genes.
     */
    public static void extractExonSkippings() throws IOException {
        boolean hasStart = false;
        boolean hasEnd = false;
        Integer transcrIntronEnd = -1;
        boolean foundExSkipEvent = false;

        File outputFile = new File(outputFilePath);
        outputFile.createNewFile();
        FileOutputStream fileOutStream = new FileOutputStream(outputFile);
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fileOutStream));
        bw.write("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon\t" +
                "max_skipped_exon\tmin_skipped_bases\tmax_skipped_bases");
        bw.newLine();

        perGene: for(String currentGeneID : genes.keySet()) {
            // doing this for each gene
            Gene currentGene = genes.get(currentGeneID);
            TreeSet<Region> intronsAllSorted = sortCoordinatesWithOverlaps(currentGene.intronsAll.regions);

            if(currentGene.intronsAll.regions.size() > 0 && currentGene.transcriptsAll.size() > 0) {
                // only makes sense if gene has introns

                if(currentGeneID.equals("ENSG00000160191")) {
                    currentGene = currentGene;  // put breakpoint here
                }

                perIntronPerGene: for(Region currentIntron : intronsAllSorted) {
                    // doing this for each intron of the gene

                    if(currentIntron.x1 == 123994512 && currentIntron.x2 == 124012305) {
                        currentIntron = currentIntron;    // put breakpoint here
                    }

                    if(currentGeneID.equals("ENSG00000160191") &&
                            currentIntron.belongsTo.equals("ENST00000398227")) {
                        currentIntron = currentIntron; // put breakpoint (A1) here
                    }

                    HashSet<Region> WT_Introns = new HashSet<>();
                    HashSet<Region> WT_Exons = new HashSet<>();
                    TreeMap<Integer, Integer> sortedExonsOfTranscr = new TreeMap<>();
                    TreeSet<String> WT_proteinIDs = new TreeSet<>();
                    TreeSet<String> SV_proteinIDs = new TreeSet<>();
                    Integer minSkippedExon = Integer.MAX_VALUE;
                    int maxSkippedExon = 0;
                    Integer minSkippedBases = Integer.MAX_VALUE;
                    int maxSkippedBases = 0;

                    wildtypeSearch: for(String currentTranscriptID : currentGene.transcriptsAll.keySet()) {
                        // find out if this is a wildtype for a splicevariant
                        Transcript currentTranscript = currentGene.transcriptsAll.get(currentTranscriptID);

                        if(currentTranscript.proteinID.equals("ENSP00000381291")) {
                            currentTranscript = currentTranscript;  // put breakpoint (A2) here
                        }

                        if(currentTranscript.cdss.regions.size() >= 2) {
                            // only makes sense to look if transcript contains at least 2 CDSs (1 intron)

                            if(currentTranscript.introns.get(currentIntron.x1) != null &&
                                    currentTranscript.introns.get(currentIntron.x1) == currentIntron.x2) {
                                // this transcript contains the exact currentIntron -> is a splicevariant
                                SV_proteinIDs.add(currentTranscript.proteinID);
                            } else {
                                // this transcript is a wildtype candidate
                                inner: for(Integer transcrIntronStart : currentTranscript.introns.keySet()) {
                                    transcrIntronEnd = currentTranscript.introns.get(transcrIntronStart);
                                    if(!hasStart) {
                                        // not reached the start yet
                                        if (currentIntron.x1 < transcrIntronStart) break inner;  // cant find a start
                                        if (currentIntron.x1 == transcrIntronStart &&
                                                currentIntron.x2 > transcrIntronEnd) {
                                            // found an intron that has the same start but not the same end
                                            hasStart = true;
                                        }
                                    } else {
                                        // looking for the end now
                                        if(currentIntron.x2 < transcrIntronEnd) break inner;    // cant find an end
                                        if(currentIntron.x2 == transcrIntronEnd
                                                && currentIntron.x1 < transcrIntronStart) {
                                            // found an intron that has the same end but not the same start
                                            hasEnd = true;
                                            break inner;
                                        }
                                    }
                                }
                            }
                            if(hasStart && hasEnd) {
                                // found an exonskipping event! currentTranscript is a wildtype for currentIntron
                                WT_proteinIDs.add(currentTranscript.proteinID);

                                // look for all skipped introns
                                for (Integer skippedIntronStart : currentTranscript.introns.keySet()) {
                                    Integer skippedIntronEnd = currentTranscript.introns.get(skippedIntronStart);
                                    if (skippedIntronStart >= currentIntron.x1 &&
                                            skippedIntronEnd <= currentIntron.x2) {
                                        // found a skipped intron from wildtype
                                        Region WT_intron = new Region(skippedIntronStart, skippedIntronEnd,
                                                currentTranscriptID);
                                        WT_Introns.add(WT_intron);
                                    }
                                }

                                // look for all skipped exons
                                sortedExonsOfTranscr = currentTranscript.exons.getCoordinatesSorted();
                                int numSkippedExon = 0; // for this concrete WT
                                int numSkippedBases = 0; // for this concrete WT
                                for(Integer skippedExonStart : sortedExonsOfTranscr.keySet()) {
                                    Integer skippedExonEnd = sortedExonsOfTranscr.get(skippedExonStart);
                                    if(skippedExonStart > currentIntron.x1 && skippedExonEnd < currentIntron.x2) {
                                        // found a skipped exon from wildtype
                                        Region WT_exon = new Region(skippedExonStart, skippedExonEnd,
                                                currentTranscriptID);
                                        WT_Exons.add(WT_exon);
                                        numSkippedExon++;
                                        numSkippedBases += skippedExonEnd - skippedExonStart + 1;
                                    }
                                }
                                if(maxSkippedExon < numSkippedExon) {
                                    maxSkippedExon = numSkippedExon;
                                }
                                if(minSkippedExon > numSkippedExon) {
                                    minSkippedExon = numSkippedExon;
                                }
                                if(maxSkippedBases < numSkippedBases) {
                                    maxSkippedBases = numSkippedBases;
                                }
                                if(minSkippedBases > numSkippedBases) {
                                    minSkippedBases = numSkippedBases;
                                }

                                foundExSkipEvent = true;
                            }
                            hasEnd = hasStart = false;
                        }
                    }
                    if(foundExSkipEvent) {
                        // start creating output entry here

                        /*
                        if(currentGeneID.equals("ENSG00000131018")) {
                            String outputLine = createOutputEntry(currentGene, currentIntron.x1, currentIntron.x2,
                                    WT_Introns, WT_proteinIDs, SV_proteinIDs,
                                    minSkippedExon, maxSkippedExon, minSkippedBases, maxSkippedBases);
                            System.out.println(outputLine);
                        }
                        */
                        bw.write(createOutputEntry(currentGene, currentIntron.x1, currentIntron.x2,
                                WT_Introns, WT_proteinIDs, SV_proteinIDs,
                                minSkippedExon, maxSkippedExon, minSkippedBases, maxSkippedBases));
                        bw.newLine();
                    }
                    foundExSkipEvent = false;
                }
            }
        }
        bw.close();
    }


    private static String createOutputEntry(Gene gene, int intronSVstart, int intronSVend,
                                            HashSet<Region> WT_introns, TreeSet<String> WT_proteinIDs,
                                            TreeSet<String> SV_proteinIDs,
                                            int minSkippedExons, int maxSkippedExons,
                                            int minSkippedBases, int maxSkippedBases) {
        String line = "";
        // basic info
        line += gene.id + "\t" + gene.name + "\t" + gene.chromosome + "\t" + gene.strand;
        // number of annotated CDS and transcripts in the gene
        line += "\t" + gene.numberOfCDStanscr;
        line += "\t" + gene.transcriptsAll.size();
        // the SV intron as start:end
        line += "\t" + intronSVstart + ":" + intronSVend + "\t";
        // the WT introns within the SV intron as start:end, separated by |
        TreeSet<Region> WT_introns_sorted = sortCoordinatesWithOverlaps(WT_introns);
        for(Region WT_intron : WT_introns_sorted) {
            line += WT_intron.x1 + ":" + WT_intron.x2 + "|";
        }
        if(line.endsWith("|")) {
            // cut last |
            line = line.substring(0, line.length() - 1);
        }
        // ids of the WT CDSs, separated by |
        line += "\t";
        for(String WT_proteinID : WT_proteinIDs) {
            line += WT_proteinID + "|";
        }
        if(line.endsWith("|")) {
            // cut last |
            line = line.substring(0, line.length() - 1);
        }
        // ids of the SV CDSs, separated by |  (protein_id)
        line += "\t";
        for(String SV_proteinID : SV_proteinIDs) {
            line += SV_proteinID + "|";
        }
        if(line.endsWith("|")) {
            // cut last |
            line = line.substring(0, line.length() - 1);
        }
        // the minimal number of skipped exons in any WT/SV pair
        line += "\t" + minSkippedExons;
        // the maximum number of skipped exons
        line += "\t" + maxSkippedExons;
        // the minimal number of skipped bases (joint lengths of skipped exons) in any WT/SV pair
        line += "\t" + minSkippedBases;
        // the maximum number of skipped bases
        line += "\t" + maxSkippedBases;
        return line;
    }

    /**
     * Compares sample output to my output.
     * @throws IOException
     */
    private static void compareToSampleOutput() throws IOException {
        // collect info from sample file
        ArrayList<String> WT = new ArrayList<>();
        ArrayList<String[]> WT_prots = new ArrayList<>();
        ArrayList<String[]> SV_prots = new ArrayList<>();
        BufferedReader sampleBR = new BufferedReader(new FileReader("sample-output.txt"));
        String sampleLine = sampleBR.readLine();
        while((sampleLine = sampleBR.readLine()) != null) {
            String[] sampleFields = sampleLine.split("\t");
            WT.add(sampleFields[7]);
            WT_prots.add(sampleFields[8].split("\\|"));
            SV_prots.add(sampleFields[9].split("\\|"));
        }

        // collect info from my file
        ArrayList<String> rWT = new ArrayList<>();
        ArrayList<String[]> rWT_prots = new ArrayList<>();
        ArrayList<String[]> rSV_prots = new ArrayList<>();
        BufferedReader rBR = new BufferedReader(new FileReader("test2.txt"));
        String rLine = rBR.readLine();
        while((rLine = rBR.readLine()) != null) {
            String[] rFields = rLine.split("\t");
            rWT.add(rFields[7]);
            rWT_prots.add(rFields[8].split("\\|"));
            rSV_prots.add(rFields[9].split("\\|"));
        }

        // compare to info from my output

        // same WT coordinates?
        String wrongWTs = "";
        for(int i = 0; i < rWT.size(); i++) {
            if(!rWT.get(i).equals(WT.get(i))) {
                wrongWTs += i + ", ";
            }
        }
        System.out.println(wrongWTs);
        // see if WT_prots of raw are missing in sample
        for(int i = 0; i < rWT_prots.size(); i++) {
            for(int j = 0; j < rWT_prots.get(i).length; j++) {
                String rWT_prot = rWT_prots.get(i)[j];
                boolean foundProt = false;
                for(int k = 0; k < WT_prots.get(i).length; k++) {
                    if(rWT_prot.equals(WT_prots.get(i)[k])) {
                        foundProt = true;
                    }
                }
                if(!foundProt) {
                    System.err.println("Couldn't find " + rWT_prot + " in sample!");
                }
            }
        }
        // see if WT_prots of sample are missing in my output
        for(int i = 0; i < WT_prots.size(); i++) {
            for(int j = 0; j < WT_prots.get(i).length; j++) {
                String WT_prot = WT_prots.get(i)[j];
                boolean foundProt = false;
                for(int k = 0; k < rWT_prots.get(i).length; k++) {
                    if(WT_prot.equals(rWT_prots.get(i)[k])) {
                        foundProt = true;
                    }
                }
                if(!foundProt) {
                    System.err.println("Couldn't find " + WT_prot + " in my output!");
                }
            }
        }
    }

    /**
     *
     * @param search the concerned attribute
     * @param attributes last field
     * @return the attribute value or the empty string, if nothing could be found
     */
    private static String findArgument(String search, String[] attributes){
        // finds argument value in an array of * argument_name "value" * blocks
        for(int i = 0; i < attributes.length; i++) {
            if(attributes[i].startsWith(search)) {
                if(search.equals("exon_number")) {
                    return attributes[i].substring(attributes[i].length());
                } else {
                    return attributes[i].split("\"")[1];
                }
            }
        }
        return "";
    }

    /**
     * Converts an unsorted HashSet of Regions into a sorted TreeMap of Integers.
     * @param regions
     * @return TreeMap
     */
    private static TreeMap<Integer, Integer> sortCoordinates(HashSet<Region> regions) {
        TreeMap<Integer, Integer> sortedRegions = new TreeMap<>();
        for(Region reg : regions) {
            sortedRegions.put(reg.x1, reg.x2);
        }
        return sortedRegions;
    }


    /**
     * Converts an unsorted HashSet of Regions into a sorted ArrayList of Integer[]s.
     * @param regions
     * @return TreeSet
     */
    private static TreeSet<Region> sortCoordinatesWithOverlaps(HashSet<Region> regions) {
        TreeSet<Region> sortedRegions = new TreeSet<>();
        for(Region reg : regions) {
            sortedRegions.add(reg);
        }
        return sortedRegions;
    }

    /**
     * Checks if exons/CDSs of the same regions have different start/end positions (yes).
     * @param blocks all fields
     * @param attributes last field
     * @param startAndEndOfExons collection
     */
    /*
    private static void analyzeCDSandExons(String[] blocks, String[] attributes,
                                           HashMap<String, Integer[]> startAndEndOfExons) {
        int exonStart = -1;
        int exonStop = -2;
        String geneID = "";
        String transcriptID = "";
        String exonNumberStr = "";
        if (blocks[2].equals("exon") || blocks[2].equals("CDS")) {
            exonStart = Integer.parseInt(blocks[3]);
            exonStop = Integer.parseInt(blocks[4]);

            geneID = attributes[0].split("\"")[1];
            transcriptID = attributes[1].split("\"")[1];
            exonNumberStr = attributes[2].split("\"")[1];

            // extra for analysis reasons: create exon key
            String newID = geneID + "_" + transcriptID + "_" + exonNumberStr;
            if(startAndEndOfExons.containsKey(newID)) {
                if (startAndEndOfExons.get(newID)[0] != exonStart || startAndEndOfExons.get(newID)[1] != exonStop) {
                    System.out.println(newID + " " + blocks[2] + ": " +
                            startAndEndOfExons.get(newID)[0] + " vs " + exonStart + ", " +
                            startAndEndOfExons.get(newID)[1] + " vs " + exonStop);
                }
            } else {
                Integer[] startAndEnd = new Integer[2];
                startAndEnd[0] = exonStart;
                startAndEnd[1] = exonStop;
                startAndEndOfExons.put(newID, startAndEnd);
            }
        }
    }
     */

    /**
     * Checks if the fake ID of gene_id, trancript_id and exon_number I'm giving all exons/CDSs is unique (yes).
     * @param blocks all fields
     * @param attributes last field
     * @param exons set of all exons
     * @param cdss set of all introns
     */
    /*
    private static void analyzeExonFakeIDs(String[] blocks, String[] attributes,
                                           HashMap<String, Region> exons, HashMap<String, Region> cdss) {
        // collect id
        String geneID = findArgument("gene_id", attributes);
        String transcriptID = findArgument("transcript_id", attributes);
        String exonNumber = findArgument("exon_number", attributes);
        String exonIDfake = geneID + "_" + transcriptID + "_" + exonNumber;
        if(blocks[2].equals("exon")) {
            // check if already saved
            if(!exons.containsKey(exonIDfake)) {
                exons.put(exonIDfake, new Region(Integer.parseInt(blocks[3]),
                        Integer.parseInt(blocks[4]), transcriptID));
            } else {
                Region exon = exons.get(exonIDfake);
                // check if regions are different
                int newx1 = Integer.parseInt(blocks[3]);
                int newx2 = Integer.parseInt(blocks[4]);
                if(exon.x1 != newx1 || exon.x2 != newx2) {
                    System.err.println("Same exon fake ID ( " + exonIDfake + " ), different region: " +
                            "o=" + exon.x1 + ", n=" + newx1 + " | o=" + exon.x2 + ", n=" + newx2);
                }
            }
        } else if(blocks[2].equals("CDS")) {
            // check if already saved
            if(!cdss.containsKey(exonIDfake)) {
                cdss.put(exonIDfake, new Region(Integer.parseInt(blocks[3]),
                        Integer.parseInt(blocks[4]), transcriptID));
            } else {
                Region cds = cdss.get(exonIDfake);
                // check if regions are different
                int newx1 = Integer.parseInt(blocks[3]);
                int newx2 = Integer.parseInt(blocks[4]);
                if(cds.x1 != newx1 || cds.x2 != newx2) {
                    System.err.println("Same exon fake ID ( " + exonIDfake + " ), different region: " +
                            "o=" + cds.x1 + ", n=" + newx1 + " | o=" + cds.x2 + ", n=" + newx2);
                }
            }
        }
    }
     */

    /**
     * Checks if there is an entry where CDS comes before its transcript.
     * @param blocks all fields
     * @param attributes last field
     */
    /*
    private static void analyzeCDSbeforeTranscript(String[] blocks, String[] attributes,
                                                   HashMap<String, String> allTranscriptsToGene) {
        String transcriptID = findArgument("transcript_id", attributes);
        if(blocks[2].equals("transcript")) {
            String geneID = findArgument("gene_id", attributes);
            if(allTranscriptsToGene.containsKey(transcriptID)) {
                if(allTranscriptsToGene.get(transcriptID) != geneID) {
                    System.err.println("Same transcript_id, different gene!");
                }
            } else {
                allTranscriptsToGene.put(transcriptID, geneID);
            }
        } else if(blocks[2].equals("CDS")) {
            if(!allTranscriptsToGene.containsKey(transcriptID)) {
                System.err.println("CDS before transcript!!!");
            }
        }
    }
     */

    /**
     * Analyzes if there are transcript entries that appear before their gene entries
     * (if you only enter protein_coding genes, than yes, that may happen!).
     * @param blocks all fields
     * @param attributes last field
     */
    /*
    private static void analyzeTranscriptBeforeGene(String[] blocks, String[] attributes) {
        String transcriptID = findArgument("transcript_id", attributes);
        String geneID = findArgument("gene_id", attributes);
        if(!genes.containsKey(geneID)) {
            System.err.println("Transcript " + transcriptID + " appears before its gene!");
        }
    }
    */

    /**
     * Cuts annoying empty spaces at the start of the semicolon-separated attributes.
     * @param words attributes/last field
     * @return the same field but without empty spaces at the beginning of certain attributes fields
     */
    private static String[] cutEmptySpace(String[] words) {
        for(int i = 0; i < words.length; i++) {
            if(words[i].startsWith(" ")) {
                words[i] = words[i].substring(1);
            }
        }
        return words;
    }

    /**
     * Prints usage if invoked without/with missing/with wrong parameters.
     */
    private static void printUsage() {
        System.out.println("Following parameters (in any order) are both necessary: \n" +
                "-o <output file path>: path to the output where the table (defined below)" +
                "containing all exon skippings will be written \n" +
                "-gtf <GTF file>: genomic annotation necessary for the task");
        System.out.println("If the parameters are correct, the program will calculate all ExonSkipping events " +
                "the given GTF file contains and print them in a tabular manner.");

    }

    public static void main (String[]args){
        try {
            /*
            readFile("Homo_sapiens.GRCh37.75.gtf");
            calculateAllIntrons();
            extractExonSkippings();
            compareToSampleOutput();
            Gene searchedGene = genes.get("ENSG00000131018");
            searchedGene = searchedGene;
            */

            // Argparse
            if(args.length == 0) {
                printUsage();
            } else if(args.length == 4){
                if(args[0].equals("-gtf")) {
                    if(args[1].endsWith(".gtf")) {
                        readFile(args[1]);
                        outputFilePath = args[3];
                        calculateAllIntrons();
                        extractExonSkippings();
                    } else {
                        printUsage();
                    }
                } else if(args[0].equals("-o")) {
                    // give args[1] to output method
                    if(args[3].endsWith(".gtf")) {
                        readFile(args[3]);
                        outputFilePath = args[1];
                        calculateAllIntrons();
                        extractExonSkippings();
                    } else {
                        printUsage();
                    }
                } else {
                    printUsage();
                }

            } else {
                printUsage();
            }
        } catch (IOException e) {
            System.err.println("Error: Something went wrong.");
        }
        /*
        Tests:

        SortedSet<Region> testset = new TreeSet<>();
        testset.add(new Region(1, 2, "blub"));
        testset.add(new Region(1, 2, "piep"));
        testset.add(new Region(1, 2, "blub"));

        TreeMap<Integer, Integer> treeMap = new TreeMap<>();
        treeMap.put(3, 30);
        treeMap.put(1, 10);
        treeMap.put(2, 20);
        for(Integer first : treeMap.keySet()) {
            System.out.println(first);
        }
        */
    }
}

