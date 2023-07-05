import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;


public class ReadSimulator {
    private static File gtf;
    private static File fasta;
    private static RandomAccessFile raf;
    private static File fidx;   // fasta file index
    private static File readcons;
    private static HashMap<String, String[]> fidxMap;
    private static HashMap<String, Integer> rankInFidx;
    private static HashMap<String, Long> based;      // possible offset
    private static ArrayList<Integer> current_tMutPos;
    private static HashMap<String, String> transcripts;
    private static HashMap<String, Boolean> myStrands;

    /**
     * Simulates reads.
     * @throws IOException
     */
    private static void simulateReads(float mean, float sd, String outputDir, int readLength, float mutRate)
            throws IOException {
        BufferedReader brReadcons = new BufferedReader(new FileReader(readcons));
        fidxMap = readFastaIndex();
        myStrands = new HashMap<>();
        transcripts = new HashMap<>();
        File fwFastq = new File(outputDir + "/fw.fastq");
        fwFastq.delete();   // delete if already existing
        fwFastq = new File(outputDir + "/fw.fastq");
        File rwFastq = new File(outputDir + "/rw.fastq");
        rwFastq.delete();   // delete if already existing
        rwFastq = new File(outputDir + "/rw.fastq");
        FileOutputStream fosFwFastq = new FileOutputStream(fwFastq);
        FileOutputStream fosRwFastq = new FileOutputStream(rwFastq);
        BufferedWriter bwFwFastq = new BufferedWriter(new OutputStreamWriter(fosFwFastq));
        BufferedWriter bwRwFastq = new BufferedWriter(new OutputStreamWriter(fosRwFastq));
        File mappingInfo = new File(outputDir + "/read.mappinginfo");
        mappingInfo.delete();
        mappingInfo = new File(outputDir + "/read.mappinginfo");
        FileOutputStream fosMapInfo = new FileOutputStream(mappingInfo);
        BufferedWriter bwMapInfo = new BufferedWriter(new OutputStreamWriter(fosMapInfo));
        bwMapInfo.write("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec" +
                "\tfw_mut\trw_mut");
        bwMapInfo.newLine();

        String line = brReadcons.readLine();
        Transcript t = null;
        String[] splitLine;
        String transcriptID;
        String geneID;
        int tStart = -1;
        int tEnd = -1;
        String chromosome = "blub";
        int readID = -1;

        while((line = brReadcons.readLine()) != null) {
            BufferedReader brGtf = new BufferedReader(new FileReader(gtf));

            // 1) get transcript sequence as string
            // ------------------------------------

            splitLine = line.split("\t");
            transcriptID = splitLine[1];
            geneID = splitLine[0];

            // 1.1) find coordinates made up of exons in gtf file

            String lineGtf;
            String[] splitLineGtf;
            String[] attributes;
            boolean foundGene = false;

            while((lineGtf = brGtf.readLine()) != null) {
                if(lineGtf.startsWith("#")) continue;
                splitLineGtf = lineGtf.split("\t");
                if(splitLineGtf[2].equals("transcript") && foundGene) {
                    attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                    if(findArgument("transcript_id", attributes).equals(transcriptID)) {
                        // found transcript
                        tStart = Integer.parseInt(splitLineGtf[3]);
                        tEnd = Integer.parseInt(splitLineGtf[4]);
                        t = new Transcript(transcriptID);
                        if(splitLineGtf[6].charAt(0) == '+') {
                            myStrands.put(transcriptID, true);
                        } else {
                            myStrands.put(transcriptID, false);
                        }
                        chromosome = splitLineGtf[0];
                    }
                } else if(splitLineGtf[2].equals("gene")) {
                    attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                    if(findArgument("gene_id", attributes).equals(geneID)) {
                        // found gene of transcript, can look for transcript
                        foundGene = true;
                    } else {
                        if(foundGene) {
                            // already found relevant gene, can stop going through gtf
                            break;
                        }
                    }
                } else if(splitLineGtf[2].equals("exon") && t != null) {
                    attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                    if(findArgument("transcript_id", attributes).equals(transcriptID)) {
                        // found CDS of transcript
                        t.exons.regions.add(new Region(Integer.parseInt(splitLineGtf[3]),
                                Integer.parseInt(splitLineGtf[4]), transcriptID));
                    }
                }
            }

            // set start and end of transcript according to exons

            t.start = t.exons.regions.first().x1;
            // need to find end manually
            int biggestEnd = 0;
            for(Region r : t.exons.regions) {
                if(r.x2 > biggestEnd) {
                    biggestEnd = r.x2;
                }
            }
            t.end = biggestEnd;

            // 1.2) get sequence of exons of transcripts in fasta by using index

            String seq = getSequence(t, chromosome);

            // check if antisense strand needed

            if(!myStrands.get(transcriptID)) {
                seq = reverseComplement(seq);
            }

            /*System.out.println(seq);*/
            transcripts.put(transcriptID, seq);

            // 6) simulate mutations with required rate
            // ----------------------------------------

            String originalSeq = seq;
            seq = mutate(seq, mutRate);

            // 2) sample fragment length
            // -------------------------

            org.apache.commons.math3.distribution.NormalDistribution nd = new NormalDistribution(mean, (double)sd);

            // simulate given number of fragments (2 reads per fragment)

            for(int i = 0; i < Integer.parseInt(splitLine[2]); i++) {
                readID += 1;
                double sampleProp = nd.sample();
                int fl = Math.max(readLength, (int) sampleProp);
                // check if fragment is bigger than transcript
                if (fl > seq.length()) {
                    fl = seq.length();
                }

                // 3) select random position (from 0, length(t) - FL)
                // --------------------------------------------------

                int startPos = (int)(Math.random() * (seq.length() - fl));  // 0 based !

                // 4) get fragment sequence
                // ------------------------

                String fragSeq = seq.substring(startPos, startPos + fl); // 0 based, substring: end exclusive

                // 5) get read sequences of length readlength
                // ------------------------------------------

                String fwRead = fragSeq.substring(0, readLength);   // 0 based, substring: end exclusive
                String orgRead = originalSeq.substring(startPos, startPos + readLength);
                String revComp = reverseComplement(fragSeq);    // "umdrehen und umkehren"
                String rwRead = revComp.substring(0, readLength);   // 0 based, substring: end exclusive
                String orgRwRead = reverseComplement(originalSeq.substring(startPos, startPos + fl)).substring(0,
                        readLength);
                current_tMutPos = current_tMutPos; // breakpoint

                // 7-8) recreate positions (strands /!\)
                // -------------------------------------

                RegionVector fwReadGenome;
                Region fwReadTranscriptome;
                ArrayList<Integer> fwReadMutPos;
                RegionVector rwReadGenome;
                Region rwReadTranscriptome;
                ArrayList<Integer> rwReadMutPos;
                int rwReadStartPos;
                if(!myStrands.get(transcriptID)){
                    // - strand transcript
                    int firstReadEnd = startPos + readLength - 1;
                    int lastTranscrIndex = seq.length() - 1;
                    int trueStartPos = lastTranscrIndex - firstReadEnd;

                    fwReadGenome = recreateGenome(t, readLength, trueStartPos);

                    int trueRwStartPos = seq.length() - (startPos + fl);

                    rwReadGenome = recreateGenome(t, readLength, trueRwStartPos);
                    rwReadStartPos = (startPos + fl - 1) - readLength + 1; // for later
                } else {
                    // + strand transcript
                    fwReadGenome = recreateGenome(t, readLength, startPos);

                    rwReadStartPos = (startPos + fl - 1) - readLength + 1;
                    rwReadGenome = recreateGenome(t, readLength, rwReadStartPos);
                }
                // same transcriptome indices for - transcripts?
                fwReadTranscriptome = new Region(startPos, startPos + readLength - 1, transcriptID);
                fwReadMutPos = transPosToReadPos(current_tMutPos, startPos, readLength);
                rwReadTranscriptome = new Region(rwReadStartPos, (startPos + fl - 1), transcriptID);
                rwReadMutPos = revTransPosToReadPos(current_tMutPos, rwReadStartPos, readLength);

                // 9) write files
                // --------------

                writeFastq(bwFwFastq, readID, fwRead);
                writeFastq(bwRwFastq, readID, rwRead);
                writeMapInfo(bwMapInfo, readID, chromosome, geneID, transcriptID,
                        fwReadTranscriptome, rwReadTranscriptome,
                        fwReadGenome, rwReadGenome,
                        fwReadMutPos, rwReadMutPos);

                // debug: write reads

                /*
                System.out.println(readID);
                String muts = ",";
                for(Integer index : fwReadMutPos) {
                    muts += index + ",";
                }
                muts = muts.substring(0, muts.length() - 1);
                System.out.println("\tfw: " + fwRead + "\t" + muts);
                System.out.println("\t    " + orgRead);
                muts = ",";
                for(Integer index : rwReadMutPos) {
                    muts += index + ",";
                }
                muts = muts.substring(0, muts.length() - 1);
                System.out.println("\trw: " + rwRead + "\t" + muts);
                System.out.println("\t    " + orgRwRead);
                muts = muts;
                */
            }
        }
        bwFwFastq.close();
        bwRwFastq.close();
        bwMapInfo.close();
    }

    /**
     * Writes one line for the tabular output file read.mappinginfo.
     * @param bw
     * @param readID
     * @param chr
     * @param geneID
     * @param transcriptID
     * @param t_fw_regvec
     * @param t_rw_regvec
     * @param fw_regvec
     * @param rw_regvec
     * @throws IOException
     */
    private static void writeMapInfo(BufferedWriter bw, int readID, String chr, String geneID, String transcriptID,
                                     Region t_fw_regvec, Region t_rw_regvec,
                                     RegionVector fw_regvec, RegionVector rw_regvec,
                                     ArrayList<Integer> fw_mut, ArrayList<Integer> rw_mut
                                     ) throws IOException {
        String outputLine = readID + "\t" + chr + "\t" + geneID + "\t" + transcriptID + "\t";
        int endExclusive;
        endExclusive = t_fw_regvec.x2 + 1;
        outputLine += t_fw_regvec.x1 + "-" + endExclusive + "\t";   // 0 based
        endExclusive = t_rw_regvec.x2 + 1;
        outputLine += t_rw_regvec.x1 + "-" + endExclusive + "\t";   // 0 based
        for(Region reg : fw_regvec.regions) {                       // 1 based
            endExclusive = reg.x2 + 1;
           outputLine += reg.x1 + "-" + endExclusive + "|";
        }
        outputLine = outputLine.substring(0, outputLine.length() - 1);
        outputLine += "\t";
        for(Region reg : rw_regvec.regions) {                       // 1 based
            endExclusive = reg.x2 + 1;
            outputLine += reg.x1 + "-" + endExclusive + "|";
        }
        outputLine = outputLine.substring(0, outputLine.length() - 1);
        outputLine += "\t";
        for(Integer index : fw_mut) {
            outputLine += index + ",";
        }
        if(!fw_mut.isEmpty()) {
            outputLine = outputLine.substring(0, outputLine.length() - 1);
        }
        outputLine += "\t";
        for(Integer index : rw_mut) {
            outputLine += index + ",";
        }
        if(!rw_mut.isEmpty()) {
            outputLine = outputLine.substring(0, outputLine.length() - 1);
        }
        bw.write(outputLine);
        bw.newLine();
    }

    /**
     * Writes one block of 4 lines per read in respective output file.
     * @param bw
     * @param count
     * @param readSeq
     */
    private static void writeFastq(BufferedWriter bw, int count, String readSeq) throws IOException {
        bw.write("@" + count);
        bw.newLine();
        bw.write(readSeq);
        bw.newLine();
        bw.write("+" + count);
        bw.newLine();
        for(int i = 0; i < readSeq.length(); i++) {
            bw.write("I");
        }
        bw.newLine();
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
     * Retrieves the corresponding sequence for t post-splicing from fasta file and index.
     * @param t
     * @return t's sequence
     */
    private static String getSequence(Transcript t, String chromosome) throws IOException {
        //! NOTE:   this is based on new line characters being counted when stepping through bc of RandomAccessFile!

        String[] info = fidxMap.get(chromosome);
        long lineLengthWithNL = Long.parseLong(info[4]);
        long lineLength = Long.parseLong(info[3]);
        long trueEntryStart = Long.parseLong(info[2]);

        // find chromosome and transcript:

        long nextEntryStart = trueEntryStart;
        long lineStart;
        long inThisLine;

        StringBuilder tsb = new StringBuilder();
        for(Region exon : t.exons.regions) {
            int b;
            byte c;
            lineStart = ((exon.x1 - 1)/lineLength)*lineLengthWithNL;
            inThisLine = (exon.x1 - 1)%lineLength;
            nextEntryStart = trueEntryStart + lineStart + inThisLine;
            raf.seek(nextEntryStart);

            int exonLength = exon.x2 - exon.x1 + 1; // end exclusive! (true length is -1)

            String line = raf.readLine();
            int restLineLength = line.length();
            if(exonLength <= restLineLength) {
                tsb.append(line.substring(0, exonLength));
            } else {
                // read first rest of line
                tsb.append(line);
                int restExonLength = exonLength - restLineLength;
                // read full lines as long as possible
                while(restExonLength >= lineLength) {
                    line = raf.readLine();
                    tsb.append(line);
                    restExonLength -= lineLength;
                }
                line = raf.readLine();
                if(line != null && (exonLength > 0)) {
                    // read last rest of line
                    tsb.append(line.substring(0, restExonLength));
                }
            }
        }
        String transcrSeq = tsb.toString();
        return transcrSeq;
    }

    /**
     * Calculates the reverse complement of a (genomic) sequence.
     * @param s
     * @return reverse complement of s
     */
    private static String reverseComplement(String s) {
        // reverse string
        String revs = new StringBuilder(s).reverse().toString();
        // calculate complement
        StringBuilder rc = new StringBuilder();
        for(int i = 0; i < revs.length(); i++) {
            if(revs.charAt(i) == 'A') {
                rc.append("T");
            } else if(revs.charAt(i) == 'C') {
                rc.append("G");
            } else if(revs.charAt(i) == 'T') {
                rc.append("A");
            } else if(revs.charAt(i) == 'G') {
                rc.append("C");
            } else if(revs.charAt(i) == 'N') {
                rc.append("N");
            } else {
                System.err.println("Found a non-genomic conform base...");
            }
        }
        return rc.toString();
    }

    /**
     * Inserts mutations in a string to a given rate.
     * @param s
     * @param rate
     * @return
     */
    private static String mutate(String s, float rate) {
        ArrayList<Integer> mutPositions = new ArrayList<>();
        double r1;
        double r2;
        StringBuilder sb = new StringBuilder(s);
        for(int i = 0; i < s.length(); i++) {
            r1 = Math.random()*100;
            if(r1 <= rate) {
                // randomly mutate
                char c = s.charAt(i);
                char cNew;
                do {
                    r2 = Math.random();
                    if(r2 < 1f/4f) {
                        cNew = 'A';
                    } else if(r2 < 1f/2f) {
                        cNew = 'C';
                    } else if(r2 < 3f/4f) {
                        cNew = 'T';
                    } else {
                        cNew = 'G';
                    }
                } while(c == cNew);
                sb.setCharAt(i, cNew);
                mutPositions.add(i);    // 0 based !
            }
        }
        current_tMutPos = mutPositions;
        return sb.toString();
    }

    /**
     * For a given number of transcript based indices, the program checks if the given read includes them and transfers
     * them into read based indices.
     * @return read based indices within the given read
     */
    private static ArrayList<Integer> transPosToReadPos(ArrayList<Integer> indices, int startPos, int readLength) {
        int indexNew = -1;
        ArrayList<Integer> positionsInReads = new ArrayList<>();
        for(Integer index : indices) {
            if(index >= startPos && index <= (startPos + readLength - 1)) { // 0 based
                indexNew = index - startPos;    // if the mutation positions on reads shall be 0 based
                positionsInReads.add(indexNew);
            }
        }
        return positionsInReads;
    }

    /**
     * For a given number of transcript based indices, the program checks if the given REVERSE-read includes them and
     * transfers them into REVERSE-read based indices.
     * @param indices
     * @param startPos
     * @param readLength
     * @return read based indices within the given REVERSE read
     */
    private static ArrayList<Integer> revTransPosToReadPos(ArrayList<Integer> indices, int startPos, int readLength) {
        int indexNew = -1;
        ArrayList<Integer> positionsInReadsDesc = new ArrayList<>();
        int fragEnd = startPos + readLength - 1;
        for(Integer index : indices) {
            if(index >= startPos && index <= fragEnd) { // 0 based
                indexNew = fragEnd - index;    // if the mutation positions on reads shall be 0 based
                positionsInReadsDesc.add(indexNew);
            }
        }
        // change to smaller indices first
        ArrayList<Integer> positionsInReadsAsc = new ArrayList<>();
        for(int i = positionsInReadsDesc.size() - 1; i > -1; i--) {
            positionsInReadsAsc.add(positionsInReadsDesc.get(i));
        }
        return positionsInReadsAsc;
    }

    /**
     * Returns the genome-based region(s) for a given read, defined by length and starting position on transcript t.
     * @param t
     * @param readLength
     * @param startPos
     * @return
     */
    private static RegionVector recreateGenome(Transcript t, int readLength, int startPos) {
        RegionVector regVec = new RegionVector();

        // 8.1) find start of read in genome

        int toBeMatched = startPos; // number of bases we still have to "go over" (if 0 based: without - 1)
                                    // if fragment starts on 0-based pos 4, we still have to go over 4 bases
        int trueStart = -1;
        int exonLength;
        int restOfReadLength = readLength;

        for(Region exon : t.exons.regions) {
            exonLength = exon.x2 - exon.x1 + 1; // end inclusive? based
            if(trueStart == -1) {
                // still looking for start
                if(toBeMatched < exonLength) {
                    // start lies within this exon, only if its really smaller than the length
                    trueStart = exon.x1 + toBeMatched;

                    // 8.2) found start, find end of read in genome and build RegionVector of exons within read

                    if((trueStart + readLength - 1) <= exon.x2) {   // end-inclusive based ofc
                        // single exon spanning read, region vector contains only 1 region
                        regVec.regions.add(new Region(trueStart, (trueStart + readLength - 1), t.id));  // add endincl.!
                        break;
                    } else {
                        // multiple exon spanning read, region vector contains at least 2 regions
                        regVec.regions.add(new Region(trueStart, exon.x2, t.id));
                        restOfReadLength = readLength - (exon.x2 - trueStart + 1);
                    }

                } else {
                    toBeMatched -= exonLength;
                }
            } else {
                // 8.3) find end for multi exon spanning read

                if((exon.x1 + restOfReadLength - 1) <= exon.x2) {
                    // end is in this exon
                    regVec.regions.add(new Region(exon.x1, (exon.x1 + restOfReadLength - 1), t.id));
                    break;
                } else {
                    // read spans whole exon and end is in another
                    regVec.regions.add(exon);
                    restOfReadLength = restOfReadLength - exonLength;
                }
            }
        }
        return regVec;
    }

    /**
     * Reads the .fai file and converts it into a usable HashMap with chromosome names as keys.
     * @return the HashMap
     * @throws IOException
     */
    private static HashMap<String, String[]> readFastaIndex() throws IOException {
        raf = new RandomAccessFile(fasta, "r");
        fidxMap = new HashMap<>();
        based = new HashMap<>();
        rankInFidx = new HashMap<>();
        BufferedReader brFidx = new BufferedReader(new FileReader(fidx));
        String line;
        String[] splitLine = new String[5];
        int rank = 0;

        while((line = brFidx.readLine()) != null) {
            splitLine = line.split("\t");
            based.put(splitLine[0], Long.parseLong(splitLine[3]) + 1); // calculate offset = headerlength + newline
            splitLine = line.split("\t");
            fidxMap.put(splitLine[0], splitLine);
            rankInFidx.put(splitLine[0], rank);
            rank++;
        }
        return fidxMap;
    }

    /**
     * Prints usage when program gets invoked with false parameters.
     *
     */
    private static void printUsage() {
        // parameters
        System.out.println("Following parameters (in any order) are necessary:\n" +
                "-length (int): read length\n" +
                "-frlength (int), -SD (int): fragment length distribution\n" +
                "-readcounts (file name): table containing simulation info (simul.readcons)\n" +
                "-mutationrate (float): rate of mutations\n" +
                "-fasta (file name), -fidx (file name), -gtf (file name): file names necessary for computation\n" +
                "-od (path): output directory for the 3 output files");
        // what will be done
        System.out.println("If the parameters are correct, the program will simulate all reads for mapping. " +
                "The sense-strand read sequences can be found fw.fastq, antisense-strand read sequences in rw.fastq. " +
                "Index and additional information for all reads are stored in read.mappinginfo.");
    }

    public static void main(String[] args) throws IOException {
        /*
        Input:
        - length (int): read length
        - frlength (int), SD (int): fragment length distribution
        - readcounts: simul.readcons, table
        - mutationrate (float): mutationsrate
        - fasta (file)
        - fidx (file)
        - gtf (file)
        - od (pwd): output directory
        */

        // Argparse
        try {
            float sd = -1;
            float mean = -1;
            int readLength = -1;
            float mt = -1;
            String outputDir = "";
            if(args.length == 0) {
                printUsage();
            } else if(args.length == 18){
                for(int i = 0; i < args.length - 1; i += 2) {
                    if(args[i].equals("-length")) {
                        readLength = Integer.parseInt(args[i + 1]);
                    } else if(args[i].equals("-frlength")) {
                        mean = Float.parseFloat(args[i + 1]);
                    } else if(args[i].equals("-SD")) {
                        sd = Float.parseFloat(args[i + 1]);
                    } else if(args[i].equals("-readcounts")) {
                        readcons = new File(args[i + 1]);
                    } else if(args[i].equals("-mutationrate")) {
                        mt = Float.parseFloat(args[i + 1]);
                    } else if(args[i].equals("-fasta")) {
                        fasta = new File(args[i + 1]);
                    } else if(args[i].equals("-fidx")) {
                        fidx = new File(args[i + 1]);
                    } else if(args[i].equals("-gtf")) {
                        gtf = new File(args[i + 1]);
                    } else if(args[i].equals("-od")) {
                        outputDir = args[i + 1];
                    } else {
                        printUsage();
                    }
                }
                simulateReads(mean, sd, outputDir, readLength, mt);
            } else {
                printUsage();
            }
        } catch (IOException e) {
            System.err.println("Error: Something went wrong.");
        }

        /*
        gtf = new File("Homo_sapiens.GRCh37.75.gtf");
        fasta = new File("Homo_sapiens.GRCh37.75.dna.toplevel.fa");
        fidx = new File("Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai");
        readcons = new File("simul.readcons");
        simulateReads(20, 80, ".", 75, 1.0f);
         */
    }
}
