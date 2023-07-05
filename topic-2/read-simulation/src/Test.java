import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.io.*;

public class Test {
    private static File gtf = new File("Homo_sapiens.GRCh37.75.gtf");
    private static HashMap<String, String[]> fidxMap;
    private static File fidx = new File("Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai");
    private static RandomAccessFile raf;
    private static File fasta = new File("Homo_sapiens.GRCh37.75.dna.toplevel.fa");
    private static HashMap<String, Integer> rankInFidx;
    private static HashMap<String, Long> based;      // possible offset
    private static HashMap<String, Boolean> refStrands;
    private static HashMap<String, Boolean> myStrands;
    private static String temp;

    private static void testFasta2() throws IOException {
        File fasta = new File("Homo_sapiens.GRCh37.75.dna.toplevel.fa");
        BufferedReader brFasta = new BufferedReader(new FileReader(fasta));
        String line;
        int counter = 0;
        int linecounter = 0;
        char c;
        while ((line = brFasta.readLine()) != null) {
            linecounter++;
            if (line.startsWith(">")) {
                System.out.println(line.charAt(1) + " " + line.length());
            }
            for(int i = 0; i < line.length(); i++) {
                counter++;
                c = line.charAt(i);
                if(counter == 253404855 || counter == 253404856 || counter == 253404854) {
                    c = line.charAt(i);
                    System.out.println(c);
                }
                if(linecounter == 4154177) {
                    System.out.println(line);
                }
            }
        }
    }

    private static void testFasta1() throws IOException {
        File fasta = new File("Homo_sapiens.GRCh37.75.dna.toplevel.fa");
        BufferedReader brFasta = new BufferedReader(new FileReader(fasta));
        String line;
        int linecounter = 0;
        int linesLength = -1;
        int linesToStart;
        int trueEntryStart;
        while ((line = brFasta.readLine()) != null) {
            linecounter++;
            int l = line.length();
            if (line.startsWith(">")) {
                trueEntryStart = (253404911 - 56*(1 + 1))/61*60
                        + 55*(1 + 1);
                /*
                trueEntryStart = (500657663 - 56*(2 + 1))/61*60
                        + 55*(2 + 1);
                trueEntryStart = (56 - 56*(0 + 1))/61*60
                        + 55*(0 + 1);
                */
                linesToStart = (int)Math.ceil(((double)trueEntryStart)/60);
                linesLength = linesLength;
                linesLength = linesToStart + 1 + (int)Math.ceil(((double)243199373)/60);
                System.out.println(line.charAt(1) + " " + line.length());
            }
        }
    }

    private static void testFasta3() throws IOException {
        File fasta = new File("Homo_sapiens.GRCh37.75.dna.toplevel.fa");
        RandomAccessFile raf = new RandomAccessFile(fasta, "r");
        raf.seek(5);
        System.out.println(raf.readLine());
        System.out.println(raf.readLine());
    }

    private static void toLines(String s) throws IOException {
        // String r = "";
        BufferedWriter bw = new BufferedWriter(new FileWriter("test1.txt"));
        for(int i = 0; i < s.length(); i++) {
            bw.write(s.charAt(i));
            bw.newLine();
            // r += s.charAt(i) + "\n";
        }
        bw.close();
        // return r;
    }

    private static HashMap<String, String> testGenomeExtractor(String wishChr) throws IOException {
        HashMap<String, String> transcripts = new HashMap<>();
        refStrands = new HashMap<>();
        GZIPInputStream inputStream = new GZIPInputStream(
                new FileInputStream("Homo_sapiens.GRCh37.75.cdna.all.fa.gz"));
        BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
        String line;
        String id = "";
        StringBuilder seq = new StringBuilder();
        String[] splitLine;
        boolean wishChrCurrent = false;
        boolean oldTranscrExists = false;
        while((line = br.readLine()) != null) {
            if(line.startsWith(">")) {
                // found a header
                if(oldTranscrExists) {
                    // save old transcript
                    transcripts.put(id, seq.toString());
                }
                // prepare new transcript
                splitLine = line.split(" ");
                if(splitLine[2].split(":")[2].equals(wishChr)) {
                    wishChrCurrent = true;
                    id = splitLine[0].substring(1);
                    seq = new StringBuilder();
                    if(splitLine[2].charAt(splitLine[2].length() - 2) != '-') {
                        refStrands.put(id, true);
                    } else {
                        refStrands.put(id, false);
                    }
                    oldTranscrExists = true;
                } else {
                    wishChrCurrent = false;
                }
            } else if(wishChrCurrent){
                seq.append(line);
            }
        }
        if(oldTranscrExists) {
            // save old transcript
            transcripts.put(id, seq.toString());
        }
        return transcripts;
    }

    private static HashMap<String, String> myGenomeExtractor(String wishChr) throws IOException {
        fidxMap = readFastaIndex();
        myStrands = new HashMap<>();
        HashMap<String, String> transcripts = new HashMap<>();
        BufferedReader brGtf = new BufferedReader(new FileReader(gtf));
        Transcript t = null;
        String transcriptID = "blub";
        String chromosome = "blub";

        String lineGtf;

        String[] splitLineGtf;
        String[] attributes;
        boolean oldTranscrExists = false;

        while((lineGtf = brGtf.readLine()) != null) {
            if(lineGtf.startsWith("#")) continue;
            splitLineGtf = lineGtf.split("\t");
            if(splitLineGtf[2].equals("transcript") && splitLineGtf[0].equals(wishChr)) {
                if(oldTranscrExists) {
                    // DEAL WITH OLD TRANSCRIPT
                    // ------------------------

                    // set start and end of transcript according to exons

                    t.start = t.exons.regions.first().x1;
                    // need to find end manually
                    int biggestEnd = 0;
                    for (Region r : t.exons.regions) {
                        if (r.x2 > biggestEnd) {
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

                    transcripts.put(transcriptID, seq);
                }

                // DEAL WITH NEW TRANSCRIPT
                // ------------------------

                oldTranscrExists = true;
                attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                // found transcript
                t = new Transcript((transcriptID = findArgument("transcript_id", attributes)));
                if(splitLineGtf[6].charAt(0) == '+') {
                    myStrands.put(transcriptID, true);
                } else {
                    myStrands.put(transcriptID, false);
                }
                chromosome = splitLineGtf[0];
            }  else if(splitLineGtf[2].equals("exon") && t != null && splitLineGtf[0].equals(wishChr)) {
                attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                if(findArgument("transcript_id", attributes).equals(transcriptID)) {
                    // found CDS of transcript
                    t.exons.regions.add(new Region(Integer.parseInt(splitLineGtf[3]),
                            Integer.parseInt(splitLineGtf[4]), transcriptID));
                }
            }
        }
        if(oldTranscrExists) {
            // DEAL WITH OLD TRANSCRIPT
            // ------------------------

            // set start and end of transcript according to exons

            t.start = t.exons.regions.first().x1;
            // need to find end manually
            int biggestEnd = 0;
            for (Region r : t.exons.regions) {
                if (r.x2 > biggestEnd) {
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

            transcripts.put(transcriptID, seq);
        }

        return transcripts;
    }

    /**
     * Retrieves the corresponding sequence for t post-splicing from fasta file and index
     * @param t
     * @return t's sequence
     */
    private static String getSequence(Transcript t, String chromosome) throws IOException {
        //! NOTE:   this is based on new line characters being counted when stepping through bc of RandomAccessFile!

        if(t.id.equals("ENST00000605759")) {
            t = t;
        }

        String[] info = fidxMap.get(chromosome);
        long lineLengthWithNL = Long.parseLong(info[4]);
        long lineLength = Long.parseLong(info[3]);
        long trueEntryStart = Long.parseLong(info[2]);
        /*
        //! NOTE:   this is based on new line characters not being counted when stepping through, but being included in
        //          entry start indices
        long trueEntryStart = (Long.parseLong(info[2]) - 56*(rankInFidx.get(chromosome) + 1))/61*60
                + 55*(rankInFidx.get(chromosome) + 1);  // maybe need rounding?
        long linesToStart = (long)Math.ceil(((double)trueEntryStart)/60)
                + (rankInFidx.get(chromosome));     // line of HEADER of chromosome
        long trueEntryEnd; // evtl to do
        long linesLength = (long)Math.ceil(((double)Long.parseLong(info[2]))/60);
        long linesEnd = linesToStart + 1 +
                (long)Math.ceil(((double)Long.parseLong(info[2]))/60); // line of LAST LINE of chromosome */

        // find chromosome and transcript:

        long nextEntryStart = trueEntryStart;
        long lineStart;
        long inThisLine;
        //String line;
        StringBuilder tsb = new StringBuilder();
        for(Region exon : t.exons.regions) {
            // nextEntryStart = trueEntryStart + exon.x1 - 1;  // -1 bc fasta file is 0 based
            int b;
            byte c;
            lineStart = ((exon.x1 - 1)/lineLength)*lineLengthWithNL;
            inThisLine = (exon.x1 - 1)%lineLength;
            nextEntryStart = trueEntryStart + lineStart + inThisLine;
            raf.seek(nextEntryStart);
            /*for(int n = 0; n < exon.x1 - 1; n++) {  // end inclusive bc it is true start
                b = raf.read();
                if(b == -1) {
                    System.err.println("File error.");
                } else {
                    c = (byte)b;
                    if(c == '\n') {
                        // skip the character that shall not be parsed
                        n--;
                    }
                }
            } */
            // raf.seek(nextEntryStart);
            // line = raf.readLine();
            int exonLength = exon.x2 - exon.x1 + 1; // end exclusive! (true length is -1)
            /*
            int i = 0;
            for(int n = 0; n < exonLength; n++) {
                if(i >= line.length()) {
                    // need to look at another line
                    line = raf.readLine();
                    i = 0;
                } else {
                    transcrSeq += line.charAt(i);
                    i++;
                }
            }
            */
            String line = raf.readLine();
            int restLineLength = line.length();
            if(exonLength <= restLineLength) {
                tsb.append(line.substring(0, exonLength)); //TODO indices...
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
            /*
            for(int n = 0; n < exonLength; n++) {
                b = raf.read();
                if(b == -1) {
                    System.err.println("File error.");
                } else {
                    c = (byte)b;
                    if(c == '\n') {
                        // skip the character that shall not be parsed
                        n--;
                    } else {
                        // append valid base
                        tsb.append((char)c);
                    }
                }
            } */
        }
        String transcrSeq = tsb.toString();
        return transcrSeq;
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

    private static void compareTranscriptMaps(HashMap<String, String> ref, HashMap<String, String> my) {
        String mySeq;
        String refSeq;
        int numberMissing = 0;
        int numberWrong = 0;
        for(String refID : ref.keySet()) {
            mySeq = my.get(refID);
            if(mySeq == null) {
                System.err.println(refID + " missing in my map.");
                numberMissing++;
            } else {
                refSeq = ref.get(refID);
                if(!refSeq.equals(mySeq)) {
                    System.err.println("Sequences for " + refID + "(" + refStrands.get(refID) + "/" + myStrands.get(refID) +
                            ") not equal.");
                    numberWrong++;
                    if(refStrands.get(refID) == myStrands.get(refID)) {
                        refSeq = refSeq;
                        mySeq = mySeq;
                        String compMySeq = reverseComplement(mySeq);
                        System.err.println("Same strand wrong.");
                    }
                }
            }
        }
        System.out.println("\tOut of " + my.size());
        System.out.println("\t" + numberMissing + " were missing in total.");
        System.out.println("\t" + numberWrong + " were wrong in total.");
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

    public static void main(String[] args) throws IOException {
        /*RandomAccessFile f = new RandomAccessFile("Homo_sapiens.GRCh37.75.dna.toplevel.fa", "r");
        f.seek(56 + 11869 - 1);
        var line = f.readLine();*/
        /*
        testFasta3();
        String s = "GTGTGTTGATTACTTAAGGGCAGAGAAACAGGATGAAAATGGGTAAATAAGGCCAAATGACTTGATTAACTTTTGTGGGAGAACAACTTCAAAGGCCCACATAATTTAACATGCTGGAAATTCTCCCCATGCAAACATTTCAGGAGTTGACTAAGATCAGAGATGTGGCTACTGCCACAGCAGCTCCCTCTGTCTGTGCATGTTAACTTAGGTAACAGATGAACTCACAGACCAGCCTTCTCTTATCCTAGCCATGTGGACAACTGGGTACAAAGCTCTTTATATATTAAAGCCCTTTGTAACAATCTGCTTGCCACTC";
        System.out.println(s.length());
        toLines("CTGAGTGACACAGCAAGACTCCATCTCAAAAAAAAAAAAAAAAAAAAGTTAAGTCTTAACAACACAAAGTGGTCCTGGCTGAATATGGTGAACTCCACATTCGTTAAAAAAGTACTGTCACCTCATGAATATCTTAGCTAGAACCTGCAGCTGGGAAGCAGAGGTTTTTATTGTTTCCCAGAGTGAAAAACGCAAATGAGACGCTATTTCGTCTAGTTGTGTGCACAACTCGAGAGGTACAGCTCTTGATCTGATAACATGACCAGCTCTCTGTCCTAGGTTGGTCTCTACCTTTTTGGTTAATTTAATAGGAGACTAAAAAACCAATGGAAGAAAGTAGCCTGCTATTAAGTCAACCAAAGAGTAGGACCTGCCTTCCGGCTCGAAAAGCCTGGGTCTGGCCCATCCCACCCTATGGTGAGCTGCCGTGTGAACCAATGCTCTGTCTGAAGGCGGCACATCCCACAGCTCCGGAGCTGGTCTATGCTGTCACTCCTCCCTGACATTCTGTGAAGCTCCGCAGTATCCAGCCTGTGCTAGGAATCACACTGCATCGTCCAAACAGAACCACTTCAACAGATGCTTGGAGGTTTCATTTTTCTCTGTGTTTTTACTACCTCCAGCTTTAACATAATGTGGGAAAACTTTACTGCACACCTGCAAAATAGCATCAGGTTAGTTACAAGAAACTCGGACAACTCACTCTTGCAAATCCCACACTACAATCAGCCAAAGGCCCTGAAAAAGCCCAAACTCAGGCAGCTGAGTTCAGCAATACACATAAGCCCTGAGTCCCAGGCCATATGAACAGCGAAACAGCCATGTCACCAGAGCTTGCTGTCCTCACCTTCCCTCGGAACCCCCACTCCCCACCCGTGCCAGGCCTTCCCCTTCCACAGCTCTTAAAAGGCTCTCCTCACCTGATTTGGTGACTCTAGACCGATTTGTAGGTCCGGGTTTCAGCTCATCCTGGTCTTTTTCTTCCTTCTTAATGTCTTCAGGCTTTTTTTCCTCCTTCTTCTCAATCTTCTCTTCCTTCTTAATTACAGTTCTATTTAAAGACATGGTTCTAACAGGAAGCACAACCATTCTTAAGACCGCAGCCCTCTTTAAAAACCCTAACGACGCCCAATCCCAACTGCTGACTTAAGAATTTAACGCACACACACATGCCAACACACAGAAGGCTGAAAATCTGGACAGCAAAAACTGAGTTTCCTTACTCAAGGAATAATCATGCCCATGTCAGTTATGAACTGAAATGTCTTTTGATGATCAGATCTCATCCAAGACGGCAGTCTCTAAAATCCTTGAAACAGTCTTCAATGAGCAGAAGACCAATGAAAAGGAGATAGCGTCTGGTCTAAAACTGAGAAAAAACAAACCTCATTTTCTGATCAGCGCCCTCACTGGTTGAGGACTCTTTAGGAGCCGGAGGGACTTCATTACACGAGTATTTTCTCATTTTCTACGGATTTAACAGGAATAGAATGGCACAAGGTTTAATTACCAGGGAAATAAAGCAATCACAAATGTGGAACTGGCACGGCTGCCACACTAACACGAAGAGAACTGCTAGCTACACAGAGATTGGAAAACCCGCGGGTACACAAAGTTATACTGGTTACACTACCTGCGCCGCGACCGCCTCTAGGTAAGCCTCTACGCTACACGGAAGAGAAAAGGCCCCCACAGATTTAAACACCTACAACCACGCGGCTGGTTCTCAATAGTCTTCTTCGAGTTTTTGTTCAAGTCTGGGTCTTCTG");
        toLines("CTGAGTGACACAGCAAGACTCCATCTCAAAAAAAAAAAAAAAAAAAAGTTAAGTCTTAACAACACAAAGTGGTCCTGGCTGAATATGGTGAACTCCACATTCGTTAAAAAAGTACTGTCACCTCATGAATATCTTAGCTAGAACCTGCAGCTGGGAAGCAGAGGTTTTTATTGTTTCCCAGAGTGAAAAACGCAAATGAGACGCTATTTCGTCTAGTTGTGTGCACAACTCGAGAGGTACAGCTCTTGATCTGATAACATGACCAGCTCTCTGTCCTAGGTTGGTCTCTACCTTTTTGGTTAAGTTAATAGGAGACTAAAAAACCAATGGAAGAAAGTAGCCTGCTATTAAGTCAACCACAGAGCAGGACCTGCCTTCCGGCTCGAAAAGCCTGGGTCTGGCCCATCCCACCCTATGGTGAGCTGCCGTGTGAACCAATGCTCTGTCTGAAGGCGGCACATCCCACAGCTCCGGAGCTGGTCTTTGCTGTCACTCCTCCCTGACATTCTGTGAAGCTCCGCAGTATCCAGCCTGTGCTAGGAATCGCACTGCATCGTCCAAACAGAACCACTTCAACAGATGCTTGGAGGTTTCATTTTTCTCTGTGTTTTTACTACCTCCAGCTTTAACATAATGTGGGAAAACTTTACTGCACACCTGCAAAAAAGCATCAGTTTAGTTACAAGAAACTCGGACAACTCACTCTTGCAAATCCCACACTACAATCAGCCAAAGACCCTGAAAAAGCCCAAACTCAGGCAGCTGAGTTCAGCAATTCACATAAGCCCTGAGTCCCAGGCCATATGAAGAGCGAAACAGCCATGTCACCAGAGATTGCTGTCCTCACCTTCCCTCGAAACCCCCACTCCCCACCCGTGCCAGGCCTTCCCCTTCCACAGCTCTTAAAAGGCTCTCCTCACCTGATTTGGTGACTCTAGACCGATTTGTAGGTCCGGGTTTCAGCTCTTCCTGGTCTTTTTCTTCCTTCTCAATGTCTTCAGGCTTTTTTTCCTCCTTCTTCTCAATCTTCTCTTCCTTCTTAATTACAGTTCTATTTAAAGACATGGTTCTAACAGGAAGCACAACCATTCTTAAGACCGCAGCCCTCTTTAAAAACCCTAACGACGCCCAATCCCAACTGCTGACTAAAGAATTTAACGCACACACACATGCCAACACACAGAAGGCTGAAAATCTGGACAGCAAAAACTGAGTTTCCTTACTCAAGGAATAATCATGCCCATGTCAGTTATGAACTGAAATGTCTTTTGATGATCAGATCTCATCCAAGACGGTAGTCTCTAAAATCCTTGAAACAGTCTTCAATGAGCAGAAGACCAATGAAAAGGAGATACCGTCTGGTCTAAAACTGAGAAAAAACAAACCTCATTTTCTGATCAGCGCCCTCACTGGTTGAGGACTCTTTAGGAGCCGGAGGGACTTCATTACACGAGTATTTTCTCATTTTCTACGGATTTAACAGGAATAGAATGGCACAAGGTTTAATTACCAGGGAAATAAAGCAATCACAAATGTGGAACTGGCACGGCTGCCACACTAACACGAAGAGAACTGCTAGCTACACAGAGATTGGAAAACCCGCGGGTACACAAAGTTATCCTGGTTACACTACCTGCGCCGCGACCGCCTCTAGGTAAGCCTCTACGCTACACGGAAGAGAAAAAGCCCCCACAGATTTAAACACCTACAACCACGCGGCTGGTTCTCAATAGTCTTCTTCGAGTTTTTGTTCAAGTCTGGGTCTTCTG");
        while(true) {
            org.apache.commons.math3.distribution.NormalDistribution nd = new NormalDistribution(20, (double) 80);
            double sampleProp = nd.sample();
            int fl = Math.max(75, (int) sampleProp);
            fl = fl;
        }
        */
        temp = reverseComplement("AGATGATCCCAATTTTGTTACAACATCGAAAGCATCATAATCAGGAGCAAGTCGAACATATGCCTTGTTCTCTTTATCAGGACAAATCAGGGTGGTGACCTTGGCCACATCACTGTCATAGAGCTTCTTCACAGCCTGTCTGATCTGGTGCTTGTTGGCTTTAACATCCACAGTGAACACAAGCGTGTTGTTTTCTTCTATCTTCTTCCGGCCGACTCAGTGGTCAGCGGAAACTTGATGATAGCATAGTGGCCAAGCTTGTTTCTCCTGGGGGTGCTCTTCCGAGGATATCTGGGCTGCCTCCGGAGTCGCAGTGTCTTGGGCCGCCTGAAGGTGGGTGACATGCGGATCTTCTTTTTTGCGTGTGGCTGCGGACACCTTTCAACACTGCCTTCTTGGCCTTTAAAACCTTCACTTTGGCTTCGGCTTTAGGAGGAGCAGGAGCTTCCTTCGC");
        if(temp.equals("GCGAAGGAAGCTCCTGCTCCTCCTAAAGCCGAAGCCAAAGTGAAGGTTTTAAAGGCCAAGAAGGCAGTGTTGAAAGGTGTCCGCAGCCACACGCAAAAAAGAAGATCCGCATGTCACCCACCTTCAGGCGGCCCAAGACACTGCGACTCCGGAGGCAGCCCAGATATCCTCGGAAGAGCACCCCCAGGAGAAACAAGCTTGGCCACTATGCTATCATCAAGTTTCCGCTGACCACTGAGTCGGCCGGAAGAAGATAGAAGAAAACAACACGCTTGTGTTCACTGTGGATGTTAAAGCCAACAAGCACCAGATCAGACAGGCTGTGAAGAAGCTCTATGACAGTGATGTGGCCAAGGTCACCACCCTGATTTGTCCTGATAAAGAGAACAAGGCATATGTTCGACTTGCTCCTGATTATGATGCTTTCGATGTTGTAACAAAATTGGGATCATCT")){
            System.out.println("meh.");
        }
        String[] a = {"X", "Y", "MT"};
        for(String i : a) {
            System.out.println("WishChr: " + i);
            HashMap<String, String> transcriptsRef = testGenomeExtractor("" + i);
            HashMap<String, String> transcriptsMy = myGenomeExtractor("" + i);
            transcriptsMy = transcriptsMy;
            compareTranscriptMaps(transcriptsRef, transcriptsMy);
        }
    }
}
