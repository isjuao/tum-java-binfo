import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;


public class Test {
    private static void checkPcrIndex(String annotFileName) throws IOException {
        String line;
        String[] splitLine;
        BufferedReader br = new BufferedReader(new FileReader(new File(annotFileName)));
        String firstColumn;
        String lastColumn;

        File out = new File(("pcr-" + annotFileName));
        out.delete();
        FileOutputStream outStream = new FileOutputStream(out);
        BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(outStream));

        while((line = br.readLine()) != null) {
            splitLine = line.split("\t");
            if(splitLine.length <= 2) {
                bwOut.write(line);
                bwOut.newLine();
            } else {
                firstColumn = splitLine[0];
                lastColumn = splitLine[splitLine.length - 1];
                bwOut.write((firstColumn + "\t" + lastColumn));
                bwOut.newLine();
            }
        }
        bwOut.close();
    }

    private static void scanBams(String bamFileName) throws IOException {
        SAMFileReader samReader = new SAMFileReader(new File(bamFileName), false);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> it = samReader.iterator();
        int readcounter = 0;
        int unmappedcounter = 0;
        int mateunmappedcounter = 0;
        int unpairedcounter = 0;

        HashMap<Integer, Integer> qualities = new HashMap<>();
        while(it.hasNext()) {
            readcounter++;
            SAMRecord sr = it.next();
            boolean strand = !sr.getReadNegativeStrandFlag();
            boolean mateStrand = !sr.getMateNegativeStrandFlag();
            String name = sr.getReadName();

            int mapq = sr.getMappingQuality();
            String test = sr.getBaseQualityString();
            if(qualities.containsKey(mapq)) {
                int count = qualities.get(mapq);
                qualities.replace(mapq, count + 1);
            } else {
                qualities.put(mapq, 1);
            }

            /*if(sr.getReadUnmappedFlag()) {
                unmappedcounter++;
            }
            if(sr.getMateUnmappedFlag()) {
                mateunmappedcounter++;
            }
            if(!sr.getReadPairedFlag()) {
                unpairedcounter++;
            }
            System.out.println(ref);
            boolean unmapped = sr.getReadUnmappedFlag();
            boolean duplicated = sr.getDuplicateReadFlag();
            boolean mt;
            if(ref.equals("MT")) {
                mt = true;
            } else {
                mt = false;
            }
            if(unmapped || duplicated || mt) {
                System.out.println(sr.getReadName());
            }*/
        }
        for(int mapq : qualities.keySet()) {
            System.out.println(mapq + "\t" + qualities.get(mapq));
        }
    }

    public static void main(String[] args) throws IOException {
        // checkPcrIndex("h.sp.8.annot");
        // (!) maybe add "bams/"
        String fileName;
        // scanBams("shifted_GL1_esp.bam");
        // scanBams("shifted_GL1.bam");
        // scanBams("sub_sam.sam");
        scanBams("filtered_alignment.bam");
        // scanBams("shift2_sorted_GL1.bam");
        // scanBams("sort2_filter_GL1.bam");
        // scanBams("shifted_sort2_filter_GL1.bam");
        // scanBams("GL1.bam");
        /*
        for(int i = 1; i < 11; i++) {
            fileName = "h.sp." + i + ".bam";
            scanBams(fileName);
        }
        for(int i = 1; i < 11; i++) {
            fileName = "h.sn." + i + ".bam";
            scanBams(fileName);
        }
        for(int i = 1; i < 11; i++) {
            fileName = "y.ns." + i + ".bam";
            scanBams(fileName);
        }
        */
    }
}
