import java.io.*;
import java.util.ArrayList;
import java.util.TreeMap;

public class Analyze {

    /**
     * Collects information about the fragment length distribution in a read.mappinginfo file.
     * @throws IOException
     */
    private static void flDistribution() throws IOException {
        File readMapInfo = new File("output/read.mappinginfo");
        BufferedReader br = new BufferedReader(new FileReader(readMapInfo));

        TreeMap<Integer, Integer> flCounts = new TreeMap<>();
        String line = br.readLine();
        int start;
        int end;
        int length;
        String[] splitLine;
        File outputFile = new File("analysis/fl-distribution.txt");
        outputFile.createNewFile();
        FileOutputStream fileOutStream = new FileOutputStream(outputFile);
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fileOutStream));

        // collect fragment length for each line
        while((line = br.readLine()) != null) {
            splitLine = line.split("\t");
            start = Integer.parseInt(splitLine[4].split("-")[0]);
            end = Integer.parseInt(splitLine[5].split("-")[1]);
            length = end - start;
            if(flCounts.containsKey(length)) {
                flCounts.put(length, flCounts.get(length) + 1);
            } else {
                flCounts.put(length, 1);
            }
        }

        // write distribution in output file
        int old = flCounts.firstKey();
        for(Integer l : flCounts.keySet()) {
            fillIn(old, l, bw);
            bw.write(l + "\t" + flCounts.get(l));
            bw.newLine();
            old = l;
        }
        bw.close();
    }

    /**
     * Collects information about the distribution of the number of mutations.
     */
    private static void mutDistribution() throws IOException {
        File readMapInfo = new File("output/read.mappinginfo");
        BufferedReader br = new BufferedReader(new FileReader(readMapInfo));

        TreeMap<Integer, Integer> mutCounts = new TreeMap<>();
        String line = br.readLine();
        String[] splitLine;
        String[] fwMutLine;
        String[] rwMutLine;
        int fwMut;
        int rwMut;
        File outputFile = new File("analysis/mut-distribution.txt");
        outputFile.createNewFile();
        FileOutputStream fileOutStream = new FileOutputStream(outputFile);
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fileOutStream));

        // collect number of mutations
        while((line = br.readLine()) != null) {
            splitLine = line.split("\t");
            if(splitLine.length >= 9 && splitLine[8] != null) {
                fwMutLine = splitLine[8].split(",");
                fwMut = fwMutLine.length;
            } else {
                fwMut = 0;
            }
            if(splitLine.length >= 10 && splitLine[9] != null) {
                rwMutLine = splitLine[9].split(",");
                rwMut = rwMutLine.length;
            } else {
                rwMut = 0;
            }

            if(mutCounts.containsKey(fwMut)) {
                mutCounts.put(fwMut, mutCounts.get(fwMut) + 1);
            } else {
                mutCounts.put(fwMut, 1);
            }
            if(mutCounts.containsKey(rwMut)) {
                mutCounts.put(rwMut, mutCounts.get(rwMut) + 1);
            } else {
                mutCounts.put(rwMut, 1);
            }
        }

        // write distribution in output file
        int old = 0;
        for(Integer c : mutCounts.keySet()) {
            fillIn(old, c, bw);
            bw.write(c + "\t" + mutCounts.get(c));
            bw.newLine();
            old = c;
        }
        bw.close();
    }

    /**
     * Fills in 0 values in files for plotting.
     * @param old
     * @param x
     * @param bw
     * @throws IOException
     */
    private static void fillIn(int old, int x, BufferedWriter bw) throws IOException {
        if(old + 1 != x) {
            for(int i = old + 1; i < x; i++) {
                bw.write(i + "\t" + 0);
                bw.newLine();
            }
        }
    }

    public static void main(String[] args) throws IOException {
        mutDistribution();
    }
}
