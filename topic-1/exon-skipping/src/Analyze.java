import javax.management.relation.Role;
import java.io.*;
import java.util.TreeMap;

public class Analyze {
    private static TreeMap<Integer, Integer> allMaxSkippedExons;
    private static TreeMap<Integer, Integer> allMaxSkippedBases;
    private static TreeMap<Integer, String> top10Exons;

    private static void analyzeMaxSkipped(String fileName, String outputFilePath) throws IOException {
        allMaxSkippedExons = new TreeMap<>();
        allMaxSkippedBases = new TreeMap<>();

        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = br.readLine();

        File outputFile = new File(outputFilePath);
        outputFile.createNewFile();
        FileOutputStream fileOutStream = new FileOutputStream(outputFile);
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fileOutStream));

        String[] splitLine;
        int currentMaxExons;
        int currentMaxBases;

        while((line = br.readLine()) != null) {
            splitLine = line.split("\t");
            currentMaxExons = Integer.parseInt(splitLine[10]);
            currentMaxBases = Integer.parseInt(splitLine[12]);
            if(allMaxSkippedExons.containsKey(currentMaxExons)) {
                allMaxSkippedExons.put(currentMaxExons, allMaxSkippedExons.get(currentMaxExons) + 1);
            } else {
                allMaxSkippedExons.put(currentMaxExons, 1);
            }
            if(allMaxSkippedBases.containsKey(currentMaxBases)) {
                allMaxSkippedBases.put(currentMaxBases, allMaxSkippedBases.get(currentMaxBases) + 1);
            } else {
                allMaxSkippedBases.put(currentMaxBases, 1);
            }
        }
        int old = 0;
        for(int maxSkippedExon : allMaxSkippedExons.keySet()) {
            fillIn(old, maxSkippedExon, bw);
            bw.write(maxSkippedExon + "\t" + allMaxSkippedExons.get(maxSkippedExon));
            bw.newLine();
            old = maxSkippedExon;
        }
        bw.write("#");
        old = 0;
        for(int maxSkippedBase : allMaxSkippedBases.keySet()) {
            fillIn(old, maxSkippedBase, bw);
            bw.write(maxSkippedBase + "\t" + allMaxSkippedBases.get(maxSkippedBase));
            bw.newLine();
            old = maxSkippedBase;
        }
        bw.close();
    }

    private static void analyzeTop10(String fileName, String outputFilePath) throws IOException {
        top10Exons = new TreeMap<>();

        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = br.readLine();

        File outputFile = new File(outputFilePath);
        outputFile.createNewFile();
        FileOutputStream fileOutStream = new FileOutputStream(outputFile);
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fileOutStream));

        String[] splitLine;
        int currentMaxExons;
        String currentGeneID;

        while((line = br.readLine()) != null) {
            splitLine = line.split("\t");
            currentMaxExons = Integer.parseInt(splitLine[12]);
            currentGeneID = splitLine[0];
            if(!top10Exons.values().isEmpty()) {
                if(top10Exons.containsValue(currentGeneID)) {
                    // gene is already in top 10, but maybe needs updating
                    for(Integer score : top10Exons.keySet()) {
                        if(top10Exons.get(score).equals(currentGeneID) && score < currentMaxExons) {
                            top10Exons.remove(score);
                            top10Exons.put(currentMaxExons, currentGeneID);
                            break;
                        }
                    }
                } else if(top10Exons.size() == 10) {
                    // gene is not yet in top 10, but map is full
                    if(currentMaxExons > top10Exons.descendingKeySet().last()) {
                        // gene should belong in top 10!
                        top10Exons.remove(top10Exons.descendingKeySet().last());
                        top10Exons.put(currentMaxExons, currentGeneID);
                    }
                } else {
                    // gene is not yet in top 10 and map is not full yet
                    top10Exons.put(currentMaxExons, currentGeneID);
                }
            } else {
                // map is empty
                top10Exons.put(currentMaxExons, currentGeneID);
            }
        }

        for(Integer score : top10Exons.descendingKeySet()) {
            bw.write((score + "\t" + top10Exons.get(score)));
            bw.newLine();
        }
        bw.close();
    }

    private static void fillIn(int old, int end, BufferedWriter bw) throws IOException {
        if(old + 1 != end) {
            for(int i = old + 1; i < end; i++) {
                bw.write(i + "\t" + 0);
                bw.newLine();
            }
        }
    }

    public static void main(String[] args) throws IOException {
        analyzeTop10("out-Homo-38.90.txt", "blub.txt");
    }
}
