import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

public class Gene {
    public String id;
    public String chromosome; // seqname field
    public char strand;
    public String name;     // geneName attribute
    public HashMap<String, Transcript> transcriptsAll; // collection of all transcripts
    public int numberOfCDStanscr;
    public RegionVector intronsAll;
    public HashSet<String> intronsWithoutRegions;

    public Gene() {
        this.transcriptsAll = new HashMap<>();
    }

    /**
     * Calculates a list of all unique introns for the gene, depending on what transcript they belong to.
     */
    public void calcIntrons() {
        for(Transcript currentTranscript : transcriptsAll.values()) {
            if(currentTranscript.containsCDS) {
                int endOfLastCDS = 0;
                TreeMap<Integer, Integer> sortedCDSs = currentTranscript.cdss.getCoordinatesSorted();
                int counter = 0;
                TreeMap<Integer, Integer> sortedIntrons = new TreeMap<>();
                for(Integer currentCDSstart : sortedCDSs.keySet()) {
                    if(counter == 0) {
                        // start of transcript
                        endOfLastCDS = sortedCDSs.get(currentCDSstart);
                    } else {
                        //TODO: -1 oder nicht
                        sortedIntrons.put(endOfLastCDS + 1, currentCDSstart);
                        endOfLastCDS = sortedCDSs.get(currentCDSstart);
                    }
                    counter++;
                }
                currentTranscript.introns = sortedIntrons;
            }
        }
    }

    /**
     *  Collects all previously calculated introns (per transcript) into a collection (per gene).
     */
    public void collectIntrons() {
        this.intronsAll = new RegionVector();
        this.intronsWithoutRegions = new HashSet<>();
        Integer intronEnd = -1;
        String newID = "blub";
        for(Transcript currentTranscript : transcriptsAll.values()) {
            if(currentTranscript.containsCDS) {
                for(Integer intronStart : currentTranscript.introns.keySet()) {
                    intronEnd = currentTranscript.introns.get(intronStart);
                    // another way of testing to see if introns are duplicates
                    if(intronsWithoutRegions.contains("" + intronStart + intronEnd)) {
                        System.err.println("found duplicate!");
                    }
                    Region possiblyNewIntron = new Region(intronStart, intronEnd, currentTranscript.id);
                    if(!intronsAll.regions.contains(possiblyNewIntron)) {
                        intronsAll.regions.add(possiblyNewIntron);
                        newID = currentTranscript.id + "_" + intronStart + "_" + intronEnd;
                        intronsAll.ids.add(newID);
                    }
                }
            }
        }
        int x = 0;  // for breakpoint
    }
}
