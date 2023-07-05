import java.util.HashMap;
import java.util.HashSet;

public class Gene {
    public Region region;
    public String id;
    public boolean strand;
    public HashMap<String, RegionVector> transcripts;
    public String biotype;

    public Gene(Region region, String id, boolean strand) {
        this.region = region;
        this.id = id;
        this.strand = strand;
        transcripts = new HashMap<>();
    }

    /**
     * Alternative constructor for "empty" genes, where only the id (and later biotype) is required.
     * @param id
     */
    public Gene(String id) {
        this.id = id;
    }

    /**
     * Calculates the merged RegionVector of all transcripts.
     * @return the merged RegionVector of all transcripts
     */
    public RegionVector mergedTranscripts() {
        RegionVector transcript;
        RegionVector mergedTranscripts = new RegionVector();
        for(String transcriptID : this.transcripts.keySet()) {
            transcript = this.transcripts.get(transcriptID);
            mergedTranscripts.regions.addAll(transcript.regions);
        }
        mergedTranscripts = mergedTranscripts.ligateNeighbors();
        return mergedTranscripts;
    }
}
