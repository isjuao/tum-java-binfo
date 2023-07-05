import java.util.TreeMap;

public class Transcript {
    public RegionVector exons;
    public String id;
    public int start;
    public int end;
    // public TreeMap<Integer, Integer> introns;
    String proteinID;

    public Transcript(String id) {
        this.id = id;
        exons = new RegionVector();
    }
}