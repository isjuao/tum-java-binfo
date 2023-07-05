import java.util.TreeMap;
import java.util.TreeSet;

public class Transcript {
    public RegionVector cdss;
    public RegionVector exons;
    public String id;
    public TreeMap<Integer, Integer> introns;
    boolean containsCDS;
    String proteinID;

    public Transcript(String id) {
        this.id = id;
        cdss = new RegionVector();
        exons = new RegionVector();
        containsCDS = false;
    }
}
