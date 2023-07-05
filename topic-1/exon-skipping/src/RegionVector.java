import java.util.*;

public class RegionVector {
    public HashSet<Region> regions; // should only allow unique entries
    public HashSet<String> ids; // should be *transcriptID_exonNumber* if there is no exons
                                                    // for introns, ID is transcriptID_x1_x2
    public RegionVector() {
        this.regions = new HashSet<>();
        this.ids = new HashSet<>();
    }

    /**
     * Sorts coordinates, only possible because no overlaps.
     * @return
     */
    public TreeMap<Integer, Integer> getCoordinatesSorted() {
        TreeMap<Integer, Integer> sortedRegions = new TreeMap<>();
        for(Region reg : regions) {
            sortedRegions.put(reg.x1, reg.x2);
        }
        return sortedRegions;
    }
}
