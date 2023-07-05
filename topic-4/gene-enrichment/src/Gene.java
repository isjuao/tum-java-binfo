import java.util.HashSet;

public class Gene {
    public String id;
    public float foldChange;
    public Boolean significant; // = null, if not measured!!
    public HashSet<GoTerm> belongsTo;

    public Gene(String id) {
        this.id = id;
        foldChange = 0;       // default
        significant = null;    // default
        belongsTo = new HashSet<>();
    }

    /**
     * To make sure that equals does what it should.
     * @param other
     * @return true if the other Gene is equal to this gene (if their ids match 100%).
     */
    public boolean equals(Object other) {
        return ((Gene)other).id.equals(this.id);
    }
}
