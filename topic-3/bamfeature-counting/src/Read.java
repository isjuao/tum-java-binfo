import net.sf.samtools.AlignmentBlock;

import java.util.List;

public class Read {
    public String readID;
    public String chromosome;
    public RegionVector regions;

    public List<AlignmentBlock> blocks;     // ?! introns ?!
    public boolean strand;

    public Read(String readID, String chromosome, boolean strand) {
        this.readID = readID;
        this.chromosome = chromosome;
        this.strand = strand;
    }
}
