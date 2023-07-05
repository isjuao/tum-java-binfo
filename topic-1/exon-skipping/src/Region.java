public class Region implements Comparable<Region> {
    public int x1;
    public int x2;
    public String belongsTo; // for transcripts: geneID, for CDSs: transcriptID

    public Region(int x1, int x2, String belongsTo) {
        this.x1 = x1;
        this.x2 = x2;
        this.belongsTo = belongsTo;
    }

    /**
     * Checks if two regions are equal, indepent of the transcript they belong to.
     * @param region to compare
     * @return boolean
     */
    @Override
    public boolean equals(Object region) {
        if(region instanceof Region) {
            if(this.x1 == ((Region) region).x1 && this.x2 == ((Region) region).x2
                    // && this.belongsTo.equals(((Region) region).belongsTo)
            ) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    /**
     * Forces HashCode of two equal Regions (same x1/x2) to be the same.
     * @return a funny number
     */
    @Override
    public int hashCode() {
        int upperBits = (x1 & 0xffff) << 16;
        int lowerBits = (x2 & 0xffff);
        return upperBits | lowerBits;
    }

    /**
     * Orders Regions so that Region may be used in a TreeSet.
     * @param region
     * @return
     */
    @Override
    public int compareTo(Region region) {
        if(region.x1 < this.x1) {
            return 1;
        } else if(region.x1 == this.x1) {
            return 1;
        } else {
            return -1;
        }
    }
}
