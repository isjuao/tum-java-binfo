import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

public class RegionVector {
    public TreeSet<Region> regions; // !! different from ExonSkipping (HashSet) !!
    // public HashSet<String> ids; // should be *transcriptID_exonNumber* if there is no exons

    public RegionVector() {
        this.regions = new TreeSet<>();
        // this.ids = new HashSet<>();
    }

    /**
     * Calculates the total length of a RV.
     * @return length
     */
    public int getLength() {
        int length = 0;
        for(Region reg : regions) {
            length += reg.x2 - reg.x1 + 1;
        }
        return length;
    }

    /**
     * Calculates the introns between the regions of the RegionVector.
     *
     * @return a RegionVector containing the introns.
     */
    public RegionVector getTntrons() {
        TreeSet<Region> introns = new TreeSet<>();
        int endOfLastCDS = 0;
        int counter = 0;
        for (Region reg : this.regions) {
            if (counter == 0) {
                // start of transcript
                endOfLastCDS = reg.x2;
            } else {
                introns.add(new Region(endOfLastCDS + 1, reg.x1 - 1, ""));  // -1 oder nicht?
                endOfLastCDS = reg.x1;
            }
            counter++;
        }
        RegionVector re = new RegionVector();
        re.regions = introns;
        return re;
    }

    /**
     * Decides whether a given RegionVector is contained in the one this method is called upon.
     *
     * @param regVec
     * @return True if the given RegionVector is contained in this RegionVector.
     */
    public boolean contains(RegionVector regVec) {
        for (Region other : regVec.regions) {
            int otherStart = other.x1;
            int otherEnd = other.x2;
            boolean foundFittingMy = false;
            for (Region my : this.regions) {
                if (my.x2 >= otherEnd) {
                    if (my.x1 <= otherStart) {
                        // found a fitting my region!
                        foundFittingMy = true;
                        break;
                    } else {
                        // my region does not conserve left end of other region
                        return false;
                    }
                }
            }
            if (!foundFittingMy) {
                // all my regions did not conserve right end of other region
                return false;
            }
        }
        return true;
    }

    /**
     * Decides whether a given RegionVector is a sub-RV of the one this method is called upon.
     *
     * @param regVec
     * @return True if the given RegionVector is a sub-RV of this RegionVector.
     */
    public boolean superRV(RegionVector regVec) {
        int firstStart = regVec.regions.first().x1;
        int lastEnd = regVec.regions.last().x2;
        RegionVector cutThisRV = this.cut(firstStart, lastEnd);
        if (cutThisRV != null) {
            boolean areEqual = cutThisRV.isEqual(regVec);
            return areEqual;
        } else {
            return false;
        }
    }

    /**
     * Cuts the current RegionVector according to the given indices (including).
     *
     * @param start
     * @param end
     * @return The cut RegionVector, or null if start or end were not in a Region.
     */
    public RegionVector cut(int start, int end) {
        boolean searchForEnd = false;
        RegionVector regVecNew = new RegionVector();
        boolean firstRegion = true;
        for (Region reg : this.regions) {
            if (!searchForEnd) {
                if (reg.x2 >= start) {
                    if (reg.x1 <= start) {
                        // found first Region
                        if (end <= reg.x2) {
                            // case: [ start end ]
                            regVecNew.regions.add(new Region(start, end, reg.belongsTo));
                            return regVecNew;
                        } else {
                            // case: [ start ] ...end?
                            regVecNew.regions.add(new Region(start, reg.x2, reg.belongsTo));
                            searchForEnd = true;
                        }
                    } else {
                        if (firstRegion) {
                            // special case
                            if (end <= reg.x2) {
                                // case: start [ end ]
                                regVecNew.regions.add(new Region(reg.x1, end, reg.belongsTo));
                                return regVecNew;
                            } else {
                                // case: start [ ] ...end?
                                regVecNew.regions.add(reg);
                                searchForEnd = true;
                            }
                        } else {
                            // start is not in a Region!
                            return null;
                        }
                    }
                }
            } else {
                if (reg.x2 >= end) {
                    if (reg.x1 <= end) {
                        // found second Region
                        regVecNew.regions.add(new Region(reg.x1, end, reg.belongsTo));
                        return regVecNew;
                    } else {
                        // end is not in a Region!
                        return null;
                    }
                } else {
                    // add in between Region
                    regVecNew.regions.add(reg);
                }
            }
            firstRegion = false;
        }
        if (searchForEnd) {
            return regVecNew;
        }
        return null;
    }

    public RegionVector cut2(int start, int end) {
        boolean searchForEnd = false;
        RegionVector regVecNew = new RegionVector();
        boolean firstRegion = true;
        for (Region reg : this.regions) {
            if (!searchForEnd) {
                if (reg.x2 >= start) {
                    if (reg.x1 <= start) {
                        // found first Region
                        if (end <= reg.x2) {
                            // case: [ start end ]
                            regVecNew.regions.add(new Region(start, end, reg.belongsTo));
                            return regVecNew;
                        } else {
                            // case: [ start ] ...end?
                            regVecNew.regions.add(new Region(start, reg.x2, reg.belongsTo));
                            searchForEnd = true;
                        }
                    } else {
                        if (firstRegion) {
                            // special case
                            if (end <= reg.x2) {
                                // case: start [ end ]
                                regVecNew.regions.add(new Region(reg.x1, end, reg.belongsTo));
                                return regVecNew;
                            } else {
                                // case: start [ ] ...end?
                                regVecNew.regions.add(reg);
                                searchForEnd = true;
                            }
                        } else {
                            // start is not in a Region!
                            return null;
                        }
                    }
                }
            } else {
                if (reg.x2 >= end) {
                    if (reg.x1 <= end) {
                        // found second Region
                        regVecNew.regions.add(new Region(reg.x1, end, reg.belongsTo));
                        return regVecNew;
                    } else {
                        // end is not in a Region!
                        return null;
                    }
                } else {
                    // add in between Region
                    regVecNew.regions.add(reg);
                }
            }
            firstRegion = false;
        }
        if (searchForEnd) {
            return regVecNew;
        }
        return null;
    }

    /**
     * Decides whether a given RegionVector has exactly the same Regions as this one and is therefore equal.
     *
     * @param regVec
     * @return True if the given RegionVector is completely equal to this RegionVector.
     */
    public boolean isEqual(RegionVector regVec) {
        Iterator<Region> it = this.regions.iterator();
        for (Region r : regVec.regions) {
            if (it.hasNext()) {
                Region next = it.next();
                if (r.x1 != next.x1 || r.x2 != next.x2) {
                    return false;
                }
            } else {
                // other RegionVector has Regions this one does not have
                return false;
            }
        }
        if (it.hasNext()) {
            // this RegionVector has Regions this one does not have
            return false;
        }
        return true;
    }

    public boolean equals(Object other) {
        // throw new RuntimeException("Plz use isEquals bc of mem");
        return isEqual((RegionVector)other);
    }

    /**
     * Ligates two neighboring or overlapping Regions in this RegionVector to one Region.
     *
     * @return the transformed version of this RegionVector
     */
    public RegionVector ligateNeighbors() {
        Region left = null;
        RegionVector transformed = new RegionVector();
        for (Region r : this.regions) {
            if (left == null) {
                transformed.regions.add(r);
                left = r;
                continue;
            } else {
                // check if left neighbor is adjoining or overlapping
                if (left.x2 >= (r.x1 - 1)) {
                    Region rNew = new Region(left.x1, Math.max(left.x2, r.x2), r.belongsTo);
                    transformed.regions.remove(left);   // remove old single Region
                    transformed.regions.add(rNew);      // add new merged Region
                    left = rNew;
                } else {
                    transformed.regions.add(r);         // add new single Region
                    left = r;
                }
            }
        }
        return transformed;
    }

    @Override
    public int hashCode() {
        int result = 0;
        for(Region r : regions) {
            result ^= r.hashCode();
        }
        return(result);
    }
}
