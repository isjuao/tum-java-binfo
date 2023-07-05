import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import java.util.HashMap;
import java.util.HashSet;

public class GoTerm {
    public String id;
    public String name;
    public String[] parentIDs;
    public GoTerm[] parents;
    public HashSet<Gene> geneSet;
    public double[] allFoldChanges;
    boolean simulTrue;
    public double hgP;
    public double fejP;
    public double kP;
    public double kStat;
    public double hgPCorr;
    public double fejPCorr;
    public double kPCorr;
    boolean discardThis;

    public GoTerm(String id, String name, String[] parentIDs) {
        this.id = id;
        this.name = name;
        this.parentIDs = parentIDs;
        this.parents = new GoTerm[0];
        this.geneSet = new HashSet<>();
        this.simulTrue = false;
        this.hgP = -1.0d;
        this.hgPCorr = -1.0d;
        this.fejP = -1.0d;
        this.fejPCorr = -1.0d;
        this.kP = -1.0d;
        this.kPCorr = -1.0d;
        this.kStat = -1.0d;
        discardThis = false;
    }

    /**
     * Recursively adds given Gene to self and parents.
     * @param gene
     */
    public void addGene(Gene gene) {
        for(GoTerm parent : parents) {
            parent.addGene(gene);
        }
        geneSet.add(gene);  // automatically only adds if not contained already
        gene.belongsTo.add(this);
    }

    public double getHgP() {
        return this.hgP;
    }

    public double getFejP() {
        return this.fejP;
    }

    public double getkP() {
        return this.kP;
    }

    /**
     * Counts the number of measured associated genes and the number of significant genes of the node.
     * @return an array of the number of genes in the gene set that are measured ([0]), signifcantly measured ([1]).
     */
    public int[] getMeasuredSize() {
        int n = 0;
        int s = 0;
        for(Gene associatedGene : geneSet) {
            if(associatedGene.significant != null) {
                // measured Gene
                n++;
                if(associatedGene.significant) {
                    s++;
                }
            }
        }
        int[] result = {n, s};
        return result;
    }

    /**
     * Calculates and sets the enrichment p-values given by the hypergeometric distribution
     * and the Fisher's Exact test using jackknifing
     * @param N population size (total number of genes (intersect))
     * @param K number of successes (total number of DE/significant genes)
     * @param n gene set sized that are also measured
     * @param k number of observed successes (this term's DE/significant genes = noverlap (small intersect))
     */
    public void calcOraPValues(int N, int K, int n, int k) {
        if(n == 0) {
            int x = 0;  // breakpoint
        }
        // hypergeometric based
        HypergeometricDistribution hg = new HypergeometricDistribution(N, K, n);
        this.hgP = hg.upperCumulativeProbability(k);
        // Fisher's Exact with jackknifing based
        HypergeometricDistribution fej = new HypergeometricDistribution(
                N - 1, K - 1, n - 1);
        this.fejP = fej.upperCumulativeProbability(k - 1);
    }

    /**
     * Collects the fold changes of all genes of the node.
     */
    public void collectFoldChanges() {
        allFoldChanges = new double[geneSet.size()];
        int i = 0;
        for(Gene g : geneSet) {
            allFoldChanges[i] = g.foldChange;
            i++;
        }
    }

    /**
     * Calculates and sets the distribution comparison based Kolmogorov-Smirnov values.
     */
    public void calcDistPValues(HashMap<String, Gene> allGenes, int N, int numMeasured) {
        KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
        double[] bgDistrib = new double[N - numMeasured];
        double[] inSetDistrib = new double[numMeasured];

        // collect in-set distribution (measured only!)
        int j = 0;
        for(Gene currentGene : geneSet) {
            if(currentGene.significant != null) {
                inSetDistrib[j] = currentGene.foldChange;
                j++;
            }
        }
        // collect background distribution
        int i = 0;
        for(Gene currentGene : allGenes.values()) {
            boolean isBackground = !geneSet.contains(currentGene);
            /*inner: for(Gene associatedGene : geneSet) {
                if(associatedGene.id.equals(gID)) {
                    isBackground = false;
                    break inner;
                }
            }*/
            if(isBackground && (currentGene.significant != null)) {
                bgDistrib[i] = currentGene.foldChange;
                i++;
            } else {
                // either not background or not measured
                if(!isBackground) {
                    int x = 0;  // breakpoint
                }
            }
        }

        // calculate statistics
        if(inSetDistrib.length == 1) {
            int x = 4;  // breakpoint
        }
        this.kP = ks.kolmogorovSmirnovTest(inSetDistrib, bgDistrib);
        // allFoldChanges oder inSetDistrib (Filter auf measured)?
        this.kStat = ks.kolmogorovSmirnovStatistic(inSetDistrib, bgDistrib);
    }

    /**
     * Decides if given GoTerm is an Ascendant of this one.
     * @param other
     * @return true if this one is a Descendant of the given one.
     */
    public boolean isDescendantOf(GoTerm other) {
        boolean result = false;
        for(GoTerm parent : this.parents) {
            if(parent.id.equals(other.id)) {
                return true;
            } else {
                result |= parent.isDescendantOf(other);
            }
        }
        return result;
    }

    /**
     * Calculates the length of the path between this node and its possibly ascendant other.
     * @param other
     * @param pathLength
     * @return -1, if other is not an ascendant of this node, else the length of the path between the 2 relatives.
     */
    public int getPathLengthToAscendant(GoTerm other, int pathLength) {
        int result = -1;
        for(GoTerm parent : this.parents) {
            if(parent.id.equals(other.id)) {
                if(result != -1) {
                    return Math.min(result, pathLength + 1);
                }
                return pathLength + 1;
            } else {
                int possPathLength = parent.getPathLengthToAscendant(other, pathLength + 1);
                if(result == -1) {
                    result = possPathLength;
                    // result = Math.min(result, possPathLength); WRONG??
                } else {
                    if(possPathLength != -1) {
                        result = Math.min(result, possPathLength);
                    }
                }
            }
        }
        return result;
    }

    /**
     * Calculates the path between this node and its possibly ascendant other.
     * @param other
     * @param path
     * @return null, if other is not an ascendant of this node, else the path of this node to its ascendant.
     */
    public String getPathToAscendant3(GoTerm other, String path) {
        String result = null;
        if(this.id.equals(other.id)) {
            result = path + "|" + this.id;
        } else {
            int minPathLength = Integer.MAX_VALUE;
            for(GoTerm parent : this.parents) {
                String testPath = parent.getPathToAscendant3(other, path + "|" + this.id);
                if(testPath != null) {
                    int testPathLength = testPath.length();
                    if(testPathLength < minPathLength) {
                        minPathLength = testPathLength;
                        result = testPath;
                    }
                }
            }
        }
        return result;
    }

    public String toString() {
        String parentIDStrings = "";
        for(String parentID : parentIDs) {
            parentIDStrings += parentID + ", ";
        }
        String fullName = id + ": " + parentIDStrings;
        return fullName;
    }
}
