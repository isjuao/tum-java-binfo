import augmentedTree.*;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.*;
import java.util.*;

public class BamFeatures {
    private static File bamFile;
    private static File gtfFile;
    private static Boolean frstrand = null;
    private static HashMap<String, Read> readsPerChromosome;    // unsure if possible bc of memory
    private static HashMap<String, SAMRecord> recordsPerChrom;
    private static IntervalTree<Region> genesFw;   // fw genes
    private static IntervalTree<Region> genesRw;    // rw genes
    private static HashMap<String, String> biotypes;
    private static String outputTsvString;
    private static boolean nonDuplicates;

    private static void iterateBamFile() throws IOException {
        SAMFileReader samReader = new SAMFileReader(bamFile, false);
        samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> it = samReader.iterator();

        File out = new File(outputTsvString);
        out.delete();
        FileOutputStream outStream = new FileOutputStream(out);
        BufferedWriter bwOut = new BufferedWriter(new OutputStreamWriter(outStream));

        // write additional info file for plots
        /*File infoOut;
        if(!nonDuplicates) {
            infoOut = new File("info_" + outputTsvString);
        } else {
            infoOut = new File("info_nd_" + outputTsvString);
        }
        infoOut.delete();
        FileOutputStream infoOutStream = new FileOutputStream(infoOut);
        BufferedWriter bwInfoOut = new BufferedWriter(new OutputStreamWriter(infoOutStream)); */

        String readID = "";
        String chromosome = "";
        boolean firstOfPair;
        Boolean queryStrand = null;    // true = forward, false = reverse
        Boolean leadingStrand = null;
        HashMap<String, Gene> genesFwMap = new HashMap<>();
        HashMap<String, Gene> genesRwMap = new HashMap<>();
        HashMap<String, Gene> genesBoth = new HashMap<>();
        HashMap<RegionVector, Integer> chrMergeRVsFw = new HashMap<>(); // if strand unspecific: use this only!
        HashMap<RegionVector, Integer> chrMergeRVsRw = new HashMap<>();
        HashMap<String, Integer> readCounts = new HashMap<>();
        HashMap<String, Integer> mergeRVlengths = new HashMap<>();
        long libSize = 0;

        while (it.hasNext()) {
            SAMRecord sr = it.next();

            // A) collect information
            // ----------------------

            readID = sr.getReadName();
            firstOfPair = sr.getFirstOfPairFlag();
            if(readID.equals("3412065")) {
                readID = readID;    // breakpoint
            }

            // B) strandness
            // -------------

            queryStrand = !sr.getReadNegativeStrandFlag();
            if(firstOfPair) {
                if(frstrand == null) {
                    // not strand specific
                } else if(frstrand) {
                    // first read maps sense to transcribed region
                    leadingStrand = queryStrand;
                } else {
                    // first read maps antisense to transcribed region
                    leadingStrand = !queryStrand;   // because OtherOfPair decides
                }
            } else {
                if(frstrand == null) {
                    // not strand specific
                } else if(frstrand) {
                    leadingStrand = !queryStrand;   // because OtherOfPair decides
                } else {
                    leadingStrand = queryStrand;
                }
            }

            // 0) clear lookup if chromosome changes
            // -------------------------------------

            if(!chromosome.equals(sr.getReferenceName())) {
                /*
                // collect readcounts
                for(String geneID : readCounts.keySet()) {
                    if(readCounts.size() != mergeRVlengths.size()) {
                        System.err.println("Unequal readCounts and mergeRVlengths sizes!");
                    }
                    int readCount = readCounts.get(geneID);
                    Integer mergeLength = mergeRVlengths.get(geneID);
                    if(mergeLength == null) {
                        System.err.println("Missing " + geneID + " in mergeRVlenghts!");
                    } else {
                        bwInfoOut.write(geneID + "\t" + readCount + "\t" + mergeLength);
                        bwInfoOut.newLine();
                    }
                } */

                /* DEBUG */
                /*for(String geneID : mergeRVlengths.keySet()) {
                    Integer readCount = readCounts.get(geneID);
                    if(readCount == null) {
                        System.err.println("Missing " + geneID + " in readCounts!");
                    }
                }*/

                chromosome = sr.getReferenceName();
                readsPerChromosome = new HashMap<>();
                recordsPerChrom = new HashMap<>();
                genesFw = loadGtfPerChromosome(sr.getReferenceName(), true);
                genesRw = loadGtfPerChromosome(sr.getReferenceName(), false);
                biotypes = new HashMap<>();
                genesBoth = new HashMap<>();
                chrMergeRVsFw = new HashMap<>();
                chrMergeRVsRw = new HashMap<>();
                readCounts = new HashMap<>();
                mergeRVlengths = new HashMap<>();

                // collect genes
                genesFwMap = generateGeneSet(genesFw, true);
                genesRwMap = generateGeneSet(genesRw, false);

                // collect transcripts
                genesFwMap = addTranscriptsToGene(genesFwMap, true, chromosome);
                genesRwMap = addTranscriptsToGene(genesRwMap, false, chromosome);

                genesBoth.putAll(genesRwMap);
                genesBoth.putAll(genesFwMap);

                if(genesBoth.containsKey("ENSG00000233653")) {
                    int z = 6;
                }
            }

            // C) decide: handle read pair or run preliminary checks?
            // ------------------------------------------------------

            Read otherRead = readsPerChromosome.get(readID);
            Read currentRead = new Read(readID, chromosome, queryStrand);
            currentRead.blocks = sr.getAlignmentBlocks();
            SAMRecord otherRecord = recordsPerChrom.get(readID);
            if(otherRead != null) {
                // C.1) handle readpair
                Read read1;
                Read read2;
                if(frstrand == null || firstOfPair) {   // not sure about first condition
                    read1 = currentRead;
                    read2 = otherRead;
                } else {
                    read1 = otherRead;
                    read2 = currentRead;
                }

                // ::: DEBUG :::
                String x = "x";
                if(readID.equals("7482309")) {
                    x = readID;
                }
                x = x.substring(0,1);
                // ::: DEBUG :::

                // D) split inconsistent check
                // ---------------------------

                // collect regions, introns, and the mergedIntrons RV
                read1.regions = blocksToRV(read1, new RegionVector()).ligateNeighbors();
                read2.regions = blocksToRV(read2, new RegionVector()).ligateNeighbors();
                RegionVector introns1 = read1.regions.getTntrons();
                RegionVector introns2 = read2.regions.getTntrons();
                RegionVector mergedIntrons;
                RegionVector mergedReadPair = new RegionVector();
                int splitCount = -1;
                int pcrIndex = 0;

                // calculate cgenes and igenes
                int alStart = Math.min(Math.min(sr.getAlignmentStart(), otherRecord.getAlignmentStart()),
                        Math.min(sr.getAlignmentEnd(), otherRecord.getAlignmentEnd()));
                int alEnd = Math.max(Math.max(sr.getAlignmentStart(), otherRecord.getAlignmentStart()),
                        Math.max(sr.getAlignmentEnd(), otherRecord.getAlignmentEnd()));
                // check if full genes are contained and no genes are spanning
                boolean doesContain = doesContainFlag(alStart, alEnd, leadingStrand);   // maybe use fopStrand????
                boolean noSpanning = isNotSpannedFlag(alStart, alEnd, leadingStrand);

                // test if inconsistent
                if(inconsistency3(read1.regions, read2.regions)) {
                    if(noSpanning && doesContain) {
                        continue;
                    }
                    // inconsistent -> skip all other annotations, directly output
                    bwOut.write(readID + "\tsplit-inconsistent:true");
                    bwOut.newLine();
                    // System.out.println(readID + "\tsplit-inconsistent:true");
                    continue;
                } else {
                    // calculate nsplits
                    mergedIntrons = introns1;
                    mergedIntrons.regions.addAll(introns2.regions);
                    mergedIntrons = mergedIntrons.ligateNeighbors();
                    splitCount = mergedIntrons.regions.size();

                    // merge readpair
                    mergedReadPair.regions.addAll(read1.regions.regions);
                    mergedReadPair.regions.addAll(read2.regions.regions);
                    mergedReadPair = mergedReadPair.ligateNeighbors();

                    if(mergedReadPair.regions.first().x1 == 58470 &&
                        mergedReadPair.regions.last().x2 == 58667) {
                        int u = 4;    // breakpoint
                    }

                    // store for pcr index

                    HashMap<RegionVector, Integer> chrMergeRVs;
                    if(frstrand == null) {
                        chrMergeRVs = chrMergeRVsFw;
                    } else {
                        if (leadingStrand) {
                            chrMergeRVs = chrMergeRVsFw;
                        } else {
                            chrMergeRVs = chrMergeRVsRw;
                        }
                    }
                    Integer pcrCount = chrMergeRVs.get(mergedReadPair);
                    if(pcrCount != null) {
                        chrMergeRVs.replace(mergedReadPair, pcrCount + 1);
                        pcrIndex = pcrCount + 1;
                    } else {
                        chrMergeRVs.put(mergedReadPair, 0);
                        pcrIndex = 0;
                    }
                }

                // E) intergenic and antisense check
                // ---------------------------------

                int minDistToGene = -1;
                Boolean antisenseFlag = false;
                String matchClassOutStr = "";
                int geneCount = 0;
                if(noSpanning) {
                    // G.1) check if intergenic
                    if(doesContain) {
                        // -> skip
                        continue;
                    } else {
                        // intergenic, collect distance
                        antisenseFlag = false;
                        if(frstrand == null) {
                            minDistToGene = Math.max(0, Math.min(distanceToGenes(alStart, alEnd, true) - 1,
                                    distanceToGenes(alStart, alEnd, false) - 1));
                        } else {
                            if(leadingStrand) {
                                minDistToGene = Math.max(0, distanceToGenes(alStart, alEnd, true) - 1);
                            } else {
                                minDistToGene = Math.max(0, distanceToGenes(alStart, alEnd, false) - 1);
                            }
                        }

                        // G.2) check if antisense
                        if(frstrand != null) {
                            if(leadingStrand) {
                                if(!checkWithinGenes(alStart, alEnd, false)) {
                                    antisenseFlag = true;
                                }
                            } else {
                                if(!checkWithinGenes(alStart, alEnd, true)) {
                                    antisenseFlag = true;
                                }
                            }
                        }
                    }
                }

                if(!antisenseFlag && minDistToGene == -1) {
                    // F) match plausability levels
                    // ----------------------------

                    // not intergenic/antisense -> find highest class

                    // collect genes & transcripts
                    StringBuilder matchClassOut = new StringBuilder();
                    HashMap<String, Gene> genes;
                    if(frstrand == null) {
                        genes = genesBoth;
                    } else {
                        if (leadingStrand) {
                            genes = genesFwMap;
                        } else {
                            genes = genesRwMap;
                        }
                    }

                    Gene currentGene;
                    HashMap<String, Gene> transcriptomic = new HashMap<>();
                    HashMap<String, Integer> potMergeRVlengths = new HashMap<>();   // if readpair is (m)transcriptomic
                    RegionVector mergedTranscripts;
                    HashMap<String, RegionVector> currentTranscripts;
                    for(String geneID : genes.keySet()) {
                        currentGene = genes.get(geneID);

                        // check if transcriptomic
                        RegionVector transcript;
                        currentTranscripts = new HashMap<>();
                        for(String transcriptID : currentGene.transcripts.keySet()) {
                            transcript = currentGene.transcripts.get(transcriptID);
                            if(transcriptID.equals("ENST00000432723")) {
                                transcript = transcript;    // breakpoint
                            }
                            if(transcript.superRV(read1.regions) && transcript.superRV(read2.regions)) {
                                // both reads are sub-RVs of the transcript
                                currentTranscripts.put(transcriptID, transcript);
                                /**NOTE: brauch eig nur transcript-ID*/
                            }
                        }
                        if(!currentTranscripts.isEmpty()) {
                            // the gene is transcriptomic!

                            // collect merge-RV lengths
                            if(!nonDuplicates || pcrIndex == 0) {
                                if (!mergeRVlengths.containsKey(geneID)) {
                                    // collect merge-RV lengths
                                    if (geneID.equals("YAL016C-B")) {
                                        int v = 5;  // breakpoint
                                    }
                                    mergedTranscripts = currentGene.mergedTranscripts();
                                    int mergedLength = mergedTranscripts.getLength();
                                    potMergeRVlengths.put(geneID, mergedLength);
                                }
                            }

                            Gene toAdd = new Gene(geneID);
                            toAdd.biotype = currentGene.biotype;
                            toAdd.transcripts = currentTranscripts;
                            transcriptomic.put(geneID, toAdd);
                        }
                    }
                    // check if the highest class is transcriptomic
                    if(transcriptomic.isEmpty()) {
                        HashMap<String, Gene> mergeTranscriptomic = new HashMap<>();
                        potMergeRVlengths = new HashMap<>();
                        for(String geneID : genes.keySet()) {
                            currentGene = genes.get(geneID);

                            // merge transcripts
                            mergedTranscripts = currentGene.mergedTranscripts();
                            int mergedLength = mergedTranscripts.getLength();

                            // check if merge-transcriptomic
                            if(mergedTranscripts.contains(mergedReadPair)) {
                                if(!nonDuplicates || pcrIndex == 0) {
                                    if (!mergeRVlengths.containsKey(geneID)) {
                                        // collect merge-RV lengths
                                        potMergeRVlengths.put(geneID, mergedLength);
                                    }
                                }

                                mergeTranscriptomic.put(geneID, currentGene);
                            }
                        }
                        if(mergeTranscriptomic.isEmpty()) {
                            // fetch intronic matches
                            HashSet<Region> spanningGenes = new HashSet<>();
                            if(frstrand == null) {
                                spanningGenes = genesFw.getIntervalsSpanning(alStart, alEnd, spanningGenes);
                                HashSet<Region> spaGenes2 = genesRw.getIntervalsSpanning(alStart, alEnd, spanningGenes);
                                spanningGenes.addAll(spaGenes2);
                            } else {
                                if(leadingStrand) {
                                    spanningGenes = genesFw.getIntervalsSpanning(alStart, alEnd, spanningGenes);
                                } else {
                                    spanningGenes = genesRw.getIntervalsSpanning(alStart, alEnd, spanningGenes);
                                }
                            }
                            // IS INTRONIC -> output:
                            for(Region geneRegion : spanningGenes) {
                                geneCount++;
                                String geneID = geneRegion.belongsTo;   // unsure if that works
                                String biotype = biotypes.get(geneID);
                                matchClassOut.append(geneID).append(",").append(biotype).append(":INTRON|");
                            }
                        } else {
                            // IS MERGE-TRANSCRIPTOMIC -> output:
                            if(!nonDuplicates || pcrIndex == 0) {
                                libSize++;
                                mergeRVlengths.putAll(potMergeRVlengths);   // collect merge-RV lengths
                            }
                            for(String geneID : mergeTranscriptomic.keySet()) {
                                if(!nonDuplicates || pcrIndex == 0) {
                                    // collect readcounts
                                    Integer readCountPerGene = readCounts.get(geneID);
                                    if (readCountPerGene != null) {
                                        readCounts.replace(geneID, readCountPerGene + 1);
                                    } else {
                                        readCounts.put(geneID, 1);
                                    }
                                }

                                geneCount++;
                                currentGene = mergeTranscriptomic.get(geneID);
                                matchClassOut.append(geneID).append(",").append(currentGene.biotype).append(":MERGED|");
                            }
                        }
                    } else {
                        // IS TRANSCRIPTOMIC -> output:
                        if(!nonDuplicates || pcrIndex == 0) {
                            libSize++;
                            mergeRVlengths.putAll(potMergeRVlengths);   // collect merge-RV lengths
                        }
                        for(String geneID : transcriptomic.keySet()) {
                            if(!nonDuplicates || pcrIndex == 0) {
                                // collect readcounts
                                Integer readCountPerGene = readCounts.get(geneID);
                                if (readCountPerGene != null) {
                                    readCounts.replace(geneID, readCountPerGene + 1);
                                } else {
                                    readCounts.put(geneID, 1);
                                }
                            }

                            geneCount++;
                            currentGene = transcriptomic.get(geneID);
                            matchClassOut.append(geneID).append(",").append(currentGene.biotype).append(":");
                            for(String transcriptID : currentGene.transcripts.keySet()) {
                                matchClassOut.append(transcriptID).append(",");
                            }
                            if(matchClassOut.charAt(matchClassOut.length() - 1) == ',') {
                                // cut redundant characters
                                matchClassOut = new StringBuilder(matchClassOut.substring(0,
                                        matchClassOut.length() - 1));
                            }
                            matchClassOut.append("|");
                        }
                    }
                    matchClassOutStr = matchClassOut.toString();
                    if(matchClassOutStr.charAt(matchClassOut.length() - 1) == '|') {
                        // cut redundant characters
                        matchClassOutStr = matchClassOutStr.substring(0, matchClassOut.length() - 1);
                    }
                }

                // G) clipping
                // -----------

                int clipSum = Math.abs(otherRecord.getAlignmentStart() - otherRecord.getUnclippedStart());
                clipSum += Math.abs(otherRecord.getAlignmentEnd() - otherRecord.getUnclippedEnd());
                clipSum += Math.abs(sr.getAlignmentStart() - sr.getUnclippedStart());
                clipSum += Math.abs(sr.getAlignmentEnd() - sr.getUnclippedEnd());

                // H) mismatches
                // -------------

                Integer nm = (Integer) sr.getAttribute("NM");   // different tagging options
                if (nm == null) nm = (Integer) sr.getAttribute("nM");
                if (nm == null) nm = (Integer) sr.getAttribute("XM");
                int mismatchSum = nm;
                nm = (Integer) otherRecord.getAttribute("NM");
                if (nm == null) nm = (Integer) otherRecord.getAttribute("nM");
                if (nm == null) nm = (Integer) otherRecord.getAttribute("XM");
                mismatchSum += nm;

                // I) output
                // ---------

                // writePcrOnly(bwOut, readID, pcrIndex);
                writeOutput(bwOut, readID, mismatchSum, clipSum, splitCount, geneCount,
                            matchClassOutStr, minDistToGene, antisenseFlag, pcrIndex);
            } else {
                // C.2) run preliminary checks

                if(sr.getNotPrimaryAlignmentFlag()) {   // getSupplementaryAlignmentFlag()?
                    // read is secondary -> skip
                    continue;
                }
                if(sr.getReadUnmappedFlag()) {
                    // read is unmapped -> skip
                    continue;
                }
                if(sr.getMateUnmappedFlag()) {
                    // mate is unmapped -> skip
                    continue;
                }
                if(!chromosome.equals(sr.getMateReferenceName())) {
                    // mate on different chromosome -> skip
                    continue;
                }
                if(sr.getReadNegativeStrandFlag() == sr.getMateNegativeStrandFlag()) {
                    // mate on the same strand -> skip
                    continue;
                }

                // pre-intergenic check: very unsure
                /*
                int start1 = sr.getAlignmentStart();        // start
                int start2 = sr.getMateAlignmentStart();    // start of mate
                if(frstrand == null) {
                    if((checkContainsGenes(start1, start2, true) && checkWithinGenes(start1, start2, true))
                            || (checkContainsGenes(start1, start2, false)
                            && checkWithinGenes(start1, start2, false))) {
                        // read is PRE-intergenic -> skip
                        continue;
                    }
                } else {
                    if(queryStrand) {   // right strand? incorporate firstofpair? what if secondofpair? ueberhaupt?
                        if(checkContainsGenes(start1, start2, true)
                                && checkWithinGenes(start1, start2, true)) {
                            // read is PRE-intergenic -> skip
                            continue;
                        }
                    } else {
                        if(checkContainsGenes(start1, start2, false)
                                && checkWithinGenes(start1, start2, false)) {
                            // read is PRE-intergenic -> skip
                            continue;
                        }
                    }
                } */

                // C.3) add new read
                currentRead.blocks = sr.getAlignmentBlocks();
                readsPerChromosome.put(readID, currentRead);
                recordsPerChrom.put(readID, sr);
            }
        }
        // collect readcounts for last chromosome genes
        /*
        for(String geneID : readCounts.keySet()) {
            if(readCounts.size() != mergeRVlengths.size()) {
                System.err.println("Unequal readCounts and mergeRVlengths sizes!");
            }
            int readCount = readCounts.get(geneID);
            Integer mergeLength = mergeRVlengths.get(geneID);
            if(mergeLength == null) {
                System.err.println("Missing " + geneID + " in mergeRVlenghts!");
            } else {
                bwInfoOut.write(geneID + "\t" + readCount + "\t" + mergeLength);
                bwInfoOut.newLine();
            }
        } */

        /* DEBUG */
        /*for(String geneID : mergeRVlengths.keySet()) {
            Integer readCount = readCounts.get(geneID);
            if(readCount == null) {
                System.err.println("Missing " + geneID + " in readCounts!");
            }
        }*/
        // collect libsize
        // bwInfoOut.write("#\t" + libSize);
        bwOut.close();
        // bwInfoOut.close();
    }

    private static void writePcrOnly(BufferedWriter bw, String readID, int pcrIndex) throws IOException {
        bw.write(readID + "\tpcrindex: " + pcrIndex);
        bw.newLine();
    }

    /**
     * Writes one line for the tabular output file of outputTsvString.
     * @param bw
     * @param readID
     * @param mm
     * @param clipping
     * @param splits
     * @param gcount
     * @param matchClassOutStr
     * @param gdist
     * @param antisense
     * @throws IOException
     */
    private static void writeOutput(BufferedWriter bw, String readID, int mm, int clipping, int splits,
                                    int gcount, String matchClassOutStr, int gdist, Boolean antisense, int pcrIndex)
                                    throws IOException {
        String outString = readID + "\tmm:" + mm + "\tclipping:" + clipping +
                            "\tnsplit:" + splits + "\tgcount:" + gcount;
        if(gdist == -1) {
            outString += "\t" + matchClassOutStr;
        } else {
            outString += "\tgdist:" + gdist + "\tantisense:" + antisense.toString();
        }
        outString += "\tpcrindex: " + pcrIndex;

        // System.out.println(outString);
        bw.write(outString);
        bw.newLine();
    }

    /**
     * Transforms a IntervalTree of genes as Regions into a HashMap of Genes that can contain transcripts.
     * @param geneTree
     * @param strand
     * @return HashMap of genes referenced by their ids.
     */
    private static HashMap<String, Gene> generateGeneSet(IntervalTree<Region> geneTree, boolean strand) {
        HashMap<String, Gene> genes = new HashMap<>();
        Gene g;
        for(Region reg : geneTree) {
            g = new Gene(reg, reg.belongsTo, strand);
            genes.put(g.id, g);
        }
        return genes;
    }

    /**
     * Collects (strand-specific) all transcripts for a given set of genes to their respective genes, using a gtf.
     * @param genes
     * @param strand
     */
    private static HashMap<String, Gene> addTranscriptsToGene(HashMap<String, Gene> genes, boolean strand, String chr)
            throws IOException {
        char s;
        if(strand) {
            s = '+';
        } else {
            s = '-';
        }
        int eStart = -1;
        int eEnd = -1;
        String lineGtf;
        String[] splitLineGtf;
        String[] attributes;
        String transcriptID = "blub";
        String geneID = "blub";
        Gene currentGene = null;
        RegionVector currentTranscript = null;
        boolean firstExons = true;
        boolean isFirstTranscript = true;
        boolean geneHandled = true;

        BufferedReader brGtf = new BufferedReader(new FileReader(gtfFile));
        while((lineGtf = brGtf.readLine()) != null) {
            if(lineGtf.startsWith("#")) continue;

            splitLineGtf = lineGtf.split("\t");
            if(splitLineGtf[0].equals(chr)) {
                if(splitLineGtf[2].equals("gene") && splitLineGtf[6].charAt(0) == s) {
                    if(currentGene != null) {
                        if(!firstExons) {
                            // store old transcript
                            currentGene.transcripts.put(transcriptID, currentTranscript);
                            firstExons = true;
                        }
                        // handle old gene
                        genes.remove(geneID);
                        genes.put(geneID, currentGene);
                        geneHandled = true;
                    }

                    // handle new gene
                    attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                    geneID = findArgument("gene_id", attributes);
                    currentGene = genes.get(geneID);
                    currentGene.biotype = splitLineGtf[1];
                    biotypes.put(geneID, splitLineGtf[1]);
                    isFirstTranscript = true;
                    geneHandled = false;

                    currentTranscript = new RegionVector();
                }
                if(splitLineGtf[2].equals("exon") && splitLineGtf[6].charAt(0) == s) {
                    // found an exon of a transcript
                    attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                    String transcriptIDnew = findArgument("transcript_id", attributes);
                    if(transcriptIDnew.equals("YJR159W")) {
                        int c = 5;  // breakpoint
                    }
                    if(!transcriptID.equals(transcriptIDnew) && !isFirstTranscript) {
                        // store old transcript
                        currentGene.transcripts.put(transcriptID, currentTranscript);

                        // handle new transcript
                        currentTranscript = new RegionVector();
                    }
                    isFirstTranscript = false;
                    transcriptID = transcriptIDnew;
                    firstExons = false;
                    eStart = Integer.parseInt(splitLineGtf[3]);
                    eEnd = Integer.parseInt(splitLineGtf[4]);
                    Region reg = new Region(eStart, eEnd, geneID);
                    currentTranscript.regions.add(reg);
                }
            } else {
                if(!geneID.equals("blub")) {
                    // already seen all entries of chromosome chr
                    break;
                }
            }
        }
        if(!geneHandled) {
            // maybe a gene got left unhandled
            if(currentGene != null) {
                if(!firstExons) {
                    // store old transcript
                    currentGene.transcripts.put(transcriptID, currentTranscript);
                }
                // handle old gene
                genes.remove(geneID);
                genes.put(geneID, currentGene);
            }
        }
        currentGene = genes.get("YJR159W");
        return genes;
    }

    /**
     * Decides whether the introns of one read are inconsistent with the other.
     * @param introns1
     * @param regions2
     * @return True if the readpair is inconsistent for introns of readx, exons of read(x^(-1)).
     */
    private static boolean inconsistency(RegionVector introns1, RegionVector regions2) {
        boolean test1 = false;  // there exists a region with i.x1 == r.x2
        boolean test2 = false;  // there exists a region with i.x2 == r.x1
        boolean test3 = false;  // r2 follows directly after r1
        int counter = 0;
        int counterForR1 = -4;
        int counterForR2 = -2;
        RegionVector editedRegions;
        for(Region intron : introns1.regions) {
            test1 = false;
            test2 = false;
            test3 = false;
            editedRegions = new RegionVector();
            for(Region region : regions2.regions) {
                if((region.x2 >= intron.x1 && region.x2 <= intron.x2) ||       // ende ragt rein
                        (region.x1 >= intron.x1 && region.x1 <= intron.x2) ||   // anfang ragt rein
                        (region.x1 <= intron.x1 && region.x2 >= intron.x2)) {   // geht komplett drüber
                    // is overlapping
                    editedRegions.regions.add(region);
                }
            }
            if(editedRegions.regions.isEmpty()) {
                test1 = true;
                test2 = true;
                test3 = true;
            } else {
                for(Region region : regions2.regions) {
                    if(region.x2 == intron.x1 - 1) {
                        counterForR1 = counter;
                        test1 = true;
                    }
                    if(region.x1 == intron.x2 + 1) {
                        test2 = true;
                        counterForR2 = counter;
                    }
                    if(counterForR2 == counterForR1 + 1) {
                        test3 = true;
                    }
                    counter++;
                }
            }
            if(!test1 && !test2 && !test3) {
                return true;
            }
            counter = 0;
        }
        return false;
    }

    /**
     * Decides whether the introns of one read are inconsistent with the other.
     * @param introns1
     * @param regions2
     * @return True if the readpair is inconsistent for introns of readx, exons of read(x^(-1)).
     */
    private static boolean inconsistency2(RegionVector introns1, RegionVector regions2) {
        for(Region intron : introns1.regions) {
            boolean intronNotInconsistent = true;
            // test if overlapping
            RegionVector overlapping = new RegionVector();
            for(Region region : regions2.regions) {
                if((region.x2 >= intron.x1 && region.x2 <= intron.x2) ||       // ende ragt rein
                        (region.x1 >= intron.x1 && region.x1 <= intron.x2) ||   // anfang ragt rein
                        (region.x1 <= intron.x1 && region.x2 >= intron.x2)) {   // geht komplett drüber
                    // is overlapping
                    overlapping.regions.add(region);
                    intronNotInconsistent = false;
                }
            }
            for(Region region : overlapping.regions) {
                // test if right one
                if(region.x2 == intron.x1 + 1) {
                    intronNotInconsistent = true;
                }
                if(region.x1 == intron.x2 - 1) {
                    intronNotInconsistent = true;
                }
            }
            if(overlapping.regions.size() == 2) {
                // r2 follows immediately after r1
                intronNotInconsistent = true;
            }
            if(!intronNotInconsistent) {
                return true;
            }
        }
        return false;
    }

    private static boolean inconsistency3(RegionVector regions1, RegionVector regions2) {
        // init
        int firstStart1 = regions1.regions.first().x1;
        int lastEnd1 = regions1.regions.last().x2;
        int firstStart2 = regions2.regions.first().x1;
        int lastEnd2 = regions2.regions.last().x2;
        // overlap?
        if(lastEnd1 < firstStart2 || lastEnd2 < firstStart1) {
            // no overlap
            return false;
        }
        // cut
        if(firstStart1 <= firstStart2) {
            if(lastEnd1 >= lastEnd2) {
                // case 1: cutting once
                regions1 = regions1.cut(firstStart2, lastEnd2);
            } else {
                // case 3: cutting twice
                regions1 = regions1.cut(firstStart2, lastEnd1);
                regions2 = regions2.cut(firstStart2, lastEnd1);
            }
        } else {
            if(lastEnd1 <= lastEnd2) {
                // case 2: cutting once
                regions2 = regions2.cut(firstStart1, lastEnd1);
            } else {
                // case 4: cutting twice
                regions1 = regions1.cut(firstStart1, lastEnd2);
                regions2 = regions2.cut(firstStart1, lastEnd2);
            }
        }
        boolean equal;
        if(regions1 == null || regions2 == null) {
            equal = false;
        } else {
            equal = regions1.isEqual(regions2);
        }
        return !equal;
    }


        /**
         * Transforms the AlignmentBlocks of a read into a region vector.
         * @param read
         * @param rv
         * @return RegionVector
         */
    private static RegionVector blocksToRV(Read read, RegionVector rv) {
        for(AlignmentBlock ab : read.blocks) {
            int refStart = ab.getReferenceStart();      // start on ref
            // int readStart = ab.getReadStart();          // start on read
            int refEnd = refStart + ab.getLength() - 1;     // end on ref
            // int readEnd = readStart + ab.getLength() - 1;   // end on read
            rv.regions.add(new Region(refStart, refEnd, ""));
        }
        read.regions = rv;
        return rv;
    }

    /**
     * Calculates the distance to the nearest gene for intergenic readpairs.
     * @param start
     * @param end
     * @param fwFlag
     * @return the minimal distance to the nearest gene
     */
    private static int distanceToGenes(int start, int end, boolean fwFlag) {
        HashSet<Region> neighborsRight = new HashSet<>();
        HashSet<Region> neighborsLeft = new HashSet<>();
        // collect the nearest neighbors
        if(fwFlag) {
            neighborsRight = genesFw.getIntervalsRightNeighbor(start, end, neighborsRight);
            neighborsLeft = genesFw.getIntervalsLeftNeighbor(start, end, neighborsLeft);
        } else {
            neighborsRight = genesRw.getIntervalsRightNeighbor(start, end, neighborsRight);
            neighborsLeft = genesRw.getIntervalsLeftNeighbor(start, end, neighborsLeft);
        }
        // find the closest start/stop indices
        Integer distLeft = null;
        for(Region rightNeighbor : neighborsRight) {
            if(distLeft == null) {
                distLeft = Math.min(Math.max(0, rightNeighbor.x1 - end), Math.max(0, rightNeighbor.x2 - end));
            } else {
                distLeft = Math.min(distLeft,
                        Math.min(Math.max(0, rightNeighbor.x1 - end), Math.max(0, rightNeighbor.x2 - end)));
            }
        }
        Integer distRight = null;
        for(Region leftNeighbor : neighborsLeft) {
            if(distRight == null) {
                distRight = Math.min(Math.max(0, start - leftNeighbor.x2), Math.max(0, start - leftNeighbor.x1));
            } else {
                distRight = Math.min(distRight,
                        Math.min(Math.max(0, start - leftNeighbor.x2), Math.max(0, start - leftNeighbor.x1)));
            }
        }
        // calculate distance
        if(distRight != null) {
            if(distLeft == null) {
                // no left neighbors
                return distRight;
            } else {
                // both neighbors
                return Math.min(distLeft, distRight);
            }
        } else {
            if(distLeft != null) {
                // no right neighbors
                return distLeft;
            } else {
                // no neighbors ??
                return -1;
            }
        }
    }

    /**
     *
     * Decides if a given region, depending on the given strandness, contains a gene.
     * @param alStart
     * @param alEnd
     * @param leadingStrand
     * @return True if the given region contains a gene.
     */
    private static boolean doesContainFlag(int alStart, int alEnd, Boolean leadingStrand) {
        boolean doesContain = false;
        if(frstrand == null) {  // ???
            if(checkContainsGenes(alStart, alEnd, true) ||      //TODO: maybe && ?
                    checkContainsGenes(alStart, alEnd, false)) {
                doesContain = true;
            }
        } else {
            if(leadingStrand) {
                if(checkContainsGenes(alStart, alEnd, true)) {
                    doesContain = true;
                }
            } else {
                if(checkContainsGenes(alStart, alEnd, false)) {
                    doesContain = true;
                }
            }
        }
        return doesContain;
    }

    /**
     * Calculates part 1 of the pre check for (possibly?) intergenic readpairs.
     * @param start start of read
     * @param end end of read
     * @param fwFlag searching for which strand
     * @return true if the check alerts, full genes are contained and the readpair could be discarded
     */
    private static boolean checkContainsGenes(int start, int end, boolean fwFlag) {
        HashSet<Region> containedGenes = new HashSet<>();
        // special region contains genes?
        if(fwFlag) {
            containedGenes = genesFw.getIntervalsSpannedBy(start, end, containedGenes);
        } else {
            containedGenes = genesRw.getIntervalsSpannedBy(start, end, containedGenes);
        }
        if(!containedGenes.isEmpty()) {
            // read could be (PRE)-intergenic -> maybe skip
            return true;
        } else {
            return false;
        }
    }

    /**
     * Decides if a given region, depending on the given strandness, is spanned by a gene.
     * @param alStart
     * @param alEnd
     * @param fopStrand
     * @return True if the given region is NOT spanned by a gene.
     */
    private static boolean isNotSpannedFlag(int alStart, int alEnd, Boolean fopStrand) {
        boolean noSpanning = false;
        if(frstrand == null) {  // ???
            if(checkWithinGenes(alStart, alEnd, true) && checkWithinGenes(alStart, alEnd, false)) {
                noSpanning = true;
            }
        } else {
            if(fopStrand) {
                if(checkWithinGenes(alStart, alEnd, true)) {
                    noSpanning = true;
                }
            } else {
                if(checkWithinGenes(alStart, alEnd, false)) {
                    noSpanning = true;
                }
            }
        }
        return noSpanning;
    }

    /**
     * Calculates part 2 of the pre check for (possibly?) intergenic readpairs.
     * @param start start of read
     * @param end end of read
     * @param fwFlag searching for which strand
     * @return true if the check alerts, no genes are spanning and the readpair could be discarded
     */
    private static boolean checkWithinGenes(int start, int end, boolean fwFlag) {
        HashSet<Region> spanningGenes = new HashSet<>();
        // special region within genes?
        if(fwFlag) {
            spanningGenes = genesFw.getIntervalsSpanning(start, end, spanningGenes);
        } else {
            spanningGenes = genesRw.getIntervalsSpanning(start, end, spanningGenes);
        }
        if(spanningGenes.isEmpty()) {
            // read could be (PRE)-intergenic -> maybe skip
            return true;
        } else {
            return false;
        }
    }

    /**
     * Loads the genes in the gtf annotation for specific chromosome into an interval tree.
     * @return an interval tree of the genes of the specified chromosome
     * @throws IOException
     */
    private static IntervalTree loadGtfPerChromosome(String chr, boolean strand) throws IOException {
        IntervalTree<Region> tree = new IntervalTree<>();
        char s;
        if(strand) {
            s = '+';
        } else {
            s = '-';
        }
        int gStart = -1;
        int gEnd = -1;
        String lineGtf;
        String[] splitLineGtf;
        String[] attributes;
        String geneID;

        BufferedReader brGtf = new BufferedReader(new FileReader(gtfFile));
        while((lineGtf = brGtf.readLine()) != null) {
            if(lineGtf.startsWith("#")) continue;

            splitLineGtf = lineGtf.split("\t");
            if(splitLineGtf[0].equals(chr)) {
                if(splitLineGtf[2].equals("gene") && splitLineGtf[6].charAt(0) == s) {
                    attributes = cutEmptySpace(splitLineGtf[8].split(";"));
                    geneID = findArgument("gene_id", attributes);
                    gStart = Integer.parseInt(splitLineGtf[3]);
                    gEnd = Integer.parseInt(splitLineGtf[4]);
                    Region reg = new Region(gStart, gEnd, geneID);
                    tree.add(reg);
                }
            } else {
                if(!tree.isEmpty()) {
                    // already seen all entries of chromosome chr
                    break;
                }
            }
        }
        return tree;
    }

    /**
     * Cuts annoying empty spaces at the start of the semicolon-separated attributes.
     * @param words attributes/last field
     * @return the same field but without empty spaces at the beginning of certain attributes fields
     */
    private static String[] cutEmptySpace(String[] words) {
        for(int i = 0; i < words.length; i++) {
            if(words[i].startsWith(" ")) {
                words[i] = words[i].substring(1);
            }
        }
        return words;
    }

    /**
     * @param search the concerned attribute
     * @param attributes last field
     * @return the attribute value or the empty string, if nothing could be found
     */
    private static String findArgument(String search, String[] attributes){
        // finds argument value in an array of * argument_name "value" * blocks
        for(int i = 0; i < attributes.length; i++) {
            if(attributes[i].startsWith(search)) {
                if(search.equals("exon_number")) {
                    return attributes[i].substring(attributes[i].length());
                } else {
                    return attributes[i].split("\"")[1];
                }
            }
        }
        return "";
    }

    /**
     * Prints usage when program gets invoked with false parameters.
     *
     */
    private static void printUsage() {
        // parameters
        System.out.println("Following parameters (in any order) are necessary:\n" +
                "-gtf (file name): gtf annotated file\n" +
                "-bam (file name): binary sam formatted file\n" +
                "-o (file name): output tsv\n" +
                "*OPTIONAL* -frstrand (boolean): strandness of the experiment (true ~ first read maps sense to" +
                "transcribed region");
        // what will be done
        System.out.println("If invoked correctly, the program will annotate information on all mapped reads pairs. " +
                "The output contains one line per read pair in the format readid<tab>annotation.");
    }

    public static void main(String[] args) {
        nonDuplicates = false;
        // argparse
        try {
            if(args.length < 6) {
                printUsage();
            } else if(args.length == 6 || args.length == 8){
                for(int i = 0; i < args.length - 1; i += 2) {
                    if(args[i].equals("-gtf")) {
                        gtfFile = new File(args[i + 1]);
                    } else if(args[i].equals("-bam")) {
                        bamFile = new File(args[i + 1]);
                    } else if(args[i].equals("-o")) {
                        outputTsvString = args[i + 1];
                    } else if(args[i].equals("-frstrand")) {
                        frstrand = Boolean.parseBoolean(args[i + 1]);
                    } else {
                        printUsage();
                    }
                }
                iterateBamFile();
            } else {
                printUsage();
            }
        } catch (IOException e) {
            System.err.println("Error: Something went wrong.");
        }
    }
}
