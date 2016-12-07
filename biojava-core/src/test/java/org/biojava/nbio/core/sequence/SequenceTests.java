package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.FastaGeneWriter;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderFormat;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * Created by ashwinsl on 12/5/16.
 */
public class SequenceTests {

    @Test
    public void testChromosomeSequences(){

        try {
            ArrayList<GeneSequence> sequences = new ArrayList<GeneSequence>();

            ChromosomeSequence seq1 = new ChromosomeSequence("ATATATATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGC");
            GeneSequence gene1 = seq1.addGene(new AccessionID("gene1"), 1, 20, Strand.POSITIVE);

            gene1.addExon(new AccessionID("t1_1_10"), 1, 10);
            gene1.addExon(new AccessionID("t1_14_18"), 14, 18);
            GeneSequence gene2 = seq1.addGene(new AccessionID("gene2"), 1, 20, Strand.NEGATIVE);

            gene2.addExon(new AccessionID("t2_1_10"), 1, 10);
            gene2.addExon(new AccessionID("t2_14_18"), 14, 18);
            sequences.add(gene1);
            sequences.add(gene2);

            FastaGeneWriter fastaWriter = new FastaGeneWriter(System.out, sequences, new GenericFastaHeaderFormat<GeneSequence, NucleotideCompound>(), true);
            fastaWriter.process();

            assertEquals(seq1, gene1.getParentChromosomeSequence());
            assertEquals(seq1, gene2.getParentChromosomeSequence());

            assertEquals(20, gene1.getLength());
            assertEquals(20, gene2.getLength());

            assertEquals(gene2.getSequence5PrimeTo3Prime().toString(), "TATATATATATATATATATA");

            TranscriptSequence ts1 = gene1.addTranscript(new AccessionID("tt1_1_40"), 1, 40);
            assertEquals(ts1.getLength(), 40);
            ts1.addCDS(new AccessionID("tt1_1_40-cds1"), 4, 18, 0);
            AccessionID testAccession = new AccessionID("tt1_1_40-cds2");

            ts1.addCDS(testAccession, 20, 30, 1);
            ts1.addStartCodonSequence(new AccessionID("ss1_1_3"), 1, 3);
            ts1.addStopCodonSequence(new AccessionID("ss1_19_21"), 19, 21);

            List<ProteinSequence> proteinSequences1 = ts1.getProteinCDSSequences();
            for(ProteinSequence ps: proteinSequences1){
                if(ps.toString().equals("IYI") || ps.toString().equals("YIYIY")){
                    continue;
                }else {
                    Assert.fail("Protein Sequences are wrong");
                }
            }
            assertEquals(ts1.getDNACodingSequence().toString(), "TATATATATATATATTATATATATAT");

            assertEquals(ts1.getProteinSequence().toString(), "YIYIYYIY");

            TranscriptSequence ts2 = gene2.addTranscript(new AccessionID("tt2_1_40"), 1, 40);
            assertEquals(ts2.getLength(), 40);
            ts2.addCDS(new AccessionID("tt2_1_40-cds1"), 4, 18, 0);
            ts2.addCDS(new AccessionID("tt2_1_40-cds2"), 20, 30, 1);

            List<ProteinSequence> proteinSequences2 = ts2.getProteinCDSSequences();
            for(ProteinSequence ps: proteinSequences2){
                if(ps.toString().equals("IYI") || ps.toString().equals("YIYIY")){
                    continue;
                }else {
                    Assert.fail("Protein Sequences are wrong");
                }
            }

            assertEquals(ts2.getDNACodingSequence().toString(), "TATATATATATTATATATATATATAT");

            assertEquals(ts1.removeCDS("tt1_1_40-cds2").getAccession(), testAccession);

            gene1.addIntronsUsingExons();
            gene2.addIntron(new AccessionID("t2_11_13"), 11, 13);


            List<IntronSequence> intronSequences1 = gene1.getIntronSequences();
            List<IntronSequence> intronSequences2 = gene2.getIntronSequences();

            assertEquals(intronSequences1.size(), 1);
            assertEquals(intronSequences2.size(), 1);

            for(IntronSequence intronSequence: intronSequences1){
                assertEquals(intronSequence.getLength(), 3);
            }
            for(IntronSequence intronSequence: intronSequences2){
                assertEquals(intronSequence.getLength(), 3);
            }

            IntronSequence intronSequence = gene2.removeIntron("t2_11_13");
            assertEquals(intronSequence.getAccession().getID(), "t2_11_13");

            ExonSequence exonSequence = gene2.removeExon("t2_14_18");
            assertEquals(exonSequence.getAccession().getID(), "t2_14_18");

            assertEquals(gene2.getIntronSequences().size(), 0);



        }catch (Exception e){
            System.out.print(e);
            Assert.fail("Exception when trying to initialize chromose sequence and gene sequence");
        }

    }
}
