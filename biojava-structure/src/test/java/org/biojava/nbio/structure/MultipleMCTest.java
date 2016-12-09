package org.biojava.nbio.structure;

import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcMain;
import org.biojava.nbio.structure.align.multiple.mc.MultipleMcParameters;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by ashwinsl on 12/9/16.
 */
public class MultipleMCTest {

    @Test
    public void testMultipleMC(){
        try {
            List<String> names = Arrays.asList("d1hm9a1", "d1kk6a_", "d1krra_", "d1lxaa_", "d1ocxa_", "d1qrea_", "d1xata_", "d3tdta_");

            AtomCache cache = new AtomCache();
            List<Atom[]> atomArrays = new ArrayList<Atom[]>();

            List<StructureIdentifier> ids = new ArrayList<StructureIdentifier>();
            for (String name:names)	{
                StructureIdentifier id = new StructureName(name);
                ids.add(id);
                atomArrays.add(cache.getAtoms(id));
            }

            //Here the multiple structural alignment algorithm comes in place to generate the alignment object
            MultipleMcMain algorithm = new MultipleMcMain(new CeMain());
            MultipleMcParameters params = (MultipleMcParameters) algorithm.getParameters();
            params.setMinBlockLen(15);
            params.setMinAlignedStructures(10);

            MultipleAlignment result = algorithm.align(atomArrays);
            result.getEnsemble().setStructureIdentifiers(ids);

            //Information about the alignment
            result.getEnsemble().setAlgorithmName(algorithm.getAlgorithmName());
            result.getEnsemble().setVersion(algorithm.getVersion());

            //Output the sequence alignment + transformations
            System.out.println(MultipleAlignmentWriter.toFatCat(result));
            //System.out.println(MultipleAlignmentWriter.toFASTA(result));
            System.out.println(MultipleAlignmentWriter.toTransformMatrices(result));
            System.out.println(MultipleAlignmentWriter.toXML(result.getEnsemble()));
        } catch (Exception e) {
            e.printStackTrace();
            Assert.fail("Exception when trying to run MultipleMC");
        }
    }
}
