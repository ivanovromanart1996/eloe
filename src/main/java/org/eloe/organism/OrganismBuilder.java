package org.eloe.organism;

import java.io.File;

import org.apache.log4j.Logger;

import org.eloe.extractor.Extractor;
import org.eloe.extractor.ExtractorFactory;
import org.eloe.parameters.PExtraction;

public abstract class OrganismBuilder {
	
	private static final Logger LOG = Logger.getLogger(OrganismBuilder.class);
	
	//Get organism data from directory
	public static Organism getOrganism(File orgFolder, PExtraction parameters, File resultFolder) {
		Organism organism = new Organism();
		File[] orgFileList = orgFolder.listFiles();
		
		for (File orgFile : orgFileList) {
			Extractor extractor = ExtractorFactory.getExtractor(orgFile, parameters, resultFolder);
			if (extractor == null || !orgFile.isFile()) {
				LOG.info("Skipping unsupported input: " + orgFile);
				continue;
			}
			Organism chr = extractor.extract();
			if (chr == null) {
				continue;
			}
			organism.addChromosome(chr);
		}
		
		if (organism.getGenes().size() == 0) {
			LOG.info("There are not genes in files!");
			return null;
		}
		
		//Output of extracted data to the console
		LOG.info(organism);
		
		//Annotate operons if the organism is prokaryotic.
    	organism.markCistrons();
    	
		return organism;
	}
}
