package org.eloe.dataout;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;
import org.eloe.parameters.PGeneToString;

//gene_na;_genes_and_flanks.txt
public class LineMaker2 implements ILineMaker {
	
	private Organism organism;
	private PGeneToString parameters;
	
	public LineMaker2(Organism organism, PGeneToString parameters) {
		super();
		this.organism = organism;
		this.parameters = parameters;
	}

	@Override
	public String makeLine(Gene gene) {
		return new String(gene.toString(organism, parameters)
				+gene.getEeiMass()[organism.getIndex()]+"\n"
				+String.valueOf(gene.getUserNA(parameters)));
	}

}
