package org.eloe.dataout;

import org.eloe.organism.Gene;

public class LineMaker3 implements ILineMaker {
	
	private int ind;
	private int i;

	@Override
	public String makeLine(Gene gene) {
		StringBuilder s = new StringBuilder();
		s.append((i+1)).append("\t").append(gene.getNumber())
						.append("\t").append(gene.getEeiMass()[ind])
						.append("\t").append(gene.getLocusTag())
						.append("\t").append(gene.getName())
						.append("\t").append(gene.getProteinName())
						.append("\t").append(gene.getNote());
		if (gene.isRibosomal()) 
			s.append("\trib");
		s.append("\n");
		return s.toString();
	}

}
