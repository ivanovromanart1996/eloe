package org.eloe.dataout;

import org.eloe.organism.Gene;

//num_cis_eei; _number_eei.txt
public class LineMaker1 implements ILineMaker {
	
	private int index;
	private boolean isEukaryota;

	public LineMaker1(int index, boolean isEukaryota) {
		this.index = index;
		this.isEukaryota = isEukaryota;
	}

	@Override
	public String makeLine(Gene gene) {
		StringBuilder s = new StringBuilder();
		s.append(gene.getNumber());
		if (isEukaryota) {
			if (gene.isFirstCistron() && gene.isLastCistron()) {
				s.append("\t1\t0\t0\t0");
			}
			else {
				if (gene.isMiddleCistron()) {
					s.append("\t0\t0\t1\t0");
				}
				else {
					if (gene.isFirstCistron() && !gene.isLastCistron()) {
						s.append("\t0\t1\t0\t0");
					}
					else {
						s.append("\t0\t0\t0\t1");
					}
				}
			}
		}
		s.append("\t");
		s.append(gene.getEeiMass()[index]);
		return s.toString();
	}

}
