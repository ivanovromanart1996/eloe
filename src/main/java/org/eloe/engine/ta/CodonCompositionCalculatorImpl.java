package org.eloe.engine.ta;

import java.util.Map;

import org.apache.log4j.Logger;

import org.eloe.organism.IGene;

public class CodonCompositionCalculatorImpl implements ICodonCompositionCalculator {
	
	private static final Logger LOG = Logger.getLogger(CodonCompositionCalculatorImpl.class);
	
	private Map<String, Integer> mapOfCodons;

	public CodonCompositionCalculatorImpl(Map<String, Integer> mapOfCodons) {
		this.mapOfCodons = mapOfCodons;
	}
	
	@Override
	public int[] calculate(IGene gene, int unUsedCodonsL, int unUsedCodonsR) {
		int[] codonsMatrix = new int[mapOfCodons.size()];
		char[] codon = new char[3];
		int codonsNumber = gene.getCodonsNumber();
		if (codonsNumber <= unUsedCodonsL+unUsedCodonsR) {
			LOG.info("Not enough codons in gene "+gene.getName()+" to calculate CodonComposition.");
			return codonsMatrix;
		}
		for (int i = unUsedCodonsL; i < codonsNumber-unUsedCodonsR; i++) {
			System.arraycopy(gene.getNA(), 3*i, codon, 0, 3);
			String codonString = String.valueOf(codon);
			if (mapOfCodons.containsKey(codonString)) {
				codonsMatrix[mapOfCodons.get(codonString)]++;
			}
			else {
				LOG.info("Unexpected codon "+codonString+" was found in gene "+gene.getName()+".");
			}
		}
		return codonsMatrix;
	}

}
