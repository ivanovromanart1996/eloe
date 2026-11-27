package org.eloe.engine.te.lci.lci1;

import org.apache.log4j.Logger;

import org.eloe.engine.te.lci.ILciCalculator;
import org.eloe.organism.IGeneWithFlanks;
import org.eloe.parameters.PLCI1;

public abstract class Lci1CalculatorAbstract implements ILciCalculator {
	
	private static final Logger LOG = Logger.getLogger(Lci1CalculatorAbstract.class);
	
	private PLCI1 parameters;
	private int distance;
	private int rightFlank;
	
	public Lci1CalculatorAbstract(PLCI1 parameters) {
		this.parameters = parameters;
		distance = parameters.getLmax1() - parameters.getLmin1();
		rightFlank = parameters.getSmax1() + parameters.getLmax1();
	}

	@Override
	public double calculate(IGeneWithFlanks gene) {
		
		if (rightFlank > gene.getRightFlankLength()) {
			LOG.info("The right flank length of gene "+gene.getName()+" is not enough for LCI1 calculation.");
			return 0;
		}
		
		//The length of the part of the gene under study (3 - stop codon)
		int lenght = gene.getLength()-3+rightFlank;
		//The boundary of the left edge of the sliding frame
		int edge = lenght-2*parameters.getSmax1()-parameters.getLmax1()+1;
		if (edge < 0) {
			LOG.info("The length of gene "+gene.getName()+" with flanks is not enough for LCI1 calculation.");
			return 0;
		}
		
		char[] na = new char[lenght];
		//Get the gene without flanks + h nucleotides at the 3' end
		System.arraycopy(gene.getNaWithFlanks(), gene.getLeftFlankLength(), na, 0, lenght);
		
		//Repeat counter
		int counter = 0;
		for (int j = 0; j < edge; j++) {
			for (int s = parameters.getSmin1(); s <= parameters.getSmax1(); s++) {
				
				//Buffer for the repeat
				char[] char1 = new char[s];
				//Copying the sequence
				System.arraycopy(na, j, char1, 0, s);
				//Buffer for the complementary sequence
				String con1 = String.valueOf(reverseComplement(char1));
				
				//Length of search area
				int lenght1 = s+distance;
				//Search site
				char[] char2 = new char[lenght1];
				//Extract the search region.
				System.arraycopy(na, j+s+parameters.getLmin1(), char2, 0, lenght1);
				//Search site
				String con2 = String.valueOf(char2);
				
				int poz = 0;
				//Search for repeats (complementary regions).
				int k = con2.indexOf(con1, poz);
				while (k != -1) {
					//Counting repeats
					counter++;
					poz = k + 1;
					k = con2.indexOf(con1, poz);
				}
			}
		}
		
		return (double) counter/(gene.getLength()-3);
	}

	abstract protected char[] reverseComplement(char[] char1);
}
