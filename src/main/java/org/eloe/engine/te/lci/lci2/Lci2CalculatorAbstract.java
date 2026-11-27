package org.eloe.engine.te.lci.lci2;

import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import org.eloe.engine.te.lci.ILciCalculator;
import org.eloe.organism.IGeneWithFlanks;
import org.eloe.parameters.PLCI2;

public abstract class Lci2CalculatorAbstract implements ILciCalculator {
	
private static final Logger LOG = Logger.getLogger(Lci2CalculatorAbstract.class);
	
	private PLCI2 parameters;
	private int distance;
	private int rightFlank;
	//Map for recording the energy of dinucleotides
	private Map<String, Double> map = new HashMap<String, Double>();
	//Loop energy array
	private double loopenergy[];
	
	public Lci2CalculatorAbstract(PLCI2 parameters, Map<String, Double> map, double[] loopenergy) {
		this.parameters = parameters;
		distance = parameters.getLmax2()-parameters.getLmin2();
		rightFlank = parameters.getSmax2()+parameters.getLmax2();
		this.map = map;
		this.loopenergy = loopenergy;
	}

	@Override
	public double calculate(IGeneWithFlanks gene) {
		
		if (rightFlank > gene.getRightFlankLength()) {
			LOG.info("The right flank length of gene "+gene.getName()+" is not enough for LCI2 calculation.");
			return 0;
		}
		
		//The length of the part of the gene under study (3 - stop codon)
		int lenght = gene.getLength()-3+rightFlank;
		//The boundary of the left edge of the sliding frame
		int edge = lenght-2*parameters.getSmax2()-parameters.getLmax2()+1;
		if (edge < 0) {
			LOG.info("The length of gene "+gene.getName()+" with flanks is not enough for LCI2 calculation.");
			return 0;
		}
		
		char[] na = new char[lenght];
		//Get the gene without flanks + h nucleotides at the 3' end
		System.arraycopy(gene.getNaWithFlanks(), gene.getLeftFlankLength(), na, 0, lenght);
		
		//Total energy of interaction
		double energyFull = 0;
		for (int j = 0; j < edge; j++) {
			//Hairpin energy
			double energy = 0;
			for (int s = parameters.getSmin2(); s <= parameters.getSmax2(); s++) {
				
				//Buffer for the repeat
				char[] char1 = new char[s];
				//Copying the sequence
				System.arraycopy(na, j, char1, 0, s);
				
				//Energy of the selected repeat
				double energy_char1 = 0;
				
				for (int m = 0; m < s-1; m++) {
					String dinukleotid = char1[m]+""+char1[m+1];
					
					if (map.containsKey(dinukleotid)) {
						//Repeat interaction energy
						energy_char1 += map.get(dinukleotid);
					}
				}
				
				//Buffer for the complementary sequence
				String con1 = String.valueOf(reverseComplement(char1));
				
				//Length of search area
				int lenght1 = s+distance;
				//Search site
				char[] char2 = new char[lenght1];
				//Extract the search region.
				System.arraycopy(na, j+s+parameters.getLmin2(), char2, 0, lenght1);
				//Search site
				String con2 = String.valueOf(char2);
				
				int poz = 0;
				//Loop energy
				double loop1 = 0;
				
				//Search for repeats (complementary regions).
				int k = con2.indexOf(con1, poz);
				while (k != -1) {
					int loopsize = parameters.getLmin2()+k-1;
					
					//Loop energy
					loop1 = loopenergy[(loopsize < loopenergy.length) ? loopsize : (loopenergy.length-1)];
					double e = energy_char1 - loop1;
					
					//If the energy is greater than the threshold of the existence of the loop, then take into account this loop.
					if (e >= parameters.getMinloopEnergy()) {
						energy += e;
					}
					poz = k + 1;
					k = con2.indexOf(con1, poz);
				}
			}
			//Total energy
			energyFull += energy;
		}
		
		return energyFull/(gene.getLength()-3);
	}
	
	abstract protected char[] reverseComplement(char[] char1);
}
