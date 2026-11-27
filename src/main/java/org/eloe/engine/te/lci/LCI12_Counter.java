package org.eloe.engine.te.lci;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import org.eloe.engine.ta.FreqCounter;
import org.eloe.organism.Organism;
import org.eloe.parameters.PLCI2;

// Program for calculating the local complementarity index

public class LCI12_Counter {
	
	private static final Logger LOG = Logger.getLogger(LCI12_Counter.class);

	private static final String ENERGY = "lci2double.txt";
	private PLCI2 parameters;
	//private int forSigma1; //Sigma multiplier
	//private int forSigma2; //Sigma multiplier
	private double m1 = 0; //Average LCI value
	private double sigma1 = 0; //Average deviation of LCI
	private double m2 = 0; //Average LCI value
	private double sigma2 = 0; //Average deviation of LCI
	private double sred1 = 0;
	private double disp1 = 0;
	private double sred2 = 0;
	private double disp2 = 0;
	private double loopenergy[] = null;
	private Map<String, Double> map = new HashMap<String, Double>(); //Map for recording the energy of dinucleotides
	
	//Constructor
	public LCI12_Counter(PLCI2 parameters) throws IOException, NumberFormatException {
		this.parameters = parameters;
		//Auxiliary variable
		//this.forSigma1 = parameters.forSigma1;
		//this.forSigma2 = parameters.forSigma2;
		String loopsenergy = "";
		
		BufferedReader reader = null;
		String line = new String();
		try {
			//Reading a file with energies of dinucleotides and loops
			reader = new BufferedReader(new InputStreamReader(FreqCounter.class.getClassLoader().getResourceAsStream(ENERGY)));
			while (!(line = reader.readLine()).startsWith("//")){
				String arr[] = line.split(" - ");
				double a = Double.parseDouble(arr[1]);
				map.put(arr[0], a);
			}
			while (!((line = reader.readLine()).startsWith("////"))) {
				loopsenergy += line;
			}
		} catch (NumberFormatException e) {
			e.printStackTrace();
			LOG.error("Sorry, something wrong!", e);
			throw e;
		} catch (IOException e) {
			e.printStackTrace();
			LOG.error("Sorry, something wrong!", e);
			throw e;
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		String loopsenergy1[] = loopsenergy.split(";");
		
		//Number of loop lengths
		int massivesize = loopsenergy1.length;
		
		//Array for loop energies
		loopenergy = new double[massivesize];
		try {
			for (int i = 0; i < massivesize; i++) {
				loopenergy[i] = Double.parseDouble(loopsenergy1[i].trim());
			}
		} catch (NumberFormatException e) {
			e.printStackTrace();
			LOG.error("Sorry, something wrong!", e);
			throw e;
		}
	}
	
	//TODO Calculation of LCI12
	public void countLCI(Organism organism) {
		//LCI1
		double lci1 = 0;
		//LCI2
		double lci2 = 0;
		
		for (int i = 0; i < organism.getGenes().size(); i++) {
			//Extension beyond the right end of the gene
			int h = parameters.getSmax2()+parameters.getLmax2();
			if (h > organism.getGenes().get(i).getRightFlankLength()) {
				//h = rightFl;
				continue;
			}
			
			//The length of the part of the gene under study (3 - stop codon)
			int lenght = organism.getGenes().get(i).getLength()-3+h;
			
			//A string to record the gene being studied
			char gene[] = new char[lenght];
			
			//Length of left flank
			int leftFl = organism.getGenes().get(i).getLeftFlankLength();
			
			//Get the gene without flanks + h nucleotides at the 3' end
			System.arraycopy(organism.getGenes().get(i).getNaWithFlanks(), leftFl, gene, 0, lenght);
			
			//The boundary of the left edge of the sliding frame
			int edge = lenght-2*parameters.getSmax2()-parameters.getLmax2()+1;
			if (edge < 0) continue;
			
			//Total energy of interaction
			double energyFull = 0;
			
			//Repeat counter
			int counter = 0;
			
			for (int j = 0; j < edge; j++){
				//Hairpin energy
				double energy = 0;
				
				for (int s = parameters.getSmin2(); s <= parameters.getSmax2(); s++){
					//Buffer for the repeat
					char con1[] = new char[s];
					//Copying the sequence
					System.arraycopy(gene, j, con1, 0, s);
					//Energy of the selected repeat
					double energy_con1 = 0;
					
					for (int m = 0; m < s-1; m++) {
						String dinukleotid = con1[m]+""+con1[m+1];
						if (map.containsKey(dinukleotid)) {
							//Repeat interaction energy
							energy_con1 += map.get(dinukleotid);
						}
					}
					
					//Buffer for the complementary sequence
					StringBuilder con2 = new StringBuilder("");
					
					//Complementary Sequence Building
					for (int k = s-1; k > -1 ; k--){
						if (con1[k] == 'c') {
							con2.append("g");
						}
						else {
							if (con1[k] == 'a') {
								con2.append("t");
							}
							else {
								if (con1[k] == 'g') {
									con2.append("c");
								}
								else con2.append("a");
							}
						}
					}
					String con22 = con2.toString();
					
					//Length of search area
					int lenght1 = s+this.parameters.getLmax2()-this.parameters.getLmin2();
					//Array for the search site
					char con3[] = new char[lenght1];
					//Extract the search region.
					System.arraycopy(gene, j+s+parameters.getLmin2(), con3, 0, lenght1);
					//Search site
					String gene1 = String.valueOf(con3);
					
					int poz = 0;
					//Loop energy
					double loop1 = 0;
					
					int k = gene1.indexOf(con22, poz);
					while (k != -1) {
						//Counting repeats
						counter++;
						
						int loopsize = parameters.getLmin2()+k-1;
						
						//Loop energy
						loop1 = loopenergy[(loopsize < loopenergy.length) ? loopsize : (loopenergy.length-1)];
						double e = energy_con1 - loop1;
						//If the energy is greater than the threshold of the existence of the loop, then take into account this loop.
						if (e >= parameters.getMinloopEnergy()) {
							energy += e;
						}
						poz = k+1;
						k = gene1.indexOf(con22, poz);
					}
				}
				
				//Total energy
				energyFull += energy;
			}
			
			//Number of repeats per nucleotide (subtract 3 for the stop codon)
			lci1 = (double)counter/(organism.getGenes().get(i).getLength()-3);
			//Store LCI1 in the gene record.
			organism.getGenes().get(i).setLci1(lci1);
			this.sred1 += lci1;
			this.disp1 += lci1*lci1;
			
			//Energy per nucleotide (subtract 3 for the stop codon)
			lci2 = energyFull/(organism.getGenes().get(i).getLength()-3);
			//Recording lci2
			organism.getGenes().get(i).setLci2(lci2);
			this.sred2 += lci2;
			this.disp2 += lci2*lci2;
		}
		
		this.sred1 /= organism.getGenes().size();
		this.disp1 = this.disp1/organism.getGenes().size() - this.sred1*this.sred1;
		this.disp1 = Math.pow(this.disp1, 0.5);
		
		//Fortran
		this.m1 = this.sred1 + 2*this.disp1;
		this.sigma1 = 4*this.disp1;
		
		//Articles
		//this.m1 = this.sred1;
		//this.sigma1 = forSigma1*this.disp1;
		
		this.sred2 /= organism.getGenes().size();
		this.disp2 = this.disp2/organism.getGenes().size() - this.sred2*this.sred2;
		this.disp2 = Math.pow(this.disp2, 0.5);
		
		//Fortran
		this.m2 = this.sred2 + this.disp2;
		this.sigma2 = 2*this.disp2;
		
		//Articles
		//this.m2 = this.sred2;
		//this.sigma2 = forSigma2*this.disp2;
	}
	
	//TODO Returns the average LCI
	public double getM1() {
		return this.m1;
	}
	
	//TODO Returns the mean deviation of the LCI
	public double getSigma1() {
		return this.sigma1;
	}
	
	//TODO Returns the average LCI
	public double getM2() {
		return this.m2;
	}
	
	//TODO Returns the mean deviation of the LCI
	public double getSigma2() {
		return this.sigma2;
	}
}
