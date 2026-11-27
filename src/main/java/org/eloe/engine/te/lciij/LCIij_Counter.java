package org.eloe.engine.te.lciij;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import org.eloe.engine.ta.FreqCounter;
import org.eloe.organism.Organism;
import org.eloe.other.Sorter;
import org.eloe.parameters.PLCI2;

// Program for calculating the local complementarity index

public class LCIij_Counter {
	private static final String ENERGY = "lci2double.txt";
	
	private PLCI2 parameters;
	//private int forSigma2; //Sigma multiplier
	private double loopenergy[] = null;
	private Map<String, Double> map = new HashMap<String, Double>();
	private String resultFileLCIij; //File for recording the results of calculations LCIij
	
	//Constructor
	public LCIij_Counter(PLCI2 parameters, String resultFileLCIij) {
		this.parameters = parameters;
		//this.forSigma2 = parameters.forSigma2;
		this.resultFileLCIij = resultFileLCIij;
		String loopsenergy = "";
		BufferedReader reader = null;
		String line = new String();
		try {
			//Reading a file with energies of dinucleotides and loops
			reader = new BufferedReader(new InputStreamReader(FreqCounter.class.getClassLoader().getResourceAsStream(ENERGY)));
			while (!(line = reader.readLine()).startsWith("//")){
				String arr[] = line.split(" - ");
				this.map.put(arr[0], Double.parseDouble(arr[1]));
			}
			while (!((line = reader.readLine()).startsWith("////"))) {
				loopsenergy += line;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
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
		for (int i = 0; i < massivesize; i++) {
			loopenergy[i] = Double.parseDouble(loopsenergy1[i].trim());
		}
	}
	
	//TODO Calculation of LCIij
	public void countLCIij(Organism organism) {
		PrintStream outLCIij = null;
		try {
			outLCIij = new PrintStream(resultFileLCIij);
			if (organism.isEukaryota()) {
				outLCIij.println("#\tLCIij\t(eucaryota)");
			}
			else {
				outLCIij.println("#\tfirst&last\tfirst&!last\tmiddle\tlast&!first\tLCIij\t(procaryota)");
			}
			
			for (int i = 0; i < organism.getGenes().size(); i++) {
				int d = parameters.getLmax2()-parameters.getLmin2();
				int lenght = organism.getGenes().get(i).getLenghtWithFlanks(); //Length of the part of the gene under study 
				char gene[] = organism.getGenes().get(i).getNaWithFlanks(); //Get the gene sequence with flanks.
				int edge = lenght-2*parameters.getSmax2()-parameters.getLmax2()+1; //The boundary of the left edge of the sliding frame
				if (edge < 0) continue;
				double lciij[] = new double[lenght]; // LCI individual nucleotide
				for (int j = 0; j < edge; j++){
					double energy = 0; //Energy of the codon interaction
					for (int s = parameters.getSmin2(); s <= parameters.getSmax2(); s++){
						char con1[] = new char[s]; //Buffer for the repeat
						System.arraycopy(gene, j, con1, 0, s); //Copying the sequence
						double energy1 = 0; //Energy of the codon interaction
						for (int m = 0; m < s-1; m++) {
							String dinukleotid = con1[m]+""+con1[m+1];
							if (map.containsKey(dinukleotid)) {
								energy1 += map.get(dinukleotid); //Repeat interaction energy
							}
						}
						StringBuilder con2 = new StringBuilder(""); //Buffer for the complementary sequence
						for (int k = s-1; k > -1 ; k--){ //Complementary Sequence Building
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
						int lenght1 = s+d; //Length of search area
						char con3[] = new char[lenght1]; //Search site
						System.arraycopy(gene, j+s+parameters.getLmin2(), con3, 0, lenght1); //Extract the search region.
						String gene1 = String.valueOf(con3); //Search site
						int poz = 0;
						double loop1 = 0; //Loop energy
						int k = gene1.indexOf(con22, poz);
						while (k != -1) {
							int loopsize = parameters.getLmin2()+k-1;
							//Loop energy
							loop1 = loopenergy[(loopsize < loopenergy.length) ? loopsize : (loopenergy.length-1)];
							double e = energy1 - loop1;
							//If the energy is greater than the threshold of the existence of the loop, then take into account this loop
							if (e >= parameters.getMinloopEnergy()) {
								energy += e;
							}
							poz = k+1;
							k = gene1.indexOf(con22, poz);
						}
						for (int jj = j; jj < j+s; jj++) {
							lciij[jj] += 10*energy;
						}
					}
				}
				//Record the local complementarity index of an individual nucleotide
				StringBuilder s = new StringBuilder("");
				s.append(organism.getGenes().get(i).getNumber());
				if (!organism.isEukaryota()) {
					if (organism.getGenes().get(i).isFirstCistron() && organism.getGenes().get(i).isLastCistron()) {
						s.append("\t1\t0\t0\t0");
					}
					else {
						if (organism.getGenes().get(i).isMiddleCistron()) {
							s.append("\t0\t0\t1\t0");
						}
						else {
							if (organism.getGenes().get(i).isFirstCistron() && !organism.getGenes().get(i).isLastCistron()) {
								s.append("\t0\t1\t0\t0");
							}
							else {
								s.append("\t0\t0\t0\t1");
							}
						}
					}
				}
				for (int j = 0; j < lciij.length; j++) {
					s.append("\t");
					s.append(lciij[j]);
				}
				outLCIij.println(s.toString());
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (outLCIij != null) {
				outLCIij.close();
			}
		}
	}	
}
