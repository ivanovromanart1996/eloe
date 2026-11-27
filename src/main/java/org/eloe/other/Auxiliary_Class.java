package org.eloe.other;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import org.apache.log4j.Logger;

import org.eloe.engine.ta.FreqCounter;
import org.eloe.engine.ta.M1_Counter;
import org.eloe.organism.Gene;
import org.eloe.organism.Organism;
import org.eloe.parameters.PAutoCounting;
import org.eloe.parameters.PGeneToString;

public abstract class Auxiliary_Class {
	
	private static final Logger LOG = Logger.getLogger(Auxiliary_Class.class);
	
	public static int[] ribpozition(List<Gene> l, int ribgenes, File fileRibPos) {
		int[] djmass = {0, 0};
		
		//If there are ribosomal genes
		if (ribgenes > 0) {
			//Array of positions of ribosomal genes
			int ribPozMass[] = new int[ribgenes];
			int j = 0;
			//Center of the ribosomal gene group
			int ng = 0;
			for (int i = 0; i < l.size(); i++) {
				if (l.get(i).isRibosomal()) {
					ribPozMass[j] = i+1;
					j++;
					ng += (i+1);
				}
			}
			//Center of the ribosomal gene group
			double m = (double)ng/ribgenes;
			//Average deviation of ribosomal genes from the center
			double r = 0;
			for (int i = 0; i < ribgenes; i++) {
				r += Math.pow((ribPozMass[i]-m), 2);
			}
			r = Math.sqrt(r/ribgenes);
			
			//Center of the ribosomal gene group (scaled)
			int dj = (int)Math.round(100*(2*m-l.size()-1)/(l.size()-ribgenes));
			//Average deviation of ribosomal genes from the center (scaled)
			int djr = (int)Math.round(100*(2*r-ribgenes-1)/(l.size()-ribgenes));
			
			/*
			int dj = (int)Math.round(100*(2*(m-1)/(l.size()-1)-1));
			int djr = (int)Math.round(200*(r-1)/(l.size()-1));
			*/
			
			//Fortran
			/*double djr1 = 0;
			double dj1 = (2*m-l.size()-1)/(l.size()-1);
			double hh = 0;
			for (int i = 0; i < ribgenes; i++) {
				hh = (double)(2*(ribPozMass[i]-1)-(l.size()-1))/(l.size()-1);
				djr1 += Math.pow((hh-dj1), 2);
			}
			djr1 /= ribgenes;
			int dj = (int)Math.round(100*(2*m-l.size()-1)/(l.size()-1));
			int djr = (int)Math.round(100*Math.sqrt(djr1);*/
						
			djmass[0] = dj;
			djmass[1] = djr;
			
			//Record file with positions of ribosomal genes
			PrintStream outRibPoz = null;
			try {
				outRibPoz = new PrintStream(
						new BufferedOutputStream(
			               new FileOutputStream(fileRibPos, true)));
				for (int i = 0; i < ribPozMass.length; i++) {
					outRibPoz.print(ribPozMass[i] + "\t");
				}
				outRibPoz.print("\n");
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (outRibPoz != null) {
					outRibPoz.close();
				}
			}
		} else { //If there are no ribosomal genes
			PrintStream out;
			try {
				out = new PrintStream(fileRibPos);
				out.println("//");
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return djmass;
	}
	
	//TODO Reading from the ribosomal gene positions file
	public static int[] ribosomsPositions(String ribpoz) {
		String[] as = ribpoz.split("\t");
		int[] aa = new int[as.length];
		for (int i = 0; i < as.length; i++) {
			aa[i] = Integer.parseInt(as[i].trim());
		}
		return aa;
	}
	
	//TODO Creating Complementary Sequences
	public static char[] reversComplement(final char[] na) {
		//Gene length
		int geneLength = na.length;
		//Gene sequence
		char[] naReversComplement = new char[geneLength];
		
		for (int j = 0; j < geneLength; j++) {
			int q = geneLength-j-1;
			switch (na[j]) {
				case ('c'): naReversComplement[q] = 'g'; break;
				case ('a'): naReversComplement[q] = 't'; break;
				case ('g'): naReversComplement[q] = 'c'; break;
				case ('t'): naReversComplement[q] = 'a'; break;
				default: naReversComplement[q] = 'n';
			}
		}
		return naReversComplement;
	}
	
	//TODO Computation time
	public static long time(final String s, final long t) {
		long time = System.currentTimeMillis()-t;
		LOG.info("Operating period ("+s+"): "+time+" millisec");
		return System.currentTimeMillis();
	}
	
	//TODO Automatic flank selection
	public static void flangSizeAutocounting(PAutoCounting parameters, Organism organism, int number2, FreqCounter freq) {
		if (parameters.isAutocounting()) {			
			//To select the number of codons that are not considered on the flanks of genes
			int mmm[][] = new int [parameters.getMaxNumbCod()+1][2];
			
			for (int i = 0; i <= parameters.getMaxNumbCod(); i++) {
				freq.countFreq(organism.getGenes(), i, i);
				mmm[i] = M1_Counter.count_m1(organism, number2, parameters.getMinKc(), parameters.getVe(), i, i);
			}
			
			int max = -101;
			int maxind = 0;
			for (int i = 0; i < parameters.getMaxNumbCod(); i++) {
				if (max < mmm[i][0]) {
					max = mmm[i][0];
					maxind = i;
				} else {
					if ((max == mmm[i][0]) && (mmm[maxind][1] > mmm[i][1])) {
						maxind = i;
					}
				}
			}
			
			//Number of excluded codons in the flanks
			parameters.setUnUseCodonsL1(maxind);
			parameters.setUnUseCodonsR1(maxind);
			LOG.info("Automatical count of number of unrecorded codons on the genes edges. n = " + maxind);
		}
	}
	
	//TODO Determination of an index with a maximum value of M and a minimum value of R
	public static int indexCounter(int[][] mr) {
		int index = 0;
		int m = -101;
		for (int i = 0; i < 5; i++) {
			if (m < mr[i][0]) {
				m = mr[i][0];
				index = i;
			} else {
				if ((m == mr[i][0]) && (mr[index][1] > mr[i][1])) {
					index = i;
				}
			}
		}
		return index;
	}
	
	//TODO Plot ribosomal gene positions.
	public static void ribPosGraphic(File fileRibPoz, Organism organism) throws IOException {
		//Reading a file with positions of ribosomal genes
		BufferedReader readerRibPoz = new BufferedReader(new FileReader(fileRibPoz));
		String ribpoz = readerRibPoz.readLine();
		
		if (!ribpoz.contains("//")) {
			ribpoz = readerRibPoz.readLine();
			int[] aa = Auxiliary_Class.ribosomsPositions(ribpoz);
			ribpoz = readerRibPoz.readLine();
			int[] bb = Auxiliary_Class.ribosomsPositions(ribpoz);
			ribpoz = readerRibPoz.readLine();
			int[] cc = Auxiliary_Class.ribosomsPositions(ribpoz);
			ribpoz = readerRibPoz.readLine();
			int[] dd = Auxiliary_Class.ribosomsPositions(ribpoz);
			ribpoz = readerRibPoz.readLine();
			int[] ee = Auxiliary_Class.ribosomsPositions(ribpoz);
			ribpoz = readerRibPoz.readLine();
			readerRibPoz.close();
				
			//Creating a window with a graphical display of the location of ribosomal genes
			if (fileRibPoz.canRead()) {
				//Displaying the position of ribosomal genes
				final GraphOfRibPoz f = new GraphOfRibPoz("Positions of ribosomal genes in " + organism.getName(), organism.getGenesNumber(), aa, bb, cc, dd, ee);
				f.addWindowListener(new WindowAdapter() {
					public void windowClosing(final WindowEvent ev) {
						f.dispose();
					}
				});
			}
		}
	}
	
	//Workaround for displaying LCI profiles for viruses
	//Number\teei + ("number_eei")
	public static void numberEEItoFile(Organism organism, File resultFolder1) {
		PrintStream outGenesAndFlangs6 = null;
		int index = organism.getIndex();
		//Record a file with Number and EEI
		File fileGenesAndFlangs6 = new File(resultFolder1, organism.getName() + "_number_eei.txt");
		try {
			outGenesAndFlangs6 = new PrintStream(fileGenesAndFlangs6);
			outGenesAndFlangs6.println("number\tsingle\tfirst\tmiddle\tlast\tEEI"+(index+1));
			for (int i = 0; i < organism.getGenes().size(); i++) {
				StringBuilder s = new StringBuilder();
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
				s.append("\t");
				//stub for compatibility
				s.append(0.0);//organism.getGenes().get(i).getEeiMass()[index]
				outGenesAndFlangs6.println(s.toString());
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} finally {
			if (outGenesAndFlangs6 != null) {
				outGenesAndFlangs6.close();
			}
		}
	}
	
	//TODO Record the results in files
	public static void resultsToFiles(Organism organism, PGeneToString parameters, File resultFolder1, List<String> files) {
		int index = organism.getIndex();
		
		boolean out1 = false; //">Number\tLocus_tag\tFBid\tComplement\tEEi" (Drosophila for Levitsky) ("genes_and_flangs1")
		boolean out2 = true; //"#\tlocus_tag\tprotein_id\tgene_id\tgene_name\tCOG\tEEI" (for Anatoly Bragin and me) ("eei+(index+1)") Required because files depends on this output
		boolean out3 = true; //All, l.get(i).toString()\nl.get(i).getNA ("genes_and_flangs3")
		boolean out4 = false; //"locus_tag\tEEI1\tEEI2\tEEI3\tEEI4\tEEI5\tIndex=" ("locus_tag_eei12345")
		boolean out5 = false; //("fbid_locus_tag")
		boolean out6 = true; //Number\teei + ("number_eei")
		boolean out7 = false; //>locus_tag\nstrNA2 (yeast for Levitsky)
		boolean out8 = false; //locus_tag\EEI (yeast for Levitsky)
		boolean out9 = true; //>locus_tag|EEI\nNA (.fasta for Natasha, Orlov’s student)
					
		PrintStream outGenesAndFlangs1 = null;
		PrintStream outGenesAndFlangs2 = null;
		PrintStream outGenesAndFlangs3 = null;
		PrintStream outGenesAndFlangs4 = null;
		PrintStream outGenesAndFlangs5 = null;
		PrintStream outGenesAndFlangs6 = null;
		PrintStream outGenesAndFlangs7 = null;
		PrintStream outGenesAndFlangs8 = null;
		PrintStream outGenesAndFlangs9 = null;
		
		try {
			//1
			if (out1) {
				//Record a file with gene + flank sequences with their names and product names
				File fileGenesAndFlangs1 = new File(resultFolder1, organism.getName() + "_genes_and_flanks1.txt");
				outGenesAndFlangs1 = new PrintStream(fileGenesAndFlangs1);
			}
			//2
			if (out2) {
				//Record a file with gene + flank sequences with their names and product names
				File fileGenesAndFlangs2 = new File(resultFolder1, organism.getName() + "_eei" + (index + 1) + ".txt");
				outGenesAndFlangs2 = new PrintStream(fileGenesAndFlangs2);
				outGenesAndFlangs2.println("#\t##\tlocus_tag\tprotein_id\tgene_id\tgene_name\tCOG\tEEI");
				files.add(resultFolder1 + File.separator + organism.getName() + "_eei" + (index + 1) + ".txt");
			}			
			//3
			if (out3) {
				//Record a file with gene + flank sequences with their names and product names
				File fileGenesAndFlangs3 = new File(resultFolder1, organism.getName() + "_genes_and_flanks.txt");
				outGenesAndFlangs3 = new PrintStream(fileGenesAndFlangs3);
				//outGenesAndFlangs3.println("Number\tlocus_tag\tgene_poz\t/gene=\t/product=\t/note=\t/position=\t/lenght=\t/codons=\t/complement=\tGC-compoz=\t/EEI=");
			}
			//4
			if (out4) {
				//File recording with locus_tag and EEI
				File fileGenesAndFlangs4 = new File(resultFolder1, organism.getName() + "_locus_tag_eei12345.txt");
				outGenesAndFlangs4 = new PrintStream(fileGenesAndFlangs4);
				outGenesAndFlangs4.println("locus_tag\tEEI1\tEEI2\tEEI3\tEEI4\tEEI5\tIndex=" + (index + 1));
			}
			//5
			if (out5) {
				//File recording with FBid and locus_tag
				outGenesAndFlangs5 = new PrintStream(resultFolder1, organism.getName() + "_fbid_locus_tag.txt");
				outGenesAndFlangs5.println("Number\tlocus_tag\tFBid\tgene_poz_start\tgene_poz_stop\tcdsStart\tcdsStop\tcomplement\tEEI\tgeneNA");
			}
			//6
			if (out6) {
				//Record a file with Number and EEI
				File fileGenesAndFlangs6 = new File(resultFolder1, organism.getName() + "_number_eei.txt");
				outGenesAndFlangs6 = new PrintStream(fileGenesAndFlangs6);
				outGenesAndFlangs6.println("number\tsingle\tfirst\tmiddle\tlast\tEEI"+(index+1));
			}
			//7
			if (out7) {
				//Record a file with Number and EEI
				File fileGenesAndFlangs7 = new File(resultFolder1, organism.getName() + "_locus_tag_strNA2.txt");
				outGenesAndFlangs7 = new PrintStream(fileGenesAndFlangs7);
			}
			//8
			if (out8) {
				//Record a file with Number and EEI
				File fileGenesAndFlangs8 = new File(resultFolder1, organism.getName() + "_locus_tag_EEI.txt");
				outGenesAndFlangs8 = new PrintStream(fileGenesAndFlangs8);
				outGenesAndFlangs8.println("locus_tag\tEEI");
			}
			//9
			if (out9) {
				//Record a file with Number and EEI
				File fileGenesAndFlangs9 = new File(resultFolder1, organism.getName() + "_locus_tag_EEI_2.fasta");
				outGenesAndFlangs9 = new PrintStream(fileGenesAndFlangs9);
			}
			
			String loc = "";
			for (int i = 0; i < organism.getGenes().size(); i++) {
				//1
				if (out1) {
					if (!organism.getGenes().get(i).getLocusTag().equalsIgnoreCase(loc) && 
							organism.getGenes().get(i).isFlagOut()) {
						loc = organism.getGenes().get(i).getLocusTag();
						outGenesAndFlangs1.println(">" + organism.getGenes().get(i).getNumber() +
								"\t" + organism.getGenes().get(i).getLocusTag() + 
								"\t" + organism.getGenes().get(i).getFbID() +
								"\t" + (organism.getGenes().get(i).isComplement() == true ? "-" : "+") +
								"\t" + organism.getGenes().get(i).getEeiMass()[index]);
						outGenesAndFlangs1.println(organism.getGenes().get(i).getNaWithFlanks()); //Previously NA2
					}
				}
				
				//3
				if (out3) {
					outGenesAndFlangs3.print(organism.getGenes().get(i).toString(organism, parameters));
					outGenesAndFlangs3.print(organism.getGenes().get(i).getEeiMass()[index]+"\n");
					outGenesAndFlangs3.println(organism.getGenes().get(i).getUserNA(parameters));
				}
				
				//4
				if (out4) {
					outGenesAndFlangs4.println(organism.getGenes().get(i).getLocusTag() +
							"\t" + organism.getGenes().get(i).getEeiMass()[0] + "\t" + organism.getGenes().get(i).getEeiMass()[1] +
							"\t" + organism.getGenes().get(i).getEeiMass()[2] + "\t" + organism.getGenes().get(i).getEeiMass()[3] +
							"\t" + organism.getGenes().get(i).getEeiMass()[4]);
				}
				
				//5
				if (out5) {
					outGenesAndFlangs5.print(organism.getGenes().get(i));
					outGenesAndFlangs5.print(organism.getGenes().get(i).getEeiMass()[index] + "\t");
				}
				
				//6
				if (out6) {
					StringBuilder s = new StringBuilder();
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
					s.append("\t");
					s.append(organism.getGenes().get(i).getEeiMass()[index]);
					outGenesAndFlangs6.println(s.toString());
				}
				
				//7
				if (out7 && organism.getGenes().get(i).isFlagOut()) {
					outGenesAndFlangs7.println(">" + organism.getGenes().get(i).getLocusTag());
					outGenesAndFlangs7.println(organism.getGenes().get(i).getNaWithFlanks()); //Previously NA2
				}
				
				//8
				if (out8 && organism.getGenes().get(i).isFlagOut()) {
					outGenesAndFlangs8.println(organism.getGenes().get(i).getLocusTag() +
							"\t"+organism.getGenes().get(i).getEeiMass()[index]);
				}
				
			}
			
			//2
			if (out2) {
				Sorter.sortdown(organism);
				for (int i = 0; i < organism.getGenes().size(); i++) {
					StringBuilder s = new StringBuilder();
					s.append(i + 1);
					s.append("\t");
					s.append(organism.getGenes().get(i).getNumber());
					s.append("\t");
					s.append(organism.getGenes().get(i).getLocusTag());
					s.append("\t");
					s.append(organism.getGenes().get(i).getProteinID());
					s.append("\t");
					s.append(organism.getGenes().get(i).getGeneID());
					s.append("\t");
					s.append(organism.getGenes().get(i).getName());
					s.append("\t");
					s.append(organism.getGenes().get(i).getCog());
					s.append("\t");
					s.append(organism.getGenes().get(i).getEeiMass()[index]);
					outGenesAndFlangs2.println(s.toString());
				}
			}
			
			//9
			if (out9) {
				Sorter.sortdown(organism);
				for (int i = 0; i < organism.getGenes().size(); i++) {
					StringBuilder s = new StringBuilder();
					s.append(">").append(organism.getGenes().get(i).getLocusTag())
					.append("|").append(organism.getGenes().get(i).getEeiMass()[index])
					.append("\n").append(organism.getGenes().get(i).getNAString());
					outGenesAndFlangs9.println(s.toString());
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (outGenesAndFlangs1 != null) {
				outGenesAndFlangs1.close();
			}
			if (outGenesAndFlangs2 != null) {
				outGenesAndFlangs2.close();
			}
			if (outGenesAndFlangs3 != null) {
				outGenesAndFlangs3.close();
			}
			if (outGenesAndFlangs4 != null) {
				outGenesAndFlangs4.close();
			}
			if (outGenesAndFlangs5 != null) {
				outGenesAndFlangs5.close();
			}
			if (outGenesAndFlangs6 != null) {
				outGenesAndFlangs6.close();
			}
			if (outGenesAndFlangs7 != null) {
				outGenesAndFlangs7.close();
			}
			if (outGenesAndFlangs8 != null) {
				outGenesAndFlangs8.close();
			}
			if (outGenesAndFlangs9 != null) {
				outGenesAndFlangs9.close();
			}
		}
	}
	
	//TODO Record a file with 5 EEI indexes
	public static void outEEI(File resultFolder1, Organism organism) {
		File fileIndexEEI = new File(resultFolder1, organism.getName()+"_all.txt");
		PrintStream outEEI = null;
		try {
			outEEI = new PrintStream(fileIndexEEI);
			
			for (int ind = 0; ind < 5; ind++) {
				//Sort by descending EEI index.
				Sorter.sortupEEIN(organism, ind+1);
				outEEI.println("Index = "+(ind+1));
				outEEI.println("#\tNumber\tEEI\tlocus_tag\tName\tProt\tNote");
				for (int i = 0; i < organism.getGenes().size(); i++) {
					outEEI.printf((i+1)+"\t"+organism.getGenes().get(i).getNumber()+"\t%.4e\t%s\t%s\t%s\t%s", 
							organism.getGenes().get(i).getEeiMass()[ind], 
							organism.getGenes().get(i).getLocusTag(), 
							organism.getGenes().get(i).getName(), 
							organism.getGenes().get(i).getProteinName(), 
							organism.getGenes().get(i).getNote());
					if (organism.getGenes().get(i).isRibosomal()) 
						outEEI.printf("\trib");
					outEEI.printf("\n");
				}
				if (ind < 4) outEEI.printf("\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (outEEI != null) {
				outEEI.close();
			}			
		}
	}
}
