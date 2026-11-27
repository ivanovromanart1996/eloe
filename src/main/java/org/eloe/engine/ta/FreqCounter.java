package org.eloe.engine.ta;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;

//Counts the codon frequencies.

public class FreqCounter {
	
	private static final Logger LOG = Logger.getLogger(FreqCounter.class);
	
	private static final String SECOFCODONS = " \tcga,\tcgg,\tcgc,\tcgt,\taga,\tagg,\tcta,\tctg,\tctc,\tctt,"
			+ "\ttta,\tttg,\ttca,\ttcg,\ttcc,\ttct,\tagc,\tagt,\taca,\tacg,\tacc,\tact,\tcca,\tccg,\tccc,"
			+ "\tcct,\tgca,\tgcg,\tgcc,\tgct,\tgga,\tggg,\tggc,\tggt,\tgta,\tgtg,\tgtc,\tgtt,\taaa,\taag,"
			+ "\taac,\taat,\tcaa,\tcag,\tcac,\tcat,\tgaa,\tgag,\tgac,\tgat,\ttac,\ttat,\ttgc,\ttgt,\tttc,"
			+ "\tttt,\tatc,\tatt,\tata,\tatg,\ttgg,\ttga,\ttaa,\ttag";
	private static final String CODE_TXT = "code.txt";
	
	//Map for codons and their positions in the array
	private static Map<String, Integer> mapOfCodons;
	private ICodonCompositionCalculator ccCalc;
	
	//Constructor
	public FreqCounter() throws IOException, NumberFormatException {
		mapOfCodons = new HashMap<String, Integer>();
		//Filling the map with 64 codons and their positions in the array
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(FreqCounter.class.getClassLoader()
							.getResourceAsStream(CODE_TXT)));
			String line;
			while ((line = reader.readLine()) != null) {
				String[] arr = line.split(" - ");
				mapOfCodons.put(arr[0], Integer.parseInt(arr[1]));
			}
		} catch (NumberFormatException e) {
			LOG.error("Incorrect number in file "+CODE_TXT+" for map of codons.\n"+e);
			throw e;
		} catch (IOException e) {
			LOG.error(e);
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
		ccCalc = new CodonCompositionCalculatorImpl(mapOfCodons);
	}
	
	//TODO Calculation of codon frequencies in a gene with edge codons
	public void countFreq0(List<Gene> genes) {
		for (Gene gene : genes) {
			gene.setCodonComposition0(ccCalc.calculate(gene, 0, 0));
			gene.putGC_compoz();
		}
	}

	//TODO Calculation of codon frequencies in a gene without edge codons
	public void countFreq(List<Gene> genes, int unUsedCodonsL, int unUsedCodonsR) {
		for (Gene gene : genes) {
			gene.setCodonComposition(ccCalc.calculate(gene, unUsedCodonsL, unUsedCodonsR));
		}
	}
	
	public void makeFileFreq(List<Gene> l, String organismName, File resultFolder) {
		PrintStream outFreq = null;
		try {
			outFreq = new PrintStream(new File(resultFolder, organismName + "_freq.txt"));
			outFreq.print(">"+organismName+"\n"+SECOFCODONS+"\n");
			for (int i = 0; i < l.size(); i++) {
				outFreq.println(l.get(i).getNumber()+"\t"+l.get(i).getCodonComposition0String());
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (outFreq != null)
				outFreq.close();
		}
	}
}
