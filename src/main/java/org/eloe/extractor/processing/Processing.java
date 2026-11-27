package org.eloe.extractor.processing;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;
import org.eloe.parameters.PProcessing;

//Removal from the list of “bad” genes
public abstract class Processing {
	
	private static final Logger LOG = Logger.getLogger(Processing.class);
	
	public static void processing (Organism organism, PProcessing parameters) {
		int badGenes = 0;
		int ribGenes = 0;
		
		//Remove genes marked with "n".
		for (int i = 0; i < organism.getGenes().size(); i++) {
			if (organism.getGenes().get(i).isFlagUse()) {
				if (organism.getGenes().get(i).getNAString().contains("n")) {
					//Marking a gene as unused
					organism.getGenes().get(i).setFlagUse(false);
					if (organism.getGenes().get(i).isRibosomal()) {
						ribGenes++;
						organism.delRibGene();
					}
					badGenes++;
				}
			}
		}
		
		if (parameters.isCheckStart() || parameters.isDeleteGwBC() || parameters.isCheckStop()) {
			//Filling the map with start codons
			Set<String> setStart = new HashSet<String>(parameters.getStartCodons().length);
			for (int j = 0; j < parameters.getStartCodons().length; j++) {
				setStart.add(parameters.getStartCodons()[j].trim());
			}
			//Filling the map with bad codons
			Set<String> setMiddle = new HashSet<String>(parameters.getMidlCodons().length);
			for (int j = 0; j < parameters.getMidlCodons().length; j++) {
				setMiddle.add(parameters.getMidlCodons()[j].trim());
			}
			//Filling the map with stop codons
			Set<String> setStop = new HashSet<String>(parameters.getStopCodons().length);
			for (int j = 0; j < parameters.getStopCodons().length; j++) {
				setStop.add(parameters.getStopCodons()[j].trim());
			}
			
			for (int i = 0; i < organism.getGenes().size(); i++) {
				if (organism.getGenes().get(i).isFlagUse()) {
					//Gene length without flanks
					int lenght = organism.getGenes().get(i).getLength();
					//For gene sequence without flanks
					char gene[] = organism.getGenes().get(i).getNA();
					
					char c[] = new char[3];
					String codon = "";
					
					if (parameters.isCheckStart()) {
						System.arraycopy(gene, 0, c, 0, 3);
						codon = String.valueOf(c);
						//Testing for the presence of a starting codon at the beginning of the gene
						if (!setStart.contains(codon)) {
							//Marking a gene as unused
							organism.getGenes().get(i).setFlagUse(false);
							if (organism.getGenes().get(i).isRibosomal()) {
								ribGenes++;
								organism.delRibGene();
							}
							badGenes++;
						}
					}
					
					if (parameters.isDeleteGwBC()) {
						//Number of codons in a gene without the last
						int codons = organism.getGenes().get(i).getCodonsNumber()-1;
						for (int j = 1; j < codons; j++) {
							System.arraycopy(gene, 3*j, c, 0, 3);
							codon = String.valueOf(c);
							//Checking for a bad codon in the middle of a gene
							if (setMiddle.contains(codon)) {
								//Marking a gene as unused
								organism.getGenes().get(i).setFlagUse(false);
								if (organism.getGenes().get(i).isRibosomal()) {
									ribGenes++;
									organism.delRibGene();
								}
								badGenes++;
								break;
							}
						}
					}
					
					if (parameters.isCheckStop()) {
						System.arraycopy(gene, lenght-3, c, 0, 3);
						codon = String.valueOf(c);
						//Check for the presence of a stop codon at the end of the gene
						if (!setStop.contains(codon)) {
							//Marking a gene as unused
							organism.getGenes().get(i).setFlagUse(false);
							if (organism.getGenes().get(i).isRibosomal()) {
								ribGenes++;
								organism.delRibGene();
							}
							badGenes++;
						}
					}
				}
			}
		}
		
		//Removal from the list of all tagged genes
		int counterDeletGenes = 0;
		for (int i = 0; i < organism.getGenes().size(); i++) {
			if (!organism.getGenes().get(i).isFlagUse()) {
				organism.getGenes().remove(i);
				counterDeletGenes++;
				i--;
			}
		}
		((ArrayList<Gene>) organism.getGenes()).trimToSize();
		
		int pseudogenes = organism.getPseudoGenesNumber();
		if (parameters.isUsePseudoGenes()) {
			pseudogenes = 0;
		}
		//Reporting the number of removed ribosomal genes
		//Report that short genes have been removed
		//Report that pseudogenes have been removed
		//Report that bad genes have been removed
		//Report of total number of deleted genes
		//Reporting the number of genes involved in the calculations
		LOG.info(ribGenes+" ribosomal genes were deleted!\n"
				+organism.getShortGenesNumber()+" short genes were deleted!\n"
				+pseudogenes+" pseudogenes were deleted!\n"
				+badGenes+" bad genes were deleted!\n"
				+ "In all "+counterDeletGenes+" genes were deleted!\n"
				+ "Number of genes used in counting: "+organism.getGenes().size());
	}
}
