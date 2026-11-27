package org.eloe.extractor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;
import org.eloe.other.Auxiliary_Class;
import org.eloe.parameters.PExtraction;

public class EMBFileExtractor extends Extractor {
	
	private static final Logger LOG = Logger.getLogger(EMBFileExtractor.class);
	
	public EMBFileExtractor(File orgFile, PExtraction parameters, File resultFolder) {
		super(orgFile, parameters, resultFolder);
	}

	@Override
	public final Organism extract() {
		Organism organism = new Organism();
		
		//Number of genes in the list before reading this record
		int lsize = organism.getGenes().size();
		
		//Reading a genome file
		BufferedReader reader = null;
	
		try {
			//Reading a genome file
			reader = new BufferedReader(new FileReader(orgFile));
			
			//Linear or circular DNA
			boolean linear = true;
			
			String line = reader.readLine();
			
			if (!line.startsWith("ID ")) {
				LOG.error("ERROR! Wrong file format or file is empty: " + orgFile + "! Right formats: *.gbk, *.gb, *.emb");
				return null;
			}
			
			//DNA linearity
			if (line.contains(" circular")) {
				linear = false;
			}
			String[] arr3 = line.split(" ");
			
			//Genome size
			int genomeSize = 0;
			for (int i = arr3.length - 1; i >= 0; i--) {
				if (arr3[i].toLowerCase().contains("bp")) {
					genomeSize = Integer.parseInt(arr3[i - 1].trim());
					organism.setGenomeSize(genomeSize);
					break;
				}
			}
			
			while (!(line = reader.readLine()).startsWith("DE ")) {}
			
			line = reader.readLine().trim();
			if (line.contains("plasmid")) {
				LOG.info("File with plasmid - " + orgFile.getName());
				return null;
			}
			
			//Flag to identify an organism
			boolean flagOrg = true;
			
			//Determination of the name of the organism
			do {
				line = reader.readLine().trim();
				
				//Reading the name of the organism
				if (line.startsWith("OS ")) {
					String organismName = line.replace("OS   ", "");
					organismName = organismName.replaceAll(ORGNAME, "");
					organismName = organismName.replace("( |\\.|_)+", "_");
					organism.setName(organismName.toLowerCase());
					flagOrg = false;
					
					//Is eucaryota or not (operons)
					line = reader.readLine().toLowerCase();
					if (line.contains("eukaryota")) {
						organism.setEukaryota(true);
					}
					//Is archaea
					if (line.contains("archaea")) {
						organism.setArchaea(true);
					}
				}
			} while (flagOrg);
			
			//============================================================================================
	        //Read gene information from the record===================================================
	        //TODO========================================================================================
			
			//Flag for gene detection
			boolean flagGen = false;
			//Flag to identify pseudogene
			boolean flagPseudoGen = false;
			//Flag to identify gene product
			boolean flagProduct = false;
			//Flag to determine protein index
			boolean flagProtID = false;
			//Flag to determine gene index
			boolean flagDBxref = false;
			//Flag for annotation
			boolean flagNote = false;
			//Flag for gene index
			boolean flagLocus = false;
			//Flag for determining the operon number (needed only on the first pass)
			boolean flagOper = false;
			//Position of the end of the previous gene
			int lastpoz = 0;
			//Strand of the previous gene
			boolean lastpozcomplement = false;
			String context_new = "FT   " + parameters.getContext();
			
			while (!(line = reader.readLine()).startsWith("SQ   ")) {
				line = line.trim();
				
				//Gene reading
				if (line.startsWith(context_new)) {
					//Number of genes, sequence number
					organism.addGene();
					//Number of operons
					if (!organism.isEukaryota()) {
						organism.addOperon();
					}
					
					//Flag to determine the presence of introns
					boolean flagIntr = true;
					//Whether the location contains "join"
					boolean join = false;
					//Whether the location contains "complement"
					boolean complement = false;
					
					//Record of positions (gene) of exons
					String position = line;
					
					//Information on complementarity
					if (position.contains("complement")) {
						complement = true;
						position = position.replace(context_new + "            complement", "");
						if (position.contains("join")) {
							join = true;
							//Intron presence information
							position = position.replace("join", "");
						}
						flagIntr = false;
					}
					
					//Intron presence information
					if (flagIntr && position.contains("join")) {
						join = true;
						position = position.replace(context_new + "            join", "");
					}
					
					while ((complement || join) && !(line = reader.readLine()).startsWith("FT                   /")){
						line = line.replace("FT                   ", "");
						if (line.contains("complement")) {
							line = line.replace("complement", "");
						}
						position += line.trim();
					}
					
					if (!(complement || join)) {
						position = line.replace(context_new + "            ", "");
					}
					
					//Remove parentheses and angle brackets.
					position = position.replaceAll(BRACKETS.toString(), "");
					
					//If there are several exons in the gene
					if (position.contains(",")) {
						
						//Array of exon positions
						String[] pozOfExons = position.split(",");
						
						//Number of exons
						int numberOfExons = pozOfExons.length;
						
						//Array of exon start positions
						int[] begining = new int[numberOfExons];
						
						//Array of exon end positions
						int[] end = new int[numberOfExons];
						for (int i = 0; i < numberOfExons; i++) {
							
							if (pozOfExons[i].contains("..")) {
								String arr[] = pozOfExons[i].split("\\.\\.");
								begining[i] = Integer.parseInt(arr[0]);
								end[i] = Integer.parseInt(arr[1]);
							} else { //If the exon consists of a single nucleotide
								//Positions of a single nucleotide in an exon
								begining[i] = Integer.parseInt(pozOfExons[i]);
								end[i] = Integer.parseInt(pozOfExons[i]);
							}
						}
						
						//Checking the distance between neighboring genes
						if (flagOper && !organism.isEukaryota() && !(lastpozcomplement ^ complement) && (begining[0] - lastpoz - parameters.getDistBOper() - 2) < 0) {
							organism.delOperon();
						}
						
						//Position of the end of the gene
						lastpoz = end[numberOfExons - 1];
						lastpozcomplement = complement;
						
						//Needed on the first pass.
						flagOper = true;
						Gene gene = new Gene(organism.getGenesNumber(), begining, end, organism.getOperonsNumber(), complement, true);
						
						//Putting a new gene on the list
						organism.getGenes().add(gene);
					} else { //If one exon in the gene
						String[] arr = position.split("\\.\\.");
						
						//Array containing the start position of the first exon
						int[] begining = {Integer.parseInt(arr[0])};
						
						//Array containing the end position of the first exon
						int[] end = {Integer.parseInt(arr[1])};
						
						//Checking the distance between neighboring genes
						if (flagOper && !organism.isEukaryota() && !(lastpozcomplement ^ complement) && (begining[0] - lastpoz - parameters.getDistBOper() - 2) < 0) {
							organism.delOperon();
						}
						
						//Position of the end of the gene
						lastpoz = end[0];
						lastpozcomplement = complement;
						
						//Needed on the first pass.
						flagOper = true;
						Gene gene = new Gene(organism.getGenesNumber(), begining, end, organism.getOperonsNumber(), complement, true);
						
						//Putting a new gene on the list
						organism.getGenes().add(gene);
					}
					
					flagGen = true;
					flagProduct = true;
					flagProtID = true;
					flagDBxref = true;
					flagNote = true;
					flagLocus = true;
					flagPseudoGen = true;
				}
				
				//Definition of a pseudogene
				if (flagPseudoGen && line.startsWith("FT                   /pseudo")) {
					organism.addPseudoGene();
					
					//The pseudogene mark
					organism.getGenes().get(organism.getGenes().size() - 1).setPseudo(true);
					flagPseudoGen = false;
				}
				
				//Name of the gene
				if (flagGen && line.startsWith("FT                   /gene=")) {
					line = line.replace("FT                   ", "");
					
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMALG.matcher(line).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Number of ribosomal genes
							organism.addRibGene();
						}
					}
					
					//Determination of highly expressed genes by label
					if (parameters.isUseOtherGenesAsHighExpressed1() && line.contains(" highexpressed")) {
						line = line.replace(" highexpressed", "");
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Number of ribosomal genes
							organism.addRibGene();
						}
					}
					
					//The name of the gene
					organism.getGenes().get(organism.getGenes().size() - 1).setName(line.replace("FT                   ", ""));
					flagGen = false;
				}
				
				//Gene index locus_tag=
				if (flagLocus && line.startsWith("FT                   /locus_tag=")) {
					//32 is the length of "FT                   /locus_tag="
					organism.getGenes().get(organism.getGenes().size() - 1).setLocusTag(line.substring(32).replace("\"", ""));
					flagLocus = false;
				}
				
				//Name of gene product
				if (flagProduct && line.startsWith("FT                   /product=")) {
					if (line.contains("COG")) {
						Matcher m2 = COGID.matcher(line);
						if (m2.find() && organism.getGenes().size() != 0) {
							organism.getGenes().get(organism.getGenes().size() - 1).setCog(m2.group());
						}
					}
					
					line = line.replace("FT                   ", "");
					line = line.toLowerCase();
					
					//Recording the name of the gene product
					organism.getGenes().get(organism.getGenes().size() - 1).setProteinName(line);
					
					//Determination of ribosomal genes by name of gene product
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line).matches() && !NONRIBOSOMAL.matcher(line).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Number of ribosomal genes
							organism.addRibGene();
						}
					}
					flagProduct = false;
					flagGen = false;
				}
				
				//protein_id=
				if (flagProtID && line.startsWith("FT                   /protein_id=")) {
					//33 is the length of "FT                   /protein_id="
					organism.getGenes().get(organism.getGenes().size() - 1).setProteinID(line.substring(33).replace("\"", ""));
					flagProtID = false;
				}
				
				//db_xref=\"GI:
				if (flagDBxref && line.startsWith("FT                   /db_xref=\"GI:")) {
					//30 is the length of "FT                   /db_xref="
					organism.getGenes().get(organism.getGenes().size() - 1).setGeneID(line.substring(30).replace("\"", ""));
					flagDBxref = false;
				}
				
				//Recording of annotations
				if (flagNote && line.startsWith("FT                   /note=")) {
					if (line.contains("COG")) {
						Matcher m2 = COGID.matcher(line);
						if (m2.find() && organism.getGenes().size() != 0) {
							organism.getGenes().get(organism.getGenes().size() - 1).setCog(m2.group());
						}
					}
					String line1 = line;
					/*while (!line.startsWith("FT                   /")){
						line1 += line.replace("FT                   ", "").trim()+" ";
						line = reader.readLine();
					}*/
					organism.getGenes().get(organism.getGenes().size() - 1).setNote(line1.substring(27).replace("\"", "")); //27 is the length of "FT                   /note="
					
					//Determination of ribosomal genes by gene annotation
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line1).matches() && !NONRIBOSOMAL.matcher(line1).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							
							//Number of ribosomal genes
							organism.addRibGene();
						}
					}
					flagNote = false;
				}
				
				//COG
				if (line.contains("COG")) {
					Matcher m2 = COGID.matcher(line);
					if (m2.find() && organism.getGenes().size() != 0) {
						organism.getGenes().get(organism.getGenes().size() - 1).setCog(m2.group());
					}
				}
				
				if (flagProduct && line.startsWith("FT   gene          ")) {
					flagPseudoGen = false;
					flagGen = false;
					flagProduct = false;
					flagProtID = false;
					flagDBxref = false;
				}
				
				if (line.startsWith("FT   tRNA          ")) {
					organism.addtRnaGene();
				}
			}
			
			//============================================================================================
			//Read the complete nucleotide sequence from the record===============================
			//TODO========================================================================================
			
			//An array to record a complete genome
			char[] wholeGenome = new char[genomeSize];
			
			int pozstrGenom = 0;
			
			//Number of GC and CG pairs, GC-composition
			int gc = 0;
			
			//Genomic sequence reading
			while (!(line = reader.readLine()).startsWith("//")) {
				line = line.trim();
				String[] arr = line.split(" ");
				
				//Line for recording genomic sequence
				StringBuilder strGenom = new StringBuilder();
				
				for (int i = 0; i < arr.length - 1; i++) {
					strGenom.append(arr[i]);
				}
				
				char[] ch = strGenom.toString().toCharArray();
				
				int size = ch.length;
				
				for (int i = 0; i < size; i++) {
					if (ch[i] == 'g' || ch[i] == 'c') {
						gc++;
					}
				}
				
				System.arraycopy(ch, 0, wholeGenome, pozstrGenom, size);
				pozstrGenom += size;
			}
			
			//GC composition
			organism.addGC(gc);
			double gccom = Math.round(gc*10000.0/genomeSize)/100.0;
			int lFlang1;
			int rFlang1;
			
			//Recording gene sequences in instances of a class
			for (int i = lsize; i < organism.getGenes().size(); i++) {
				if (organism.getGenes().get(i).isComplement()) {
					lFlang1 = parameters.getrFlang();
					rFlang1 = parameters.getlFlang();
				} else {
					lFlang1 = parameters.getlFlang();
					rFlang1 = parameters.getrFlang();
				}
				
				//Number of introns
				int k = organism.getGenes().get(i).getStart().length;
				
				//The position of the gene onset
				int strt = organism.getGenes().get(i).getStart()[0];
				
				//Position of the end of the gene
				int stp = organism.getGenes().get(i).getStop()[k - 1];
				int stp1 = stp + rFlang1 - genomeSize;
				int dlinaGenaBezFl = 0;
				int dlinaGena;
				
				//A string to record a gene sequence
				char[] strNA;
				
				//Determination of gene length without flanks
				for (int j = 0; j < k; j++) {
					dlinaGenaBezFl += organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
				}
				
				//If DNA is linear
				if (linear) {
					//If the left flank is longer than it is on the chromosome before the gene begins
					if (strt < lFlang1) {
						lFlang1 = strt - 1;
						organism.getGenes().get(i).getStart()[0] = 1;
						organism.getGenes().get(i).setFlagUse(false);
					} else { //If the left flank is smaller than it is on the chromosome before the gene begins
						organism.getGenes().get(i).getStart()[0] -= lFlang1;
					}
					
					//If the right flank is longer than it is on the chromosome after the end of the gene
					if (stp1 > 0) {
						rFlang1 = genomeSize - stp;
						organism.getGenes().get(i).getStop()[k - 1] = genomeSize;
						organism.getGenes().get(i).setFlagUse(false);
					} else { //If the length of the right flank is smaller than it is on the chromosome after the end of the gene
						organism.getGenes().get(i).getStop()[k - 1] += rFlang1;
					}
					
					//Length of gene with flanks
					dlinaGena = lFlang1 + dlinaGenaBezFl + rFlang1;
					
					//Array for recording gene sequences with flanks
					strNA = new char [dlinaGena];
					int pozEx = 0;
					
					//Reading a gene sequence from a complete genome
					for (int j = 0; j < k; j++) {
						int dlinaEx = organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
						
						//Gene sequence string
						System.arraycopy(wholeGenome, organism.getGenes().get(i).getStart()[j] - 1, strNA, pozEx, dlinaEx);
						pozEx += dlinaEx;
					}
				} else { //If the DNA is circular
					//Does the right flank fit?
					boolean rflng = false; 
					
					//Length of gene with flanks
					dlinaGena = lFlang1 + dlinaGenaBezFl + rFlang1;
					
					//Array for recording gene sequences with flanks
					strNA = new char [dlinaGena];
					int pozEx = 0;
					
					//If the left flank is smaller than it is on the chromosome before the gene begins
					if (strt > lFlang1) {
						
						//New starting position with the flanks in mind
						organism.getGenes().get(i).getStart()[0] -= lFlang1;
					} else { //If the left flank is longer than it is on the chromosome before the gene begins
						//Length of missing sequence to full flank length
						int lenght = lFlang1 - strt + 1;
						
						//Recording the missing sequence from the other (right) edge of the ring chromosome
						//Gene sequence string
						System.arraycopy(wholeGenome, genomeSize - lenght, strNA, 0, lenght);
						organism.getGenes().get(i).getStart()[0] = 1;
						pozEx = lenght;
					}
					
					//If the length of the right flank is smaller than it is on the chromosome after the end of the gene
					if (stp1 < 0) {
						//New end position, taking into account the flanks;
						organism.getGenes().get(i).getStop()[k - 1] += rFlang1;
					} else { //If the right flank is longer than it is on the chromosome after the end of the gene
						organism.getGenes().get(i).getStop()[k - 1] = genomeSize;
						rflng = true;
					}
					
					//Reading a gene sequence from a complete genome
					for (int j = 0; j < k; j++) {
						int dlinaEx = organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
						
						//Gene sequence string
						System.arraycopy(wholeGenome, organism.getGenes().get(i).getStart()[j] - 1, strNA, pozEx, dlinaEx); 
						pozEx += dlinaEx;
					}
					
					//If the right flank does not fit on the chromosome, the missing sequence is recorded from the other (left) edge of the ring chromosome.
					if (rflng) {
						//Gene sequence string
						System.arraycopy(wholeGenome, 0, strNA, pozEx, stp1); 
					}
				}
				
				//If the gene is on the complementary strand
				if (organism.getGenes().get(i).isComplement()) {
					//Creating Complementary Sequences
					strNA = Auxiliary_Class.reversComplement(strNA);
				}
				
				//The position of the beginning of the gene without taking into account the flanks
				organism.getGenes().get(i).getStart()[0] = strt;
				
				//The position of the end of the gene without taking into account flanks
				organism.getGenes().get(i).getStop()[k - 1] = stp;
				
				//Gene sequence recording
				organism.getGenes().get(i).setNA(strNA);
				
				if (organism.getGenes().get(i).isComplement()) {
					int z = lFlang1;
					lFlang1 = rFlang1;
					rFlang1 = z;
				}
				
				//Recording the length of a gene with flanks
				organism.getGenes().get(i).putLenght(lFlang1, dlinaGenaBezFl, rFlang1, dlinaGena);
				
				//Marking to remove short genes
				if (dlinaGenaBezFl < parameters.getMinLenGen()) {
					organism.getGenes().get(i).setFlagUse(false);
					organism.addShortGene();
					//System.out.println("Gene number "+organism.genes.get(i).getNumber()+" ("+organism.genes.get(i).getName()+") was deleted because of its small length!");
				}
				
				//Marking to remove pseudogenes
				if (!parameters.isUsePseudoGenes() && organism.getGenes().get(i).isPseudo()) {
					organism.getGenes().get(i).setFlagUse(false);
				}
			}
			
			//Creating a file with complete genomic sequence
			printWholeGenome(organism, wholeGenome, gccom);
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
		
		return organism;
	}
}
