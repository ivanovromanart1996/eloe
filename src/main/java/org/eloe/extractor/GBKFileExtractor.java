package org.eloe.extractor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;
import org.eloe.other.Auxiliary_Class;
import org.eloe.parameters.PExtraction;

public class GBKFileExtractor extends Extractor {
	
	//2025.10.15 Stores the gene feature’s /locus_tag, /note and ribosomal status.
	private static class GeneNote {
		public String locus_tag;
		public String note;
		public boolean isRib;
		
		public GeneNote() {
			this.reset();
		}
		
		public void reset() {
			this.locus_tag = null;
			this.note = null;
			this.isRib = false;
		}
	}
	
	private static final Logger LOG = Logger.getLogger(GBKFileExtractor.class);
	private GeneNote geneNote;
	
	public GBKFileExtractor(File orgFile, PExtraction parameters, File resultFolder) {
		super(orgFile, parameters, resultFolder);
		this.geneNote = new GeneNote();
	}
	
	@Override
	public final Organism extract() {
		Organism organism = new Organism();
		
		//Reading a genome file
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(orgFile));
			String line = reader.readLine();
			
			if (!line.startsWith("LOCUS ")) {
				LOG.error("ERROR! Wrong file format or file is empty: " + orgFile.getName() + "! Right formats: *.gbk, *.gb, *.emb");
				return null;
			}
			
			//Linear or circular DNA
			boolean linear = true;
			if (line.contains(" circular ")) {
				linear = false;
			}
			
			String[] arr3 = line.split("\\s+");
			//Record name
			this.card = arr3[1];
			//Genome size
			int genomeSize = Integer.parseInt(arr3[2]);
			organism.setGenomeSize(genomeSize);

			line = reader.readLine();
			//If file contains plasmid NA
			if (line.contains("plasmid")) {
				LOG.info("File with plasmid - " + orgFile.getName());
				return null;
			}
			
			//Determination of the name of the organism
			boolean flagOrg = true;
			do {
				line = reader.readLine().trim();
				//Reading the name of the organism
				if (line.startsWith("ORGANISM")) {
					String organismName = line.replace("ORGANISM  ", "");
					//Formatting the name of the organism so that it can be used in the file name
					organismName = organismName.replaceAll(ORGNAME, "");
					organismName = organismName.replaceAll("ORGNAME_", "_");
					organism.setName(organismName.toLowerCase());
					flagOrg = false;
					
					//Is eucaryota or not (operons)
					line = reader.readLine().toLowerCase();
					if (line.contains("eukaryota")) {
						organism.setEukaryota(true);
					} else if (line.contains("archaea")) {//Is archaea
						organism.setArchaea(true);
					}
				}
			} while (flagOrg);
			
			//Start point of genes positions
			while (!(line = reader.readLine()).startsWith("FEATURES ")) {}
			
			//============================================================================================
	        //Read gene information from the record===================================================
	        //TODO========================================================================================
			
			//Gene Detection Flag (CDS)
			boolean flagCDS = false;
			//Flag to identify pseudogene
			boolean flagPseudoGene = false;
			//Flag to identify gene product
			boolean flagProduct = false;
			//Flag to determine protein index
			boolean flagProtID = false;
			//Flag to determine gene index
			boolean flagDBxref = false;
			//Flag for identifying the FlyBase gene ID
			//boolean flagFBid = false;
			//Flag for reading the gene note
			boolean flagNote = false;
			//Flag for reading a multiline gene note
			boolean flagNoteMultistring = false;
			//Flag for gene index
			boolean flagLocus = false;
			//Flag of the fact that we are in context (CDS)
			//UPD:2025.10.15 Commented out because the CDS belongs to the gene feature.
			//boolean flagContext = false;
			//Position of the end of the previous gene
			int lastpoz =  -2-parameters.getDistBOper();
			//Strand of the previous gene
			boolean lastpozcomplement = false;
			
			while (!(line = reader.readLine()).startsWith("ORIGIN")) {
				line = line.trim();
				
				//Gene Information Reading (CDS)
				if (line.startsWith(parameters.getContext() + "   ")) {
					//2025.10.16 If this is the very first CDS, then the code below should be skipped
					int num = organism.getGenes().size();
					//2025.10.16 The block below is needed only for genes with multiple CDS features. Otherwise,
					//the locus tags will not match, as expected.
					if (0 != num) {
						//2025.10.16 A new CDS feature has started. Save the information collected 
						//while processing the gene feature associated with the entire series of CDS features.
						//Save the /note and ribosomal status.
						Gene g = organism.getGenes().get(num - 1);
						if (null != this.geneNote.note && this.geneNote.locus_tag.equals(g.getLocusTag())) {
							g.setNote(this.geneNote.note);
							if (this.geneNote.isRib && !g.isRibosomal()) {
								g.setRibosomal(true);
								//Number of ribosomal genes
								organism.addRibGene();
							}
						}
					}
					//flagContext = true;
					//Number of genes, sequence number
					organism.addGene();
					//Number of operons
					if (!organism.isEukaryota()) {
						organism.addOperon();
					}
					
					boolean useGeneInCounting = true;
					//Whether the location contains "complement"
					boolean complement = false;
					//Number of exons in the gene
					int numberOfExons = 0;
					
					//Record of positions (gene) of exons
					String position = line.replace(parameters.getContext() + "            ", "");
					//Reading all gene positions
					while (!(line = reader.readLine().trim()).startsWith("/")) {
						position += line;
					}
					//Removal of brackets () and characters <>
					position = position.replaceAll(BRACKETS, "");
					
					//Array of exon start positions
					int[] start = null;
					//Array of exon end positions
					int[] stop = null;
					
					//If the position begins with the "join" identifier
					if (position.startsWith("join")) {
						if (position.contains("complement")) {
							position = position.replace("complement", "");
							complement = true;
							useGeneInCounting = false;
						}
						
						position = position.replace("join", "");
						
						if (position.contains(",")) {
							//Array of exon positions
							String[] pozOfExons = position.split(",");
							//Number of exons
							numberOfExons = pozOfExons.length;
							//Array of exon start positions
							start = new int[numberOfExons];
							//Array of exon end positions
							stop = new int[numberOfExons];
							
							for (int i = 0; i < numberOfExons; i++) {
								//If both the beginning and the end are indicated
								if (pozOfExons[i].contains("..")) {
									String[] arr = pozOfExons[i].split("\\.\\.");
									start[i] = Integer.parseInt(arr[0]);
									stop[i] = Integer.parseInt(arr[1]);
								}
								//If the position consists of one nucleotide
								else {
									start[i] = Integer.parseInt(pozOfExons[i]);
									stop[i] = Integer.parseInt(pozOfExons[i]);
								}
							}
						}
					} 
					//If the position begins with the "complement" identifier
					else {
						if (position.startsWith("complement")) {
							position = position.replace("complement", "");
							complement = true;
							
							//If the positions contain the "join" identifier
							if (position.startsWith("join")) {
								position = position.replace("join", "");
								
								if (position.contains(",")) {
									//Array of exon positions
									String pozOfExons[] = position.split(",");
									//Number of exons
									numberOfExons = pozOfExons.length;
									//Array of exon start positions
									start = new int[numberOfExons];
									//Array of exon end positions
									stop = new int[numberOfExons];
									
									for (int i = 0; i < numberOfExons; i++) {
										//If both the beginning and the end are indicated
										if (pozOfExons[i].contains("..")) {
											String arr[] = pozOfExons[i].split("\\.\\.");
											start[i] = Integer.parseInt(arr[0]);
											stop[i] = Integer.parseInt(arr[1]);
										}
										//If the position consists of one nucleotide
										else {
											start[i] = Integer.parseInt(pozOfExons[i]);
											stop[i] = Integer.parseInt(pozOfExons[i]);
											useGeneInCounting = false;
										}
									}
								}
							} 
							//If the positions do not contain the "join" identifier
							else {
								//Number of exons
								numberOfExons = 1;
								//Array containing the start position of the first exon
								start = new int[numberOfExons];
								//Array containing the end position of the first exon
								stop = new int[numberOfExons];
								
								//If both the beginning and the end are indicated
								if (position.contains("..")) {
									String[] arr = position.split("\\.\\.");
									start[0] = Integer.parseInt(arr[0]);
									stop[0] = Integer.parseInt(arr[1]);
								} 
								//If the position consists of one nucleotide
								else {
									start[0] = Integer.parseInt(position);
									stop[0] = Integer.parseInt(position);
									useGeneInCounting = false;
								}
							}
						}
						//If the positions do not contain the identifiers "join" and "complement"
						else {
							//Number of exons
							numberOfExons = 1;
							//Array containing the start position of the first exon
							start = new int[numberOfExons];
							//Array containing the end position of the first exon
							stop = new int[numberOfExons];
							
							//If both the beginning and the end are indicated
							if (position.contains("..")) {
								String[] arr = position.split("\\.\\.");
								start[0] = Integer.parseInt(arr[0]);
								stop[0] = Integer.parseInt(arr[1]);
							} 
							//If the position consists of one nucleotide
							else {
								start[0] = Integer.parseInt(position);
								stop[0] = Integer.parseInt(position);
								useGeneInCounting = false;
							}
						}
					}
					
					//Checking the distance between neighboring genes
					if (!organism.isEukaryota() && !(lastpozcomplement ^ complement) && (start[0] - lastpoz - parameters.getDistBOper() - 2 < 0)) {
						organism.delOperon();
					}
					
					//Position of the end of the gene
					lastpoz = stop[numberOfExons - 1];
					//Gene location on the strand
					lastpozcomplement = complement;
					
					//Recording an organism's gene
					Gene gene = new Gene(organism.getGenesNumber(), start, stop, organism.getOperonsNumber(), complement, useGeneInCounting);
					//Putting a new gene on the list
					organism.getGenes().add(gene);
					
					flagCDS = true;
					flagProduct = true;
					flagProtID = true;
					flagDBxref = true;
					//flagFBid = true;
					flagNote = true;
					flagLocus = true;
					flagPseudoGene = true;
				}
				
				//Name of the gene
				if (flagCDS && line.startsWith("/gene=")) {
					/*if (p4.matcher(line).matches()) {
						l.get(l.size()-1).putR(true);
						//Number of ribosomal genes
						this.ribgenes++;
					}*/
					
					//Determination of ribosomal genes by gene name
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMALG.matcher(line).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Number of ribosomal genes
							organism.addRibGene();
						}
					}
					
					//Determination of highly expressed genes by label
					if (parameters.isUseOtherGenesAsHighExpressed1() && line.contains(" highexpressed")) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Number of ribosomal genes
							organism.addRibGene();
						}
					}
					
					//The name of the gene
					organism.getGenes().get(organism.getGenes().size() - 1).setName(line.substring("/gene=".length()).replace("\"", ""));
					flagCDS = false;
				}
				
				//Definition of a pseudogene
				if (flagPseudoGene && line.toLowerCase().startsWith("/pseudo")) {
					organism.addPseudoGene();
					//The pseudogene mark
					organism.getGenes().get(organism.getGenes().size() - 1).setPseudo(true);
					flagPseudoGene = false;
				}
				
				//Gene index
				if (line.startsWith("/locus_tag=")) {
					if (flagLocus) {
						organism.getGenes().get(organism.getGenes().size() - 1).setLocusTag(line.substring("/locus_tag=".length()).replace("\"", ""));
						flagLocus = false;
					} else /*if (flagGene) */{
						this.geneNote.locus_tag = line.substring("/locus_tag=".length()).replace("\"", "");
					}
				}
				
				//Name of gene product
				if (flagProduct && line.startsWith("/product=")) {
					if (line.contains("COG")) {
						Matcher m2 = COGID.matcher(line);
						if (m2.find() && organism.getGenes().size() != 0) {
							organism.getGenes().get(organism.getGenes().size() - 1).setCog(m2.group());
						}
					}
					
					String line1 = line;
					/*while (!line.startsWith("                     /")){
						line1 += line.trim()+" ";
						line = reader.readLine();
					}*/
					
					//Recording the name of the gene product
					organism.getGenes().get(organism.getGenes().size() - 1).setProteinName(line1.substring("/product=".length()).replace("\"", ""));
					
					//Determination of ribosomal genes by name of gene product
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line1).matches() && !NONRIBOSOMAL.matcher(line1).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Number of ribosomal genes
							organism.addRibGene();
						}
					}
					flagProduct = false;
					flagCDS = false;
				}
				
				//protein_id=
				if (flagProtID && line.startsWith("/protein_id=")) {
					String prID = line.substring("/protein_id=".length()).replace("\"", "");
					organism.getGenes().get(organism.getGenes().size() - 1).setProteinID(prID);
					flagProtID = false;
				}
				
				//db_xref=\"GI:
				if (flagDBxref && line.startsWith("/db_xref=\"GI:")) {
					organism.getGenes().get(organism.getGenes().size() - 1).setGeneID(line.substring("/db_xref=".length()).replace("\"", ""));
					flagDBxref = false;
				}
				
				//Multi-line /note recording
				if (flagNoteMultistring) {
					this.geneNote.note += line.replace("\"", "");
					flagNoteMultistring = false;
				}
				
				//Recording of annotations
				if (line.startsWith("/note=")) { 
					if (flagNote){
						if (line.contains("COG")) {
							Matcher m2 = COGID.matcher(line);
							if (m2.find() && organism.getGenes().size() != 0) {
								organism.getGenes().get(organism.getGenes().size() - 1).setCog(m2.group());
							}
						}
						String line1 = line;
						/*while (!line.startsWith("                     /")){
							line1 += line.trim()+" ";
							line = reader.readLine();
						}*/
						
						organism.getGenes().get(organism.getGenes().size() - 1).setNote(line1.substring("/note=".length()).replace("\"", ""));
						
						//Determination of ribosomal genes by gene annotation
						if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line1).matches() && !NONRIBOSOMAL.matcher(line1).matches()) {
							if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
								organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
								//Number of ribosomal genes
								organism.addRibGene();
							}
						}
						flagNote = false;
					} else /*if (flagGene)*/ {
						//TODO: Parse COG in the gene feature using the same approach as above.
						//Definition of the beginning of the multiline /note
						String s = line.substring("/note=".length());
						this.geneNote.note = s.replace("\"", "");
						int countBrackets = s.length() - this.geneNote.note.length();
						if (countBrackets == 1) {
							flagNoteMultistring = true;
						}
						//Determination of ribosomal genes by gene annotation
						if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line).matches() && !NONRIBOSOMAL.matcher(line).matches()) {
							this.geneNote.isRib = true;
						}
					}
				}
				
				//FlyBase gene ID
				/*if (flagFBid && line.startsWith("/db_xref=\"FLYBASE:FB")) {
					organism.genes.get(organism.genes.size()-1).putFBid(line.substring(18, 29));
					flagFBid = false;
				}
				else {
					if (flagFBid && (organism.genes.size() > 0)) organism.genes.get(organism.genes.size()-1).putFBid("none");
				}*/
				
				//COG
				if (line.contains("COG")) {
					Matcher m2 = COGID.matcher(line);
					if (m2.find() && organism.getGenes().size() != 0) {
						organism.getGenes().get(organism.getGenes().size() - 1).setCog(m2.group());
					}
				}
				
				if (/*flagContext &&*/ line.startsWith("gene          ")) {
					flagPseudoGene = false;
					flagCDS = false;
					flagProduct = false;
					flagProtID = false;
					flagDBxref = false;
					//flagContext = false;
					flagNote = false;
					//2025.10.16 The last CDS for this gene has ended. Save the 
					//information collected from the gene feature that preceded the CDS.
					//Save the /note and ribosomal status.
					int num = organism.getGenes().size();
					//The code below does not apply to the first gene.
					if (num != 0) {
						Gene g = organism.getGenes().get(num - 1);
						if (null != this.geneNote.note && this.geneNote.locus_tag.equals(g.getLocusTag())) {
							g.setNote(this.geneNote.note);
							if (this.geneNote.isRib && !g.isRibosomal()) {
								g.setRibosomal(true);
								//Number of ribosomal genes
								organism.addRibGene();
							}
						}
						//2025.10.15 Reset the information about the previous gene.
						this.geneNote.reset();
					}
				}
				
				//Determination of tRNA quantity
				if (line.startsWith("tRNA          ")) {
					organism.addtRnaGene();
				}
			}
			
			//2025.10.16 Process the final CDS of the final gene. 
			//Save the /note and ribosomal status.
			int num = organism.getGenes().size();
			//Normally, the test below always passes.
			if (num != 0) {
				Gene g = organism.getGenes().get(num - 1);
				if (null != this.geneNote.note && this.geneNote.locus_tag.equals(g.getLocusTag())) {
					g.setNote(this.geneNote.note);
					if (this.geneNote.isRib && !g.isRibosomal()) {
						g.setRibosomal(true);
						//Number of ribosomal genes
						organism.addRibGene();
					}
				}
			}
			
			//============================================================================================
			//Read the complete nucleotide sequence from the record===============================
			//TODO========================================================================================
			
			//An array to record a complete genome
			char[] wholeGenome = new char[genomeSize];
			
			int pozstrGenome = 0;
			
			//Number of GC and CG pairs, GC-composition
			int gc = 0;
			
			//Genomic sequence reading
			while (!(line = reader.readLine()).startsWith("//")) {
				line = line.trim();
				
				String[] arr = line.split(" ");
				
				//Line for recording genomic sequence
				StringBuilder strGenome = new StringBuilder();
				for (int i = 1; i < arr.length; i++) {
					strGenome.append(arr[i]);
				}
				
				char[] ch = strGenome.toString().toCharArray();
				
				for (int i = 0; i < ch.length; i++) {
					if (ch[i] == 'g' || ch[i] == 'c') {
						gc++;
					}
				}
				
				System.arraycopy(ch, 0, wholeGenome, pozstrGenome, ch.length);
				pozstrGenome += ch.length;
			}
			
			//GC composition
			organism.addGC(gc);
			double gccom = Math.round(gc*10000.0/genomeSize)/100.0;
			
			int lFlang1;
			int rFlang1;
			
			//TODO Recording gene sequences in instances of a class
			int smax = Math.max(parameters.getSmax1(), parameters.getSmax2());
			int lmax = Math.max(parameters.getLmax1(), parameters.getLmax2());
			int minRFlang = smax+lmax;
			for (int i = 0; i < organism.getGenes().size(); i++) {
				if (organism.getGenes().get(i).isComplement()) {
					lFlang1 = parameters.getrFlang() < minRFlang ? minRFlang : parameters.getrFlang();
					rFlang1 = parameters.getlFlang();
				} else {
					lFlang1 = parameters.getlFlang();
					rFlang1 = parameters.getrFlang() < minRFlang ? minRFlang : parameters.getrFlang();
				}
				
				//Number of introns
				int numderOfIntrons = organism.getGenes().get(i).getStart().length;
				
				//The position of the gene onset
				int strt = organism.getGenes().get(i).getStart()[0];
				
				//Position of the end of the gene
				int stp = organism.getGenes().get(i).getStop()[numderOfIntrons - 1];
				int stp1 = stp + rFlang1 - genomeSize;
				int dlinaGenaBezFl = 0;
				int dlinaGena;
				
				//A string to record a gene sequence
				char[] strNA;
				
				//char strNA2[] = new char[2003];
				
				//Determination of gene length without flanks
				for (int j = 0; j < numderOfIntrons; j++) {
					dlinaGenaBezFl += organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
				}
				
				//If DNA is linear
				if (linear) {
					//If the left flank is longer than it is on the chromosome before the gene begins
					if (strt < lFlang1) {
						lFlang1 = strt - 1;
						organism.getGenes().get(i).getStart()[0] = 1;
						organism.getGenes().get(i).setFlagUse(false);
					} 
					//If the left flank is smaller than it is on the chromosome before the gene begins
					else {
						organism.getGenes().get(i).getStart()[0] -= lFlang1;
					}
					
					//If the right flank is longer than it is on the chromosome after the end of the gene
					if (stp1 > 0) {
						rFlang1 = genomeSize - stp;
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] = genomeSize;
						organism.getGenes().get(i).setFlagUse(false);
					} 
					//If the length of the right flank is smaller than it is on the chromosome after the end of the gene
					else {
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] += rFlang1;
					}
					
					//Length of gene with flanks
					dlinaGena = lFlang1 + dlinaGenaBezFl + rFlang1;
					
					//Array for recording gene sequences with flanks
					strNA = new char [dlinaGena];
					int pozEx = 0;
					
					//Reading a gene sequence from a complete genome
					for (int j = 0; j < numderOfIntrons; j++) {
						int dlinaEx = organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
						//Gene sequence string
						System.arraycopy(wholeGenome, organism.getGenes().get(i).getStart()[j] - 1, strNA, pozEx, dlinaEx);
						pozEx += dlinaEx;
					}
					
					/*if (organism.genes.get(i).getC()) {
						if (stp > 1003 && stp+1000 <= genomeSize) {
							System.arraycopy(wholeGenome, stp-1003, strNA2, 0, strNA2.length);
						}
						else {
							organism.genes.get(i).putFlagOut(false);
						}
					} else {
						if (strt > 1000 && strt+1003 <= genomeSize) {
							System.arraycopy(wholeGenome, strt-1000, strNA2, 0, strNA2.length);
						}
						else {
							organism.genes.get(i).putFlagOut(false);
						}
					}*/
					
				}
				//If the DNA is circular
				else {
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
					}
					//If the left flank is longer than it is on the chromosome before the gene begins
					else {
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
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] += rFlang1;
					}
					//If the right flank is longer than it is on the chromosome after the end of the gene
					else {
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] = genomeSize;
						rflng = true;
					}
					
					//Reading a gene sequence from a complete genome
					for (int j = 0; j < numderOfIntrons; j++) {
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
					strNA = Auxiliary_Class.reversComplement(strNA);
					//strNA2 = Auxiliary_Class.complement(strNA2);
				}
				
				//The position of the beginning of the gene without taking into account the flanks
				organism.getGenes().get(i).getStart()[0] = strt;
				//The position of the end of the gene without taking into account flanks
				organism.getGenes().get(i).getStop()[numderOfIntrons - 1] = stp;
				//Gene sequence recording
				organism.getGenes().get(i).setNA(strNA);
				//organism.genes.get(i).putNA2(strNA2);
				
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
