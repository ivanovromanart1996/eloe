package org.eloe.extractor;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.regex.Pattern;

import org.eloe.organism.Organism;
import org.eloe.parameters.*;

// Extracts the genome sequence from the file specified in "file.txt" 
// and writes it to "filename - WholeGenom.fasta". It also extracts 
// gene sequences with flanks of a specified size and writes their position,  
// length, strand and ribosomal status  
// to "filename - GenesAndFlanks.fasta". It also writes "filename - Freq.txt"
// with codon frequencies in the genes.

public abstract class Extractor {
	//Pattern for organism names
	protected static final String ORGNAME = Pattern.compile("[:;|<>\\*\\?\\/\\\\]").pattern();
	//Pattern for organism names
	protected static final String ORGNAME_ = Pattern.compile("( |\\.|_)+").pattern();
	//Pattern for gene positions
	protected static final String BRACKETS = Pattern.compile("[\\(\\)<>]").pattern();
	//Regular expression to search for ribosomal proteins
	protected static final Pattern RIBOSOMAlP = Pattern.compile("(?i)^.*ribosomal\\s*(subunit)?\\s*protei.?(\\w*).*");
	//Regular expression for weeding out nonribosomal proteins
	protected static final Pattern NONRIBOSOMAL = Pattern.compile("(?i).*(mitochondrial|modification|methyltransferase|methylthiotransferase|acetyltransferase|methylation|modulation).*");
	//Regular expression to search for ribosomal protein genes
	protected static final Pattern RIBOSOMALG = Pattern.compile("(?i).*/gene=\"rp(l|s)\\w+\".*");
	//Pattern for COG identifier
	protected static final Pattern COGID = Pattern.compile("COG\\d+");
	//Regular expression to search for specific genes
	//protected static final Pattern p4 = Pattern.compile("(?i).*/gene=\"(RPL41A|RPS27A|RPL41B|RPL22B|RPL15B|RPS29B|RPS29A|RPS27B|RPS22B|RPL37B|RPL37A|RPS28A|RPS14B|RPL43B|RPL27B|RPS9A|RPS28B|RPS21A|RPS11A|RPS21B|RPL20A|RPS7B|RPS11B|RPL29|RPL20B|RPL42A|RPL6B|RPS18B|RPL23A|RPS26B|RPS25B|RPL36A|RPL28|RPL21A|RPS4A|RPS0A|RPL40A|RPL33A|RPL18B|RPS4B|RPL12A|RPS18A|RPL21B|RPL11B|RPL14A|RPL9A|RPL40B|PL33B|RPL31B|RPL9B|RPL14B|RPL17B|RPL42B|RPL26A|RPL32|RPS26A|RPL6A|RPL15A|RPS1A|RPL3|RPL43A|RPL34A|RPS31|RPL19B|RPL23B|RPL19A|RPS8B|RPL34B|RPL13A|RPS13|RPS8A|RPL26B|RPL10|RPS30A|RPL16A|RPL4A|RPS15|RPL2A|RPL2B|RPL1A|RPL4B|RPL1B|RPL39|RPS10B|RPS0B|RPS10A|RPL27A|RPS17B|RPS22A|RPL16B|RPL31A|RPS25A|RPL11A|RPL7B|RPL7A|RPS16A|RPS20|RPS2|RPS16B|RPS17A|RPS23A|RPS19B|RPL5|RPS1B|RPL22A|RPL18A|RPL13B|RPS7A|RPL36B|RPS14A|RPS30B|RPS19A|RPL17A|RPL25|RPS23B|RPL24B|RPL35B|RPS3|RPS24B|RPS5|RPL12B|PL24A|RPL8A|RPS9B|RPL8B|RPS24A|RPS6A|RPS6B|RPL35A|RPL38|RPS12|RPL30)\".*");
	
	//File name with data
	protected File orgFile;
	//Parameters of extraction
	protected PExtraction parameters;
	//Folder for results
	protected File resultFolder;
	//Record name
	protected String card;
	
	protected Extractor(File orgFile, PExtraction parameters, File resultFolder) {
		this.orgFile = orgFile;
		this.parameters = parameters;
		this.resultFolder = resultFolder;
	}
	
	public abstract Organism extract();
	
	//Returns the record name
	public final String getCard() {
		return this.card;
	}
	
	//Creating a file with complete genomic sequence
	protected void printWholeGenome(Organism organism, char[] wholeGenome, double gccom) {
		if (parameters.isMakeFileWholeGenome()) {
			File fileOut = new File(resultFolder, organism.getName() + "_whole_genome_" + organism.getChromosomesNumber() + ".fasta");
			PrintStream out = null;
			try {
				out = new PrintStream(fileOut);
				//Write organism information.
				out.println(">" + organism.getName() + "\t/Number of genes: " + organism.getGenesNumber() + "\t/Number of ribosomal genes: " 
						+ organism.getRibGenesNumber() + "\t/Number of tRNA: " + organism.gettRnaGenesNumber() + "\t/G+C composition (%): " + gccom);
				//Genome recording
				out.println(wholeGenome);
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (out != null) {
					out.close();
				}
			}
		}
	}
}
