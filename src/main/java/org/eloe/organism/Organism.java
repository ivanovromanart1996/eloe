package org.eloe.organism;
import java.util.List;
import java.util.ArrayList;

public class Organism { //Class describing the organism
	
	private String name; //Name of the organism
	private int genomeSize = 0; //Genome size in nucleotides
	private int chromosomesNumber = 0; //Number of chromosomes in the organism
	private int genesNumber = 0; //Number of genes in the genome
	private int ribGenesNumber = 0; //Number of ribosomal genes in the genome
	private int pseudoGenesNumber = 0; //Number of pseudogenes in the genome
	private int shortGenesNumber = 0; //Number of short genes in the genome
	private int tRnaGenesNumber = 0; //Number of copies of tRNA genes in the genome
	private int operonsNumber = 0; //Number of operons in the genome
	private int gcNumber = 0; //Number of GCs in the genome
	private double gcComposition; //GC composition
	private boolean archaea = false;
	private boolean eukaryota = false; //Eukaryotic organism or not
	private int index = 1; //Organism index
	private List<Gene> genes = new ArrayList<Gene>(100); //List of genes in the organism
	
	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		str.append("Organism: "+this.name+"\n");
		str.append("Genome size: "+this.genomeSize/1000.0+" Kbp\n");
		str.append("Chromosomes number: "+this.chromosomesNumber+"\n");
		str.append("Genes number: "+this.genesNumber+"\n");
		str.append("Ribosomal genes number: "+this.ribGenesNumber+"\n");
		str.append("Pseudo genes number: "+this.pseudoGenesNumber+"\n");
		str.append("Short genes number: "+this.shortGenesNumber+"\n");
		str.append("tRNA genes number: "+this.tRnaGenesNumber+"\n");
		str.append("GC-composition: "+this.gcComposition+"%");
		return str.toString();
	}
	
	public void addChromosome(Organism chr) {
		this.name = chr.name;
		this.genomeSize += chr.genomeSize;
		this.chromosomesNumber++;
		this.genesNumber += chr.getGenesNumber();
		this.ribGenesNumber += chr.getRibGenesNumber();
		this.pseudoGenesNumber += chr.getPseudoGenesNumber();
		this.shortGenesNumber += chr.getShortGenesNumber();
		this.tRnaGenesNumber += chr.gettRnaGenesNumber();
		this.operonsNumber += chr.getOperonsNumber();
		this.gcNumber += chr.gcNumber;
		this.gcComposition = this.gcNumber*100.0/this.genomeSize;
		this.archaea = chr.isArchaea();
		this.eukaryota = chr.isEukaryota();
		this.genes.addAll(chr.getGenes());
	}
		
	//Marking of operons in a prokaryotic organism
	public void markCistrons() {
		if (!this.eukaryota) {
			//The first gene is the first cistron.
			this.genes.get(0).setFirstCistron(true);
			for (int i = 1; i < this.genes.size()-1; i++) {
				if (this.genes.get(i).getnOper() == this.genes.get(i-1).getnOper()) {
					if (this.genes.get(i).getnOper() == this.genes.get(i+1).getnOper()) {
						this.genes.get(i).setMiddleCistron(true);
					}
					else {
						this.genes.get(i).setLastCistron(true);
						this.genes.get(i+1).setFirstCistron(true);
						i++;
					}
				}
				else {
					this.genes.get(i).setFirstCistron(true);
					this.genes.get(i-1).setLastCistron(true);
				}
			}
			
			if (this.genes.get(this.genes.size()-1).getnOper() != this.genes.get(this.genes.size()-2).getnOper()) {
				this.genes.get(this.genes.size()-1).setFirstCistron(true);
				this.genes.get(this.genes.size()-2).setLastCistron(true);
			}
			//The last gene is the last cistron.
			this.genes.get(this.genes.size()-1).setLastCistron(true);
		}
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getGenomeSize() {
		return genomeSize;
	}

	public void setGenomeSize(int genomeSize) {
		this.genomeSize = genomeSize;
	}

	public int getChromosomesNumber() {
		return chromosomesNumber;
	}

	public void setChromosomesNumber(int chromosomesNumber) {
		this.chromosomesNumber = chromosomesNumber;
	}

	public int getGenesNumber() {
		return genesNumber;
	}

	public void setGenesNumber(int genesNumber) {
		this.genesNumber = genesNumber;
	}
	
	public void addGene() {
		this.genesNumber++;
	}
	
	public void delGene() {
		this.genesNumber--;
	}

	public int getRibGenesNumber() {
		return ribGenesNumber;
	}

	public void setRibGenesNumber(int ribGenesNumber) {
		this.ribGenesNumber = ribGenesNumber;
	}
	
	public void addRibGene() {
		this.ribGenesNumber++;
	}
	
	public void delRibGene() {
		this.ribGenesNumber--;
	}
	
	public int getPseudoGenesNumber() {
		return pseudoGenesNumber;
	}

	public void setPseudoGenesNumber(int pseudoGenesNumber) {
		this.pseudoGenesNumber = pseudoGenesNumber;
	}
	
	public void addPseudoGene() {
		this.pseudoGenesNumber++;
	}
	
	public void delPseudoGene() {
		this.pseudoGenesNumber--;
	}

	public int getShortGenesNumber() {
		return shortGenesNumber;
	}

	public void setShortGenesNumber(int shortGenesNumber) {
		this.shortGenesNumber = shortGenesNumber;
	}
	
	public void addShortGene() {
		this.shortGenesNumber++;
	}
	
	public void delShortGene() {
		this.shortGenesNumber--;
	}
	
	public int gettRnaGenesNumber() {
		return tRnaGenesNumber;
	}

	public void settRnaGenesNumber(int tRnaGenesNumber) {
		this.tRnaGenesNumber = tRnaGenesNumber;
	}
	
	public void addtRnaGene() {
		this.tRnaGenesNumber++;
	}
	
	public void deltRnaGene() {
		this.tRnaGenesNumber--;
	}
	
	public int getOperonsNumber() {
		return operonsNumber;
	}

	public void setOperonsNumber(int operonsNumber) {
		this.operonsNumber = operonsNumber;
	}
	
	public void addOperon() {
		this.operonsNumber++;
	}
	
	public void delOperon() {
		this.operonsNumber--;
	}

	public int getGcNumber() {
		return gcNumber;
	}

	public void setGcNumber(int gcNumber) {
		this.gcNumber = gcNumber;
	}
	
	public void addGC(int gc) {
		this.gcNumber += gc;
	}
	
	public double getGcComposition() {
		return gcComposition;
	}

	public void setGcComposition(double gcComposition) {
		this.gcComposition = gcComposition;
	}

	public boolean isArchaea() {
		return archaea;
	}

	public void setArchaea(boolean archaea) {
		this.archaea = archaea;
	}

	public boolean isEukaryota() {
		return eukaryota;
	}

	public void setEukaryota(boolean eukaryota) {
		this.eukaryota = eukaryota;
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public List<Gene> getGenes() {
		return genes;
	}

	public void setGenes(List<Gene> genes) {
		this.genes = genes;
	}
	
}
