package org.eloe.organism;

import org.eloe.parameters.PGeneToString;

// Class representing a nucleotide sequence, with length, sequence, start position, 
// end position, forward or complementary strand, and ribosomal status.
public class Gene implements IGeneWithFlanks {
	private int number; //Gene sequence number in the genome
	private String name = "-"; //Name of the gene
	private String locusTag; //Gene index
	private String geneID; //Gene index
	private String cog = "-"; //COG index
	private String fbID; //FlyBase gene ID
	private String note; //Gene annotation
	private String proteinName; //Gene product
	private String proteinID; //protein index
	private int length; //Gene length without flanks
	private int lengthWithFlanks; //Length of gene with flanks
	private int leftFlankLength = 0; //Length of left flank
	private int rightFlankLength = 0; //Length of right flank
	private int codonsNumber; //Number of codons without flanks
	private int codonsWithFlanksNumber; //Number of codons with flanks
	private int[] start; //Gene start
	private int[] stop; //End of gene
	private boolean complement; //Complement
	private boolean ribosomal = false; //Whether the gene is ribosomal
	private boolean pseudo = false;
	private boolean highExpressed = false; //Whether the gene is highly expressed
	private double gcComposition; //GC composition
	private int[] codonComposition0 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
										0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
										0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
										0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	private int[] codonComposition = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			   							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			   							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			   							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; //Array for recording codon frequencies
	private double tai = 0; //Codon frequency index
	private double lci1 = 0; //Local Complementarity Index
	private double lci2 = 0; //Local Complementarity Index
	private double pi1 = 0; //Probability in the formula for Te(i) excluding energy
	private double pi2 = 0; //Probability in the formula for Te(i) taking into account energy
	private double tei = 0; //Local Complementarity Index
	private double eei = 0; //Elongation efficiency index
	private double[] eeiMass = {0, 0, 0, 0, 0}; //All eei1-eei5 indices
	private int repeats = 0; //Number of repetitions in the gene
	private double[] lciij; //LCI individual nucleotide
	private boolean flagUse = true; //Flag of use
	private boolean flagOut = true; //File output flag
	private int nOper; //Operon number
	private boolean firstCistron = false;
	private boolean middleCistron = false;
	private boolean lastCistron = false;
	private char[] naWithFlanks; //nucleotide sequence
	
	//Constructor
	public Gene(int genes, int[] begining, int[] end, int nOper, boolean complement, boolean useGeneInCounting) { //Constructor
		this.number = genes;
		this.start = begining;
		this.stop = end;
		this.nOper = nOper;
		this.complement = complement;
		this.flagUse = useGeneInCounting;
	}
	
	//TODO
	//Data output
	public String toString(Organism organism, PGeneToString parameters) {
		int rightFlang1 = (this.rightFlankLength <= parameters.getrFlangTrue()) ? this.rightFlankLength : parameters.getrFlangTrue();
		int leftFlang1 = (this.leftFlankLength <= parameters.getlFlangTrue()) ? this.leftFlankLength : parameters.getlFlangTrue();
		
		StringBuilder a = new StringBuilder(">");
		a.append(this.number).append("\t");
//		if (this.rightFl != this.rFlangTrue || this.leftFl != this.lFlangTrue) { //Mark with "*" if the flank length differs
//			a+="*\t";
//		}
//		if (this.ribosomal) a.append("rib\t"); //Mark with "rib"
		if (!organism.isEukaryota()) {
			a.append("N operona: ").append(this.nOper).append("\t");
			if (this.firstCistron && this.lastCistron) a.append("first-last\t");
			else {
				if (this.firstCistron) a.append("first\t");
				else if (this.lastCistron) a.append("last\t");
					else a.append("middle\t");
			}
		}
		a.append(this.locusTag).append("\t");
		a.append("gene_poz=(").append(leftFlang1+1).append("-").append(this.length+leftFlang1).append(")\t");
		a.append(this.name).append("\t");
		a.append(this.proteinName).append("\t");
		a.append(this.note).append("\t");
		a.append("/position=");
		if (this.complement) {
			a.append(this.start[0]-rightFlang1).append(" - ").append(this.stop[this.stop.length-1]+leftFlang1);
			a.append(" (").append(this.start[0]).append(" - ").append(this.stop[this.stop.length-1]).append(")\t");
			a.append("/lenght=").append(this.lengthWithFlanks-this.rightFlankLength+rightFlang1);
			a.append(" (").append(this.length).append(")").append("\t");
			a.append("/codons=").append((this.lengthWithFlanks-this.rightFlankLength+rightFlang1)/3);
			a.append(" (").append(this.codonsNumber).append(")\t");
			a.append("/complement=true\t");
		}
		else {
			a.append(this.start[0]-leftFlang1).append(" - ").append(this.stop[this.stop.length-1]+rightFlang1);
			a.append(" (").append(this.start[0]).append(" - ").append(this.stop[this.stop.length-1]).append(")\t");
			a.append("/lenght=").append(this.lengthWithFlanks-this.rightFlankLength+rightFlang1);
			a.append(" (").append(this.length).append(")");
			a.append("\t/codons=").append((this.lengthWithFlanks-this.rightFlankLength+rightFlang1)/3);
			a.append(" (").append(this.codonsNumber).append(")\t");
			a.append("/complement=false\t");
		}
		a.append("GC-compoz=").append(this.gcComposition).append("%\t");
		a.append("/LCI1=").append(this.lci1).append("\t");
		a.append("/LCI2=").append(this.lci2).append("\t");
		a.append("p_i1=").append(this.pi1).append("\t");
		a.append("p_i2=").append(this.pi2).append("\t");
		a.append("/Tai=").append(this.tai).append("\t");
		a.append("/Tei=").append(this.tei).append("\t");
		a.append("/EEI=");
//		for (int i = 0; i < this.lciij.length; i++) {
//			a.append(lciij[i]).append("\t");
//		}
		return a.toString();
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//TODO Obtaining data
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	public int getNumber() {
		return number;
	}

	public void setNumber(int number) {
		this.number = number;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getLocusTag() {
		return locusTag;
	}

	public void setLocusTag(String locusTag) {
		this.locusTag = locusTag;
	}

	public String getGeneID() {
		return geneID;
	}

	public void setGeneID(String geneID) {
		this.geneID = geneID;
	}

	public String getCog() {
		return cog;
	}

	public void setCog(String cog) {
		this.cog = cog;
	}

	public String getFbID() {
		return fbID;
	}

	public void setFbID(String fbID) {
		this.fbID = fbID;
	}

	public String getNote() {
		return note;
	}

	public void setNote(String note) {
		this.note = note;
	}

	public String getProteinName() {
		return proteinName;
	}

	public void setProteinName(String proteinName) {
		this.proteinName = proteinName;
	}

	public String getProteinID() {
		return proteinID;
	}

	public void setProteinID(String proteinID) {
		this.proteinID = proteinID;
	}

	public int getLength() {
		return length;
	}

	public void setLenght(int lenght) {
		this.length = lenght;
	}
	
	//Records sequence length and number of codons
	public void putLenght(int a, int b, int c, int d) {
		this.leftFlankLength = a;
		this.length = b;
		this.codonsNumber = b/3;
		this.rightFlankLength = c;
		this.lengthWithFlanks = d;
		this.codonsWithFlanksNumber = d/3;
	}

	public int getLenghtWithFlanks() {
		return lengthWithFlanks;
	}

	public void setLenghtWithFlanks(int lenghtWithFlanks) {
		this.lengthWithFlanks = lenghtWithFlanks;
	}

	public int getLeftFlankLength() {
		return leftFlankLength;
	}

	public void setLeftFlankLength(int leftFlankLength) {
		this.leftFlankLength = leftFlankLength;
	}

	public int getRightFlankLength() {
		return rightFlankLength;
	}

	public void setRightFlankLength(int rightFlankLength) {
		this.rightFlankLength = rightFlankLength;
	}

	public int getCodonsNumber() {
		return codonsNumber;
	}

	public void setCodonsNumber(int codonsNumber) {
		this.codonsNumber = codonsNumber;
	}

	public int getCodonsWithFlanksNumber() {
		return codonsWithFlanksNumber;
	}

	public void setCodonsWithFlanksNumber(int codonsWithFlanksNumber) {
		this.codonsWithFlanksNumber = codonsWithFlanksNumber;
	}

	public int[] getStart() {
		return start;
	}

	public void setStart(int[] start) {
		this.start = start;
	}

	public int[] getStop() {
		return stop;
	}

	public void setStop(int[] stop) {
		this.stop = stop;
	}

	public boolean isComplement() {
		return complement;
	}

	public void setComplement(boolean complement) {
		this.complement = complement;
	}

	public boolean isRibosomal() {
		return ribosomal;
	}

	public void setRibosomal(boolean ribosomal) {
		this.ribosomal = ribosomal;
	}

	public boolean isPseudo() {
		return pseudo;
	}

	public void setPseudo(boolean pseudo) {
		this.pseudo = pseudo;
	}

	public boolean isHighExpressed() {
		return highExpressed;
	}

	public void setHighExpressed(boolean highExpressed) {
		this.highExpressed = highExpressed;
	}

	public double getGcComposition() {
		return gcComposition;
	}

	public void setGcComposition(double gcComposition) {
		this.gcComposition = gcComposition;
	}

	//Records the GC composition of the gene as a percentage
	public void putGC_compoz() {
		double gc_compos1 = 0;
		//Number of G bases in the triplet
		int gc[] = {2, 3, 3, 2, 1, 2, 1, 2, 2, 1, 0, 1, 1, 2, 2, 1,
					2, 1, 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 2, 3,
					3, 2, 1, 2, 2, 1, 0, 1, 1, 0, 1, 2, 2, 1, 1, 2,
					2, 1, 1, 0, 2, 1, 1, 0, 1, 0, 0, 1, 2, 1, 0, 1};
		
		for (int i = 0; i < this.codonComposition0.length; i++) {
			gc_compos1 += this.codonComposition0[i]*gc[i];
		}
		
		this.gcComposition = gc_compos1/this.length*100.0;
	}

	public int[] getCodonComposition0() {
		return codonComposition0;
	}
	
	//Returns a semicolon-separated string of codon frequencies.
	public String getCodonComposition0String() {
		StringBuilder freq = new StringBuilder();
		for (int i = 0; i < codonComposition0.length; i++) {
			freq.append(codonComposition0[i]+";\t");
		}
		return freq.toString();
	}

	public void setCodonComposition0(int[] codonComposition0) {
		this.codonComposition0 = codonComposition0;
	}

	public int[] getCodonComposition() {
		return codonComposition;
	}
	
	//Returns a semicolon-separated string of codon frequencies.
	public String getCodonCompositionString() {
		StringBuilder freq = new StringBuilder();
		for (int i = 0; i < codonComposition.length; i++) {
			freq.append(codonComposition[i]+";\t");
		}
		return freq.toString();
	}

	public void setCodonComposition(int[] codonComposition) {
		this.codonComposition = codonComposition;
	}

	public double getTai() {
		return tai;
	}

	public void setTai(double tai) {
		this.tai = tai;
	}

	public double getLci1() {
		return lci1;
	}

	public void setLci1(double lci1) {
		this.lci1 = lci1;
	}

	public double getLci2() {
		return lci2;
	}

	public void setLci2(double lci2) {
		this.lci2 = lci2;
	}

	public double getPi1() {
		return pi1;
	}

	public void setPi1(double pi1) {
		this.pi1 = pi1;
	}

	public double getPi2() {
		return pi2;
	}

	public void setPi2(double pi2) {
		this.pi2 = pi2;
	}

	public double getTei() {
		return tei;
	}

	public void setTei(double tei) {
		this.tei = tei;
	}

	public double getEei() {
		return eei;
	}

	public void setEei(double eei) {
		this.eei = eei;
	}

	public double[] getEeiMass() {
		return eeiMass;
	}

	public void setEeiMass(double[] eeiMass) {
		this.eeiMass = eeiMass;
	}
	
	//Recording EEI1.
	public void putEEI1(double eei1) {
		this.eeiMass[0] = ((eei1 != 0) ? (1/eei1) : 10);
	}
	
	//Recording EEI2.
	public void putEEI2(double eei2) {
		this.eeiMass[1] = ((eei2 != 0) ? (1/eei2) : 10);
	}
	
	//Recording EEI3.
	public void putEEI3(double eei3) {
		this.eeiMass[2] = ((eei3 != 0) ? (1/eei3) : 10);
	}
	
	//Recording EEI4.
	public void putEEI4(double eei4) {
		this.eeiMass[3] = ((eei4 != 0) ? (2/eei4) : 10);
	}
	
	//Recording EEI5.
	public void putEEI5(double eei5) {
		this.eeiMass[4] = ((eei5 != 0) ? (2/eei5) : 10);
	}

	public int getRepeats() {
		return repeats;
	}

	public void setRepeats(int repeats) {
		this.repeats = repeats;
	}

	public double[] getLciij() {
		return lciij;
	}

	public void setLciij(double[] lciij) {
		this.lciij = lciij;
	}

	public boolean isFlagUse() {
		return flagUse;
	}

	public void setFlagUse(boolean flagUse) {
		this.flagUse = flagUse;
	}

	public boolean isFlagOut() {
		return flagOut;
	}

	public void setFlagOut(boolean flagOut) {
		this.flagOut = flagOut;
	}

	public int getnOper() {
		return nOper;
	}

	public void setnOper(int nOper) {
		this.nOper = nOper;
	}

	public boolean isFirstCistron() {
		return firstCistron;
	}

	public void setFirstCistron(boolean firstCistron) {
		this.firstCistron = firstCistron;
	}

	public boolean isMiddleCistron() {
		return middleCistron;
	}

	public void setMiddleCistron(boolean middleCistron) {
		this.middleCistron = middleCistron;
	}

	public boolean isLastCistron() {
		return lastCistron;
	}

	public void setLastCistron(boolean lastCistron) {
		this.lastCistron = lastCistron;
	}
	
	public char[] getNA() {
		char dest[] = new char[length];
		System.arraycopy(naWithFlanks, leftFlankLength, dest, 0, dest.length);
		return dest;
	}

	public char[] getNaWithFlanks() {
		return naWithFlanks;
	}
	
	//Returns the sequence
	public String getNAString() {
		return String.valueOf(this.getNA());
	}
	
	//Returns the sequence as required by the user
	public char[] getUserNA(PGeneToString parameters) {
		if (parameters.getrFlangTrue() < this.rightFlankLength) {
			char dest[] = new char[this.leftFlankLength+this.length+parameters.getrFlangTrue()];
			System.arraycopy(this.naWithFlanks, 0, dest, 0, dest.length);
			return dest;
		}
		return this.naWithFlanks;
	}

	public void setNA(char[] na) {
		this.naWithFlanks = na;
	}
	
	//Data output
	/*public String toString() { 
		int rightFlang1 = (this.rightFl <= this.rFlangTrue) ? this.rightFl : this.rFlangTrue;
		int leftFlang1 = (this.leftFl <= this.lFlangTrue) ? this.leftFl : this.lFlangTrue;
		StringBuilder stroka = new StringBuilder();
		stroka.append(this.number);
		stroka.append("\t");
		stroka.append(this.locus_tag);
		stroka.append("\t");
		stroka.append(this.lFlangTrue+1);
		stroka.append("-");
		stroka.append(leftFlang1+this.lenghtBezFl);
		stroka.append("\t");
		stroka.append(this.geneName);
		stroka.append("\t");
		stroka.append(this.prot);
		stroka.append("\t");
		stroka.append(this.note);
		stroka.append("\t");
		if (this.complement) {
			stroka.append(this.begining[0]-rightFlang1);
			stroka.append(" - ");
			stroka.append(this.end[this.end.length-1]+leftFlang1);
			stroka.append(" (");
			stroka.append(this.begining[0]);
			stroka.append(" - ");
			stroka.append(this.end[this.end.length-1]);
			stroka.append(")\t");
			stroka.append(this.lenght-this.rightFl+rightFlang1);
			stroka.append(" (");
			stroka.append(this.lenghtBezFl);
			stroka.append(")\t");
			stroka.append((this.lenght-this.rightFl+rightFlang1)/3);
			stroka.append(" (");
			stroka.append(this.codonsBezFl);
			stroka.append(")\ttrue");
		}
		else {
			stroka.append(this.begining[0]-leftFlang1);
			stroka.append(" - ");
			stroka.append(this.end[this.end.length-1]+rightFlang1);
			stroka.append(" (");
			stroka.append(this.begining[0]);
			stroka.append(" - ");
			stroka.append(this.end[this.end.length-1]);
			stroka.append(")\t");
			stroka.append(this.lenght-this.rightFl+rightFlang1);
			stroka.append(" (");
			stroka.append(this.lenghtBezFl);
			stroka.append(")\t");
			stroka.append((this.lenght-this.rightFl+rightFlang1)/3);
			stroka.append(" (");
			stroka.append(this.codonsBezFl);
			stroka.append(")\tfalse");
		}
		stroka.append("\t");
		stroka.append(this.gc_compoz);
		stroka.append("\t");
		return stroka.toString();
	}*/
	
	//"Number\tlocus_tag\tFBid\tgene_poz\t/position=\t/complement=\t/EEI=\tgeneNA"
	//Data output 2
	/*public String toString() {
		int rightFlang1 = (this.rightFl <= this.rFlangTrue) ? this.rightFl : this.rFlangTrue;
		int leftFlang1 = (this.leftFl <= this.lFlangTrue) ? this.leftFl : this.lFlangTrue;
		String a = this.number+"\t"+this.locus_tag+"\t"+this.FBid+"\t"+(this.lFlangTrue+1)+
		"\t"+(leftFlang1+this.lenghtBezFl)+"\t"+this.begining[0]+"\t"+this.end[this.end.length-1]+"\t";
		if (this.complement) {
			a += "true\t";
		}
		else {
			a += "false\t";
		}
		return a;
	}*/
	
	
	
	
	
}
