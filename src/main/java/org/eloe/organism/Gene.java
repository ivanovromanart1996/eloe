package org.eloe.organism;

import org.eloe.parameters.PGeneToString;

// Класс нуклеотидной последовательности с параметрами: длина, последовательность, позиция начала, 
// позиция конца, расположение в прямой или комплементарной цепи, принадлежность к рибосомальным генам.
public class Gene implements IGeneWithFlanks {
	private int number; //Порядковый номер гена в геноме
	private String name = "-"; //Название гена
	private String locusTag; //Индекс гена
	private String geneID; //Индекс гена
	private String cog = "-"; //COG индекс
	private String fbID; //id гена в FLYBASE
	private String note; //Аннотация гена
	private String proteinName; //Продукт гена
	private String proteinID; //Индекс белка
	private int length; //Длина гена без флангов
	private int lengthWithFlanks; //Длина гена с флангами
	private int leftFlankLength = 0; //Длина левого фланга
	private int rightFlankLength = 0; //Длина правого фланга
	private int codonsNumber; //Число кодонов без флангов
	private int codonsWithFlanksNumber; //Число кодонов с флангами
	private int[] start; //Начало гена
	private int[] stop; //Конец гена
	private boolean complement; //Комплемент
	private boolean ribosomal = false; //Рибосомный или нет
	private boolean pseudo = false;
	private boolean highExpressed = false; //Высоко экспрессируемый или нет
	private double gcComposition; //GC-состав гена
	private int[] codonComposition0 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
										0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
										0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
										0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	private int[] codonComposition = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			   							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			   							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
			   							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; //Массив для записи частот кодонов
	private double tai = 0; //Индекс по частотам кодонов
	private double lci1 = 0; //Индекс локальной комплементарности
	private double lci2 = 0; //Индекс локальной комплементарности
	private double pi1 = 0; //Вероятность в формуле для Te(i) без учета энергии
	private double pi2 = 0; //Вероятность в формуле для Te(i) с учетом энергии
	private double tei = 0; //Индекс по локальной комплементарности
	private double eei = 0; //Индекс эффективности элонгации
	private double[] eeiMass = {0, 0, 0, 0, 0}; //Все индексы eei1-eei5
	private int repeats = 0; //Количество повторов в гене
	private double[] lciij; //LCI индивидуального нуклеотида
	private boolean flagUse = true; //Флаг использования
	private boolean flagOut = true; //Флаг вывода в файл
	private int nOper; //Номер оперона
	private boolean firstCistron = false;
	private boolean middleCistron = false;
	private boolean lastCistron = false;
	private char[] naWithFlanks; //Последовательность нуклеотидов
	
	//Конструктор
	public Gene(int genes, int[] begining, int[] end, int nOper, boolean complement, boolean useGeneInCounting) { //Конструктор
		this.number = genes;
		this.start = begining;
		this.stop = end;
		this.nOper = nOper;
		this.complement = complement;
		this.flagUse = useGeneInCounting;
	}
	
	//TODO
	//Вывод данных
	public String toString(Organism organism, PGeneToString parameters) {
		int rightFlang1 = (this.rightFlankLength <= parameters.getrFlangTrue()) ? this.rightFlankLength : parameters.getrFlangTrue();
		int leftFlang1 = (this.leftFlankLength <= parameters.getlFlangTrue()) ? this.leftFlankLength : parameters.getlFlangTrue();
		
		StringBuilder a = new StringBuilder(">");
		a.append(this.number).append("\t");
//		if (this.rightFl != this.rFlangTrue || this.leftFl != this.lFlangTrue) { //Пометка "*", если другая длина фланга
//			a+="*\t";
//		}
//		if (this.ribosomal) a.append("rib\t"); //Пометка "rib"
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
	//TODO Получение данных
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
	
	//Записывает длину последовательности и число кодонов
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

	//Записывает GC-состав гена в процентах
	public void putGC_compoz() {
		double gc_compos1 = 0;
		//Количество gс в триплете
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
	
	//Возвращает строку с частотами кодонов через ";"
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
	
	//Возвращает строку с частотами кодонов через ";"
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
	
	//Записывает EEI1
	public void putEEI1(double eei1) {
		this.eeiMass[0] = ((eei1 != 0) ? (1/eei1) : 10);
	}
	
	//Записывает EEI2
	public void putEEI2(double eei2) {
		this.eeiMass[1] = ((eei2 != 0) ? (1/eei2) : 10);
	}
	
	//Записывает EEI3
	public void putEEI3(double eei3) {
		this.eeiMass[2] = ((eei3 != 0) ? (1/eei3) : 10);
	}
	
	//Записывает EEI4
	public void putEEI4(double eei4) {
		this.eeiMass[3] = ((eei4 != 0) ? (2/eei4) : 10);
	}
	
	//Записывает EEI5
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
	
	//Возвращает последовательность
	public String getNAString() {
		return String.valueOf(this.getNA());
	}
	
	//Возвращает последовательность как требует пользователь
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
	
	//Вывод данных
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
	//Вывод данных 2
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
