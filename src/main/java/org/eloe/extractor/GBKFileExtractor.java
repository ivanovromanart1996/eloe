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
	
	//2025.10.15 контейнер для хранения /locus_tag и /note гена (gene), а также информацию о его рибосомальности
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
		
		//Считывание файла с геномом
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(orgFile));
			String line = reader.readLine();
			
			if (!line.startsWith("LOCUS ")) {
				LOG.error("ERROR! Wrong file format or file is empty: " + orgFile.getName() + "! Right formats: *.gbk, *.gb, *.emb");
				return null;
			}
			
			//Линейная ДНК или кольцевая
			boolean linear = true;
			if (line.contains(" circular ")) {
				linear = false;
			}
			
			String[] arr3 = line.split("\\s+");
			//Название карточки
			this.card = arr3[1];
			//Размер генома
			int genomeSize = Integer.parseInt(arr3[2]);
			organism.setGenomeSize(genomeSize);

			line = reader.readLine();
			//If file contains plasmid NA
			if (line.contains("plasmid")) {
				LOG.info("File with plasmid - " + orgFile.getName());
				return null;
			}
			
			//Определение названия организма
			boolean flagOrg = true;
			do {
				line = reader.readLine().trim();
				//Считывание названия организма
				if (line.startsWith("ORGANISM")) {
					String organismName = line.replace("ORGANISM  ", "");
					//Форматирование имени организма, чтобы его можно было использовать в названии файла
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
	        //Считывание из карточки информации о генах===================================================
	        //TODO========================================================================================
			
			//Флаг для определения гена (CDS)
			boolean flagCDS = false;
			//Флаг для определения псевдогена
			boolean flagPseudoGene = false;
			//Флаг для определения продукта гена
			boolean flagProduct = false;
			//Флаг для определения индекса белка
			boolean flagProtID = false;
			//Флаг для определения индекса гена
			boolean flagDBxref = false;
			//Флаг для определения id гена в FLYBASE
			//boolean flagFBid = false;
			//Флаг для нотации гена
			boolean flagNote = false;
			//Флаг для нотации гена, разбитой на много строк
			boolean flagNoteMultistring = false;
			//Флаг для индекса гена
			boolean flagLocus = false;
			//Флаг того, что мы находимся в конетксте (CDS)
			//UPD:2025.10.15 закомментировала, т.к. это CDS находится в контексте gene, а не наоборот
			//boolean flagContext = false;
			//Позиция конца предыдущего гена
			int lastpoz =  -2-parameters.getDistBOper();
			//Позиция предыдущего гена в цепи
			boolean lastpozcomplement = false;
			
			while (!(line = reader.readLine()).startsWith("ORIGIN")) {
				line = line.trim();
				
				//Считывание информации о гене (CDS)
				if (line.startsWith(parameters.getContext() + "   ")) {
					//2025.10.16 Если это самая первая CDS, то код ниже нужно пропустить
					int num = organism.getGenes().size();
					//2025.10.16 блок ниже нужен ТОЛЬКО для генов с множеством CDS. В противном случае,
					//не отработает проверка на равенство локус-тэгов, и это нормально.
					if (0 != num) {
						//2025.10.16 Началась новая CDS, поэтому необходимо сохранить информацию, полученную 
						//при обработке gene, соответствующего ВСЕЙ СЕРИИ CDS
						//Сохраняем информацию о /note и рибосомальности
						Gene g = organism.getGenes().get(num - 1);
						if (null != this.geneNote.note && this.geneNote.locus_tag.equals(g.getLocusTag())) {
							g.setNote(this.geneNote.note);
							if (this.geneNote.isRib && !g.isRibosomal()) {
								g.setRibosomal(true);
								//Число рибосомных генов
								organism.addRibGene();
							}
						}
					}
					//flagContext = true;
					//Число генов, порядковый номер
					organism.addGene();
					//Число оперонов
					if (!organism.isEukaryota()) {
						organism.addOperon();
					}
					
					boolean useGeneInCounting = true;
					//Наличие "complement"
					boolean complement = false;
					//Количество экзонов в гене
					int numberOfExons = 0;
					
					//Запись позиций (гена) экзонов
					String position = line.replace(parameters.getContext() + "            ", "");
					//Считывание всех позиций гена
					while (!(line = reader.readLine().trim()).startsWith("/")) {
						position += line;
					}
					//Удаление скобок () и знаков <>
					position = position.replaceAll(BRACKETS, "");
					
					//Массив стартов экзонов
					int[] start = null;
					//Массив концов экзонов
					int[] stop = null;
					
					//Если позиция начинается с идентификатора "join"
					if (position.startsWith("join")) {
						if (position.contains("complement")) {
							position = position.replace("complement", "");
							complement = true;
							useGeneInCounting = false;
						}
						
						position = position.replace("join", "");
						
						if (position.contains(",")) {
							//Массив позиций экзонов
							String[] pozOfExons = position.split(",");
							//Количество экзонов
							numberOfExons = pozOfExons.length;
							//Массив начал экзонов
							start = new int[numberOfExons];
							//Массив концов экзонов
							stop = new int[numberOfExons];
							
							for (int i = 0; i < numberOfExons; i++) {
								//Если указаны и начало и конец
								if (pozOfExons[i].contains("..")) {
									String[] arr = pozOfExons[i].split("\\.\\.");
									start[i] = Integer.parseInt(arr[0]);
									stop[i] = Integer.parseInt(arr[1]);
								}
								//Если позиция состоит из одного нуклеотида
								else {
									start[i] = Integer.parseInt(pozOfExons[i]);
									stop[i] = Integer.parseInt(pozOfExons[i]);
								}
							}
						}
					} 
					//Если позиция начинается с идентификатора "complement"
					else {
						if (position.startsWith("complement")) {
							position = position.replace("complement", "");
							complement = true;
							
							//Если позиции содержат идентификатор "join"
							if (position.startsWith("join")) {
								position = position.replace("join", "");
								
								if (position.contains(",")) {
									//Массив позиций экзонов
									String pozOfExons[] = position.split(",");
									//Количество экзонов
									numberOfExons = pozOfExons.length;
									//Массив начал экзонов
									start = new int[numberOfExons];
									//Массив концов экзонов
									stop = new int[numberOfExons];
									
									for (int i = 0; i < numberOfExons; i++) {
										//Если указаны и начало и конец
										if (pozOfExons[i].contains("..")) {
											String arr[] = pozOfExons[i].split("\\.\\.");
											start[i] = Integer.parseInt(arr[0]);
											stop[i] = Integer.parseInt(arr[1]);
										}
										//Если позиция состоит из одного нуклеотида
										else {
											start[i] = Integer.parseInt(pozOfExons[i]);
											stop[i] = Integer.parseInt(pozOfExons[i]);
											useGeneInCounting = false;
										}
									}
								}
							} 
							//Если позиции не содержат идентификатор "join"
							else {
								//Количество экзонов
								numberOfExons = 1;
								//Массив из 1-го элемента начал экзонов
								start = new int[numberOfExons];
								//Массив из 1-го элемента концов экзонов
								stop = new int[numberOfExons];
								
								//Если указаны и начало и конец
								if (position.contains("..")) {
									String[] arr = position.split("\\.\\.");
									start[0] = Integer.parseInt(arr[0]);
									stop[0] = Integer.parseInt(arr[1]);
								} 
								//Если позиция состоит из одного нуклеотида
								else {
									start[0] = Integer.parseInt(position);
									stop[0] = Integer.parseInt(position);
									useGeneInCounting = false;
								}
							}
						}
						//Если позиции не содержат идентификаторы "join" и "complement"
						else {
							//Количество экзонов
							numberOfExons = 1;
							//Массив из 1-го элемента начал экзонов
							start = new int[numberOfExons];
							//Массив из 1-го элемента концов экзонов
							stop = new int[numberOfExons];
							
							//Если указаны и начало и конец
							if (position.contains("..")) {
								String[] arr = position.split("\\.\\.");
								start[0] = Integer.parseInt(arr[0]);
								stop[0] = Integer.parseInt(arr[1]);
							} 
							//Если позиция состоит из одного нуклеотида
							else {
								start[0] = Integer.parseInt(position);
								stop[0] = Integer.parseInt(position);
								useGeneInCounting = false;
							}
						}
					}
					
					//Проверка расстояния между соседними генами
					if (!organism.isEukaryota() && !(lastpozcomplement ^ complement) && (start[0] - lastpoz - parameters.getDistBOper() - 2 < 0)) {
						organism.delOperon();
					}
					
					//Позиция конца гена
					lastpoz = stop[numberOfExons - 1];
					//Позиция гена в цепи
					lastpozcomplement = complement;
					
					//Запись гена организма
					Gene gene = new Gene(organism.getGenesNumber(), start, stop, organism.getOperonsNumber(), complement, useGeneInCounting);
					//Запись нового гена в список
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
				
				//Название гена
				if (flagCDS && line.startsWith("/gene=")) {
					/*if (p4.matcher(line).matches()) {
						l.get(l.size()-1).putR(true);
						//Число рибосомных генов
						this.ribgenes++;
					}*/
					
					//Определение рибосомных генов по названию гена
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMALG.matcher(line).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Число рибосомных генов
							organism.addRibGene();
						}
					}
					
					//Определение высоко экспрессирующихся генов по метке
					if (parameters.isUseOtherGenesAsHighExpressed1() && line.contains(" highexpressed")) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Число рибосомных генов
							organism.addRibGene();
						}
					}
					
					//Запись названия гена
					organism.getGenes().get(organism.getGenes().size() - 1).setName(line.substring("/gene=".length()).replace("\"", ""));
					flagCDS = false;
				}
				
				//Определение псевдогена
				if (flagPseudoGene && line.toLowerCase().startsWith("/pseudo")) {
					organism.addPseudoGene();
					//Отметка о псевдогене
					organism.getGenes().get(organism.getGenes().size() - 1).setPseudo(true);
					flagPseudoGene = false;
				}
				
				//Индекс гена
				if (line.startsWith("/locus_tag=")) {
					if (flagLocus) {
						organism.getGenes().get(organism.getGenes().size() - 1).setLocusTag(line.substring("/locus_tag=".length()).replace("\"", ""));
						flagLocus = false;
					} else /*if (flagGene) */{
						this.geneNote.locus_tag = line.substring("/locus_tag=".length()).replace("\"", "");
					}
				}
				
				//Название продукта гена
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
					
					//Запись названия продукта гена
					organism.getGenes().get(organism.getGenes().size() - 1).setProteinName(line1.substring("/product=".length()).replace("\"", ""));
					
					//Определение рибосомных генов по названию продукта гена
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line1).matches() && !NONRIBOSOMAL.matcher(line1).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Число рибосомных генов
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
				
				//Дозапись многострочной /note
				if (flagNoteMultistring) {
					this.geneNote.note += line.replace("\"", "");
					flagNoteMultistring = false;
				}
				
				//Запись аннотации
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
						
						//Определение рибосомных генов по аннотации гена
						if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line1).matches() && !NONRIBOSOMAL.matcher(line1).matches()) {
							if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
								organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
								//Число рибосомных генов
								organism.addRibGene();
							}
						}
						flagNote = false;
					} else /*if (flagGene)*/ {
						//TODO: добавить парсинг COG у gene по аналогии, как в блоке выше
						//Определение начала многострочного /note
						String s = line.substring("/note=".length());
						this.geneNote.note = s.replace("\"", "");
						int countBrackets = s.length() - this.geneNote.note.length();
						if (countBrackets == 1) {
							flagNoteMultistring = true;
						}
						//Определение рибосомных генов по аннотации гена
						if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line).matches() && !NONRIBOSOMAL.matcher(line).matches()) {
							this.geneNote.isRib = true;
						}
					}
				}
				
				//id гена в FLYBASE
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
					//2025.10.16 Закончилась последняя CDS для текущего gene, поэтому необходимо сохранить 
					//информацию, полученную при обработке gene, соответствующего ей (этот gene шёл ДО CDS)
					//Сохраняем информацию о /note и рибосомальности
					int num = organism.getGenes().size();
					//для первого gene код ниже нерелевантен
					if (num != 0) {
						Gene g = organism.getGenes().get(num - 1);
						if (null != this.geneNote.note && this.geneNote.locus_tag.equals(g.getLocusTag())) {
							g.setNote(this.geneNote.note);
							if (this.geneNote.isRib && !g.isRibosomal()) {
								g.setRibosomal(true);
								//Число рибосомных генов
								organism.addRibGene();
							}
						}
						//2025.10.15 Обнуляем инфу о предыдущем гене
						this.geneNote.reset();
					}
				}
				
				//Определение количества tRNA
				if (line.startsWith("tRNA          ")) {
					organism.addtRnaGene();
				}
			}
			
			//2025.10.16 Обработка самой последней CDS последнего gene 
			//Сохраняем информацию о /note и рибосомальности
			int num = organism.getGenes().size();
			//в норме проверка ниже всегда проходит
			if (num != 0) {
				Gene g = organism.getGenes().get(num - 1);
				if (null != this.geneNote.note && this.geneNote.locus_tag.equals(g.getLocusTag())) {
					g.setNote(this.geneNote.note);
					if (this.geneNote.isRib && !g.isRibosomal()) {
						g.setRibosomal(true);
						//Число рибосомных генов
						organism.addRibGene();
					}
				}
			}
			
			//============================================================================================
			//Считывание из карточки полной нуклеотидной последовательности===============================
			//TODO========================================================================================
			
			//Массив для записи полного генома
			char[] wholeGenome = new char[genomeSize];
			
			int pozstrGenome = 0;
			
			//Число GC и CG пар, GC-состав
			int gc = 0;
			
			//Считывание геномной последовательности
			while (!(line = reader.readLine()).startsWith("//")) {
				line = line.trim();
				
				String[] arr = line.split(" ");
				
				//Строка для записи геномной последовательности
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
			
			//GC-состав генома
			organism.addGC(gc);
			double gccom = Math.round(gc*10000.0/genomeSize)/100.0;
			
			int lFlang1;
			int rFlang1;
			
			//TODO Запись последовательностей генов в экземпляры класса
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
				
				//Число интронов
				int numderOfIntrons = organism.getGenes().get(i).getStart().length;
				
				//Позиция начала гена
				int strt = organism.getGenes().get(i).getStart()[0];
				
				//Позиция конца гена
				int stp = organism.getGenes().get(i).getStop()[numderOfIntrons - 1];
				int stp1 = stp + rFlang1 - genomeSize;
				int dlinaGenaBezFl = 0;
				int dlinaGena;
				
				//Строка для записи последовательности гена
				char[] strNA;
				
				//char strNA2[] = new char[2003];
				
				//Определение длины гена без флангов
				for (int j = 0; j < numderOfIntrons; j++) {
					dlinaGenaBezFl += organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
				}
				
				//Если ДНК линейная
				if (linear) {
					//Если длина левого фланга больше, чем есть на хромосоме до начала гена
					if (strt < lFlang1) {
						lFlang1 = strt - 1;
						organism.getGenes().get(i).getStart()[0] = 1;
						organism.getGenes().get(i).setFlagUse(false);
					} 
					//Если длина левого фланга меньше, чем есть на хромосоме до начала гена
					else {
						organism.getGenes().get(i).getStart()[0] -= lFlang1;
					}
					
					//Если длина правого фланга больше, чем есть на хромосоме после конца гена
					if (stp1 > 0) {
						rFlang1 = genomeSize - stp;
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] = genomeSize;
						organism.getGenes().get(i).setFlagUse(false);
					} 
					//Если длина правого фланга меньше, чем есть на хромосоме после конца гена
					else {
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] += rFlang1;
					}
					
					//Длина гена с флангами
					dlinaGena = lFlang1 + dlinaGenaBezFl + rFlang1;
					
					//Массив для записи последовательности гена с флангами
					strNA = new char [dlinaGena];
					int pozEx = 0;
					
					//Считывание последовательности гена из полного генома
					for (int j = 0; j < numderOfIntrons; j++) {
						int dlinaEx = organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
						//Строка с последовательностью гена
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
				//Если ДНК кольцевая
				else {
					//Помещается ли правый фланг
					boolean rflng = false;
					
					//Длина гена с флангами
					dlinaGena = lFlang1 + dlinaGenaBezFl + rFlang1;
					
					//Массив для записи последовательности гена с флангами
					strNA = new char [dlinaGena];
					int pozEx = 0;
					
					//Если длина левого фланга меньше, чем есть на хромосоме до начала гена
					if (strt > lFlang1) {
						//Новая позиция старта с учетом флангов
						organism.getGenes().get(i).getStart()[0] -= lFlang1;
					}
					//Если длина левого фланга больше, чем есть на хромосоме до начала гена
					else {
						//Длина недостающей последовательности до полной длины фланга
						int lenght = lFlang1 - strt + 1;
						
						//Запись недостающей последовательности с другого (правого) края кольцевой хромосомы
						//Строка с последовательностью гена
						System.arraycopy(wholeGenome, genomeSize - lenght, strNA, 0, lenght);
						organism.getGenes().get(i).getStart()[0] = 1;
						pozEx = lenght;
					}
					
					//Если длина правого фланга меньше, чем есть на хромосоме после конца гена
					if (stp1 < 0) {
						//Новая позиция конца с учетом флангов;
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] += rFlang1;
					}
					//Если длина правого фланга больше, чем есть на хромосоме после конца гена
					else {
						organism.getGenes().get(i).getStop()[numderOfIntrons - 1] = genomeSize;
						rflng = true;
					}
					
					//Считывание последовательности гена из полного генома
					for (int j = 0; j < numderOfIntrons; j++) {
						int dlinaEx = organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
						//Строка с последовательностью гена
						System.arraycopy(wholeGenome, organism.getGenes().get(i).getStart()[j] - 1, strNA, pozEx, dlinaEx);
						pozEx += dlinaEx;
					}
					
					//Если правый фланг не поместился на хромосоме, запись недостающей последовательности с другого (левого) края кольцевой хромосомы
					if (rflng) {
						//Строка с последовательностью гена
						System.arraycopy(wholeGenome, 0, strNA, pozEx, stp1);
					}
				}
				
				//Если ген в комплементарной цепи
				if (organism.getGenes().get(i).isComplement()) {
					strNA = Auxiliary_Class.reversComplement(strNA);
					//strNA2 = Auxiliary_Class.complement(strNA2);
				}
				
				//Позиция начала гена без учета флангов
				organism.getGenes().get(i).getStart()[0] = strt;
				//Позиция конца гена без учета флангов
				organism.getGenes().get(i).getStop()[numderOfIntrons - 1] = stp;
				//Запись последовательности гена
				organism.getGenes().get(i).setNA(strNA);
				//organism.genes.get(i).putNA2(strNA2);
				
				if (organism.getGenes().get(i).isComplement()) {
					int z = lFlang1;
					lFlang1 = rFlang1;
					rFlang1 = z;
				}
				
				//Запись длины гена с флангами
				organism.getGenes().get(i).putLenght(lFlang1, dlinaGenaBezFl, rFlang1, dlinaGena);
				
				//Пометка для удаления коротких генов
				if (dlinaGenaBezFl < parameters.getMinLenGen()) {
					organism.getGenes().get(i).setFlagUse(false);
					organism.addShortGene();
					//System.out.println("Gene number "+organism.genes.get(i).getNumber()+" ("+organism.genes.get(i).getName()+") was deleted because of its small length!");
				}
				
				//Пометка для удаления псевдо генов
				if (!parameters.isUsePseudoGenes() && organism.getGenes().get(i).isPseudo()) {
					organism.getGenes().get(i).setFlagUse(false);
				}
			}
			
			//Создание файла с полной геномной последовательностью
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
