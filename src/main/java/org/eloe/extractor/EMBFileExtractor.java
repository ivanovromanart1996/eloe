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
		
		//Количество генов в списке до считывания этой карточки
		int lsize = organism.getGenes().size();
		
		//Считывание файла с геномом
		BufferedReader reader = null;
	
		try {
			//Считывание файла с геномом
			reader = new BufferedReader(new FileReader(orgFile));
			
			//Линейная ДНК или кольцевая
			boolean linear = true;
			
			String line = reader.readLine();
			
			if (!line.startsWith("ID ")) {
				LOG.error("ERROR! Wrong file format or file is empty: " + orgFile + "! Right formats: *.gbk, *.gb, *.emb");
				return null;
			}
			
			//Линейность ДНК
			if (line.contains(" circular")) {
				linear = false;
			}
			String[] arr3 = line.split(" ");
			
			//Размер генома
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
			
			//Флаг для определения организма
			boolean flagOrg = true;
			
			//Определение названия организма
			do {
				line = reader.readLine().trim();
				
				//Считывание названия организма
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
	        //Считывание из карточки информации о генах===================================================
	        //TODO========================================================================================
			
			//Флаг для определения гена
			boolean flagGen = false;
			//Флаг для определения псевдогена
			boolean flagPseudoGen = false;
			//Флаг для определения продукта гена
			boolean flagProduct = false;
			//Флаг для определения индекса белка
			boolean flagProtID = false;
			//Флаг для определения индекса гена
			boolean flagDBxref = false;
			//Флаг для записи аннотации
			boolean flagNote = false;
			//Флаг для индекса гена
			boolean flagLocus = false;
			//Флаг для определения номера оперона (нужен только для первого прохода)
			boolean flagOper = false;
			//Позиция конца предыдущего гена
			int lastpoz = 0;
			//Позиция предыдущего гена в цепи
			boolean lastpozcomplement = false;
			String context_new = "FT   " + parameters.getContext();
			
			while (!(line = reader.readLine()).startsWith("SQ   ")) {
				line = line.trim();
				
				//Считывание информации о гене
				if (line.startsWith(context_new)) {
					//Число генов, порядковый номер
					organism.addGene();
					//Число оперонов
					if (!organism.isEukaryota()) {
						organism.addOperon();
					}
					
					//Флаг для определения наличия интронов
					boolean flagIntr = true;
					//Наличие "join"
					boolean join = false;
					//Наличие "complement"
					boolean complement = false;
					
					//Запись позиций (гена) экзонов
					String position = line;
					
					//Информация о комплементарности
					if (position.contains("complement")) {
						complement = true;
						position = position.replace(context_new + "            complement", "");
						if (position.contains("join")) {
							join = true;
							//Информация о наличии интронов
							position = position.replace("join", "");
						}
						flagIntr = false;
					}
					
					//Информация о наличии интронов
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
					
					//Удаление скобок () и значков <>
					position = position.replaceAll(BRACKETS.toString(), "");
					
					//Если несколько экзонов в гене
					if (position.contains(",")) {
						
						//Массив позиций экзонов
						String[] pozOfExons = position.split(",");
						
						//Количество экзонов
						int numberOfExons = pozOfExons.length;
						
						//Массив начал экзонов
						int[] begining = new int[numberOfExons];
						
						//Массив концов экзонов
						int[] end = new int[numberOfExons];
						for (int i = 0; i < numberOfExons; i++) {
							
							if (pozOfExons[i].contains("..")) {
								String arr[] = pozOfExons[i].split("\\.\\.");
								begining[i] = Integer.parseInt(arr[0]);
								end[i] = Integer.parseInt(arr[1]);
							} else { //Если экзон состоит из одного нуклеотида
								//Позиции единственного нуклеотида в экзоне
								begining[i] = Integer.parseInt(pozOfExons[i]);
								end[i] = Integer.parseInt(pozOfExons[i]);
							}
						}
						
						//Проверка расстояния между соседними генами
						if (flagOper && !organism.isEukaryota() && !(lastpozcomplement ^ complement) && (begining[0] - lastpoz - parameters.getDistBOper() - 2) < 0) {
							organism.delOperon();
						}
						
						//Позиция конца гена
						lastpoz = end[numberOfExons - 1];
						lastpozcomplement = complement;
						
						//Нужен был для первого прохода
						flagOper = true;
						Gene gene = new Gene(organism.getGenesNumber(), begining, end, organism.getOperonsNumber(), complement, true);
						
						//Запись нового гена в список
						organism.getGenes().add(gene);
					} else { //Если один экзон в гене
						String[] arr = position.split("\\.\\.");
						
						//Массив из 1-го элемента начал экзонов
						int[] begining = {Integer.parseInt(arr[0])};
						
						//Массив из 1-го элемента концов экзонов
						int[] end = {Integer.parseInt(arr[1])};
						
						//Проверка расстояния между соседними генами
						if (flagOper && !organism.isEukaryota() && !(lastpozcomplement ^ complement) && (begining[0] - lastpoz - parameters.getDistBOper() - 2) < 0) {
							organism.delOperon();
						}
						
						//Позиция конца гена
						lastpoz = end[0];
						lastpozcomplement = complement;
						
						//Нужен был для первого прохода
						flagOper = true;
						Gene gene = new Gene(organism.getGenesNumber(), begining, end, organism.getOperonsNumber(), complement, true);
						
						//Запись нового гена в список
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
				
				//Определение псевдогена
				if (flagPseudoGen && line.startsWith("FT                   /pseudo")) {
					organism.addPseudoGene();
					
					//Отметка о псевдогене
					organism.getGenes().get(organism.getGenes().size() - 1).setPseudo(true);
					flagPseudoGen = false;
				}
				
				//Название гена
				if (flagGen && line.startsWith("FT                   /gene=")) {
					line = line.replace("FT                   ", "");
					
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMALG.matcher(line).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Число рибосомных генов
							organism.addRibGene();
						}
					}
					
					//Определение высоко экспрессирующихся генов по метке
					if (parameters.isUseOtherGenesAsHighExpressed1() && line.contains(" highexpressed")) {
						line = line.replace(" highexpressed", "");
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Число рибосомных генов
							organism.addRibGene();
						}
					}
					
					//Запись названия гена
					organism.getGenes().get(organism.getGenes().size() - 1).setName(line.replace("FT                   ", ""));
					flagGen = false;
				}
				
				//Индекс гена locus_tag=
				if (flagLocus && line.startsWith("FT                   /locus_tag=")) {
					//32 - длина слова "FT                   /locus_tag="
					organism.getGenes().get(organism.getGenes().size() - 1).setLocusTag(line.substring(32).replace("\"", ""));
					flagLocus = false;
				}
				
				//Название продукта гена
				if (flagProduct && line.startsWith("FT                   /product=")) {
					if (line.contains("COG")) {
						Matcher m2 = COGID.matcher(line);
						if (m2.find() && organism.getGenes().size() != 0) {
							organism.getGenes().get(organism.getGenes().size() - 1).setCog(m2.group());
						}
					}
					
					line = line.replace("FT                   ", "");
					line = line.toLowerCase();
					
					//Запись названия продукта гена
					organism.getGenes().get(organism.getGenes().size() - 1).setProteinName(line);
					
					//Определение рибосомных генов по названию продукта гена
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line).matches() && !NONRIBOSOMAL.matcher(line).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							//Число рибосомных генов
							organism.addRibGene();
						}
					}
					flagProduct = false;
					flagGen = false;
				}
				
				//protein_id=
				if (flagProtID && line.startsWith("FT                   /protein_id=")) {
					//33-длина слова "FT                   /protein_id="
					organism.getGenes().get(organism.getGenes().size() - 1).setProteinID(line.substring(33).replace("\"", ""));
					flagProtID = false;
				}
				
				//db_xref=\"GI:
				if (flagDBxref && line.startsWith("FT                   /db_xref=\"GI:")) {
					//30-длина слова "FT                   /db_xref="
					organism.getGenes().get(organism.getGenes().size() - 1).setGeneID(line.substring(30).replace("\"", ""));
					flagDBxref = false;
				}
				
				//Запись аннотации
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
					organism.getGenes().get(organism.getGenes().size() - 1).setNote(line1.substring(27).replace("\"", "")); //27 - длина "FT                   /note="
					
					//Определение рибосомных генов по аннотации гена
					if (parameters.isUseRibGenesAsHighExpressed1() && RIBOSOMAlP.matcher(line1).matches() && !NONRIBOSOMAL.matcher(line1).matches()) {
						if (!organism.getGenes().get(organism.getGenes().size() - 1).isRibosomal()) {
							organism.getGenes().get(organism.getGenes().size() - 1).setRibosomal(true);
							
							//Число рибосомных генов
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
			//Считывание из карточки полной нуклеотидной последовательности===============================
			//TODO========================================================================================
			
			//Массив для записи полного генома
			char[] wholeGenome = new char[genomeSize];
			
			int pozstrGenom = 0;
			
			//Число GC и CG пар, GC-состав
			int gc = 0;
			
			//Считывание геномной последовательности
			while (!(line = reader.readLine()).startsWith("//")) {
				line = line.trim();
				String[] arr = line.split(" ");
				
				//Строка для записи геномной последовательности
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
			
			//GC-состав генома
			organism.addGC(gc);
			double gccom = Math.round(gc*10000.0/genomeSize)/100.0;
			int lFlang1;
			int rFlang1;
			
			//Запись последовательностей генов в экземпляры класса
			for (int i = lsize; i < organism.getGenes().size(); i++) {
				if (organism.getGenes().get(i).isComplement()) {
					lFlang1 = parameters.getrFlang();
					rFlang1 = parameters.getlFlang();
				} else {
					lFlang1 = parameters.getlFlang();
					rFlang1 = parameters.getrFlang();
				}
				
				//Число интронов
				int k = organism.getGenes().get(i).getStart().length;
				
				//Позиция начала гена
				int strt = organism.getGenes().get(i).getStart()[0];
				
				//Позиция конца гена
				int stp = organism.getGenes().get(i).getStop()[k - 1];
				int stp1 = stp + rFlang1 - genomeSize;
				int dlinaGenaBezFl = 0;
				int dlinaGena;
				
				//Строка для записи последовательности гена
				char[] strNA;
				
				//Определение длины гена без флангов
				for (int j = 0; j < k; j++) {
					dlinaGenaBezFl += organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
				}
				
				//Если ДНК линейная
				if (linear) {
					//Если длина левого фланга больше, чем есть на хромосоме до начала гена
					if (strt < lFlang1) {
						lFlang1 = strt - 1;
						organism.getGenes().get(i).getStart()[0] = 1;
						organism.getGenes().get(i).setFlagUse(false);
					} else { //Если длина левого фланга меньше, чем есть на хромосоме до начала гена
						organism.getGenes().get(i).getStart()[0] -= lFlang1;
					}
					
					//Если длина правого фланга больше, чем есть на хромосоме после конца гена
					if (stp1 > 0) {
						rFlang1 = genomeSize - stp;
						organism.getGenes().get(i).getStop()[k - 1] = genomeSize;
						organism.getGenes().get(i).setFlagUse(false);
					} else { //Если длина правого фланга меньше, чем есть на хромосоме после конца гена
						organism.getGenes().get(i).getStop()[k - 1] += rFlang1;
					}
					
					//Длина гена с флангами
					dlinaGena = lFlang1 + dlinaGenaBezFl + rFlang1;
					
					//Массив для записи последовательности гена с флангами
					strNA = new char [dlinaGena];
					int pozEx = 0;
					
					//Считывание последовательности гена из полного генома
					for (int j = 0; j < k; j++) {
						int dlinaEx = organism.getGenes().get(i).getStop()[j] - organism.getGenes().get(i).getStart()[j] + 1;
						
						//Строка с последовательностью гена
						System.arraycopy(wholeGenome, organism.getGenes().get(i).getStart()[j] - 1, strNA, pozEx, dlinaEx);
						pozEx += dlinaEx;
					}
				} else { //Если ДНК кольцевая
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
					} else { //Если длина левого фланга больше, чем есть на хромосоме до начала гена
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
						organism.getGenes().get(i).getStop()[k - 1] += rFlang1;
					} else { //Если длина правого фланга больше, чем есть на хромосоме после конца гена
						organism.getGenes().get(i).getStop()[k - 1] = genomeSize;
						rflng = true;
					}
					
					//Считывание последовательности гена из полного генома
					for (int j = 0; j < k; j++) {
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
					//Создание комплементарной последовательности
					strNA = Auxiliary_Class.reversComplement(strNA);
				}
				
				//Позиция начала гена без учета флангов
				organism.getGenes().get(i).getStart()[0] = strt;
				
				//Позиция конца гена без учета флангов
				organism.getGenes().get(i).getStop()[k - 1] = stp;
				
				//Запись последовательности гена
				organism.getGenes().get(i).setNA(strNA);
				
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
				
				//Пометка для удаление псевдо генов
				if (!parameters.isUsePseudoGenes() && organism.getGenes().get(i).isPseudo()) {
					organism.getGenes().get(i).setFlagUse(false);
				}
			}
			
			//Создание файла с полной геномной последовательностью
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
