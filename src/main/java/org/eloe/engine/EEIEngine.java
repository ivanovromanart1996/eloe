package org.eloe.engine;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JTextArea;

import org.apache.log4j.Logger;

import org.eloe.engine.ta.FreqCounter;
import org.eloe.engine.te.Te_Counter;
import org.eloe.engine.te.lci.LCI12_Counter;
import org.eloe.engine.te.lci.LCI1Counter;
import org.eloe.engine.te.lci.LCI2_Counter;
import org.eloe.engine.te.lciij.LCIijDataProcessing;
import org.eloe.engine.te.lciij.LCIij_Counter;
import org.eloe.extractor.processing.Processing;
import org.eloe.organism.Organism;
import org.eloe.organism.OrganismBuilder;
import org.eloe.other.Auxiliary_Class;
import org.eloe.other.Orthologs;
import org.eloe.other.Sorter;
import org.eloe.parameters.Parameters;
import org.eloe.parameters.ParametersStructure;

public abstract class EEIEngine {
	
	private static final Logger LOG = Logger.getLogger(EEIEngine.class);
	
	//утилита для подачи в listFiles и получения всех папок из текущей папки dir/name
	private static class DirFilter implements FilenameFilter{
		@Override
		public boolean accept(File dir, String name) {
			File f = new File(dir, name);
			if (f.isDirectory()) {
				return true;
			}
			return false;
		}
	}
	
	public static void engine(ParametersStructure ps, File sourceFolder, File resultFolder, JTextArea console) {
		LOG.debug("Start EEIEngine.engine().");
		
		//Стартовая точка для определения времени работы всей программы
		long timeStart = System.currentTimeMillis();
		//Переменная для определения времени работы частей программы
		long timer = System.currentTimeMillis();
		
		//Список файлов с результатами для каждого организма
		List<String> files = new ArrayList<String>();
		
		//TODO Считывание параметров==================================================================
		LOG.debug("Parameters reading...");
		Parameters parameters = new Parameters(ps);
		
		//Читаем список папок в папке "use"
		File listFolders[] = sourceFolder.listFiles(new DirFilter());
		
		//Файл с индексами организмов
		File orgIndexes = new File(resultFolder + File.separator + "organisms_index.txt");
		Path orgIndexesPath = orgIndexes.toPath();
		PrintStream orgInd = null;
		try {
			if (Files.notExists(orgIndexesPath)) {
				orgInd = new PrintStream(new FileOutputStream(orgIndexes, true));//true - чтобы аппендить
				orgInd.println("Organism\tIndx\tM1\tR1\tM2\tR2\tM3\tR3\tM4\tR4\tM5\tR5\tNumber of all genes"
						+ "\tNumber of ribos. genes\tGenome length, bp\tGC composition, %\ttRNA");
			} else {
				//N.B. Не выносить наверх, иначе проверка будет бессмысленна, ибо файл будет так и так создаваться
				orgInd = new PrintStream(new FileOutputStream(orgIndexes, true));//true - чтобы аппендить
			}
			
			//Проход по всем папкам (организмам) из папки "Use"
			for (int folderNumber = 0; folderNumber < listFolders.length; folderNumber++) {
				//Таймер для определения времени просчета одного организма
				long timerOrganism = System.currentTimeMillis();
				
				//Вывод в консоль порядкового номера и имени папки (организма) для расчетов
				String say = (folderNumber + 1) + " of " + listFolders.length + "\n" + listFolders[folderNumber].getName();
				console.append(say);
				LOG.info(say);
	        	
	        	//Создание новой папки в результатах, если её там нет. Если есть - пропускаем.
				String curResFolderName = listFolders[folderNumber].getName().replace(" ", "_");
	        	File resultFolder1 = new File(resultFolder, curResFolderName);
	        	if (!resultFolder1.exists()) {
	        		resultFolder1.mkdir();
	        	} else {
	        		continue;//если папка уже есть - пропускаем.
	        	}
	        	
	        	//Стартовая точка для определения времени работы частей программы 
	        	timer = System.currentTimeMillis();
	    		
	    		//TODO Extraction of organism data from directory
	        	LOG.debug("Extraction...");
	        	Organism organism = OrganismBuilder.getOrganism(listFolders[folderNumber], parameters, resultFolder1);
	        	if (organism == null) {
	        		continue;
	        	}
	        	//Определение времени экстракции
	        	timer = Auxiliary_Class.time("1.1st stage - Extract", timer);
	        	
	        	//TODO Processing
	        	LOG.debug("Processing...");
	        	//Класс "Processing", выбрасывание из списка "плохих" генов
				Processing.processing(organism, parameters);
				//Определение времени процессинга
				timer = Auxiliary_Class.time("1.2st stage - Processing", timer);
	        	
	        	//Количество генов в обучающей выборке
				int number2 = parameters.teachingGenes(organism.getGenes().size());
				
				//TODO Taking into account codon frequency
				LOG.debug("Codon frequency counting...");
				countCodonFrequencies(organism, parameters, resultFolder1, number2, console);
				timer = Auxiliary_Class.time("2nd stage - Freq", timer);
				
				double m1 = 0;
				double sigma1 = 0;
				double m2 = 0;
				double sigma2 = 0;
				
				//TODO Taking into account local perfect inverted repeats in mRNA
				LOG.debug("LCI counting...");
				{
					//Если параметры расчета петель для LCI1 и LCI2 совпадают
					if (parameters.getSmin1() == parameters.getSmin2() 
							&& parameters.getSmax1() == parameters.getSmax2() 
							&& parameters.getLmin1() == parameters.getLmin2() 
							&& parameters.getLmax1() == parameters.getLmax2())
					{
						//Класс "LCI12_Counter"
						LCI12_Counter lci12 = null;
						try {
							lci12 = new LCI12_Counter(parameters);
						} catch (IOException e) {
							e.printStackTrace();
							console.append(e.toString());
						} catch (NumberFormatException e) {
							e.printStackTrace();
							console.append(e.toString());
						}
						lci12.countLCI(organism);
						//Среднее значение LCI1
						m1 = lci12.getM1();
						//Среднее отклонение LCI1
						sigma1 = lci12.getSigma1();
						//Среднее значение LCI2
						m2 = lci12.getM2();
						//Среднее отклонение LCI2
						sigma2 = lci12.getSigma2();
						//Определение времени расчета LCI12
						timer = Auxiliary_Class.time("3th stage - LCI12", timer);
					}
					else
					{
						//Класс "LCI1_Counter"
						LCI1Counter lci1 = new LCI1Counter(parameters);
						lci1.countLCI(organism.getGenes());
						//Среднее значение LCI1
						m1 = lci1.getM();
						//Среднее отклонение LCI1
						sigma1 = lci1.getSigma();
						//Определение времени расчета LCI1
						timer = Auxiliary_Class.time("3.1th stage - LCI1", timer);
									
						//Класс "LCI2_Counter"
						LCI2_Counter lci2 = null;
						try {
							lci2 = new LCI2_Counter(parameters);
						} catch (IOException e) {
							e.printStackTrace();
							console.append(e.toString());
						} catch (NumberFormatException e) {
							e.printStackTrace();
							console.append(e.toString());
						}
						lci2.countLCI(organism.getGenes());
						//Среднее значение LCI2
						m2 = lci2.getM();
						//Среднее отклонение LCI2
						sigma2 = lci2.getSigma();
						//Определение времени расчета LCI2
						timer = Auxiliary_Class.time("3.2th stage - LCI2", timer);
					}
				}
				
				// For counting average number of repeats in all no ribosomal genes
				if (ps.countRepeatsInGenes) {
					countAvNumbOfRepeats(organism);
				}
							
				//TODO Класс "LCIij_Counter"
				LOG.debug("LCIij counting...");
				if (parameters.isCountLCIij())
				{
					LCIij_Counter lciij = new LCIij_Counter(parameters, resultFolder1 + File.separator + organism.getName() + "_lciij.txt");
					//Подсчет LCIij
					lciij.countLCIij(organism);
					timer = Auxiliary_Class.time("3.3th stage - LCIij", timer);
				}
				
				//TODO Класс "Te_Counter", расчет Tei
				LOG.debug("Te counting...");
				if (parameters.isUseUserParam1()) {
					m1 = parameters.getM11();
					sigma1 = parameters.getSigma11();
					m2 = parameters.getM22();
					sigma2 = parameters.getSigma22();
				}
				
				//Вычисление p(i)1 //без учета энергии (false)
				Te_Counter.countTei(organism.getGenes(), m1, sigma1, false);
				//Вычисление p(i)2 //с учетом энергии (true)
				Te_Counter.countTei(organism.getGenes(), m2, sigma2, true);
				//Определение времени расчета Tei
				timer = Auxiliary_Class.time("4th stage - Tei", timer);
				
				//Запись заголовка в файле с позициями рибосомных генов в списке
				File fileRibPoz = new File(resultFolder1, organism.getName() + "_ribpos.txt");
				writeRibPosFileHead(fileRibPoz,	organism.getName());
				
				//TODO Класс "EEI_Counter"
				LOG.debug("EEI counting...");
				EEI_Counter eei = new EEI_Counter(number2, parameters, fileRibPoz);
				//Массив для индексов (5 - типов индекса, 2 - параметры M и R)
				int[][] mr = eei.countEEI(organism, console);
				//Определение времени расчета EEI
				timer = Auxiliary_Class.time("6th stage - EEI", timer);
				
				//Считывание файла с позициями рибосомных генов
				String ribpoz = readRibPosFile(fileRibPoz);
				if (!ribpoz.contains("//")) {
					if (ps.showRibPos) {
						//Графический вывод положения рибосомных генов
						Auxiliary_Class.ribPosGraphic(fileRibPoz, organism);
					}
					
					//Определение индекса, который работает в данном организме
					int index = Auxiliary_Class.indexCounter(mr);
					organism.setIndex(index);
					
					//Вывод в консоль максимального значения M и минимального R
					say = "Max index: M" + (index + 1) + "=" + mr[index][0] + "\tR" + (index + 1) + "=" + mr[index][1];
					console.append("\n" + say + "\n");
					LOG.info(say);
					
					//Запись в лог-файл
					orgInd.println(organism.getName() + "\t" + (index + 1) + 
							"\t" + mr[0][0] + "\t" + mr[0][1] + "\t" + mr[1][0] + "\t" + mr[1][1] + 
							"\t" + mr[2][0] + "\t" + mr[2][1] + "\t" + mr[3][0] + "\t" + mr[3][1] + 
							"\t" + mr[4][0] + "\t" + mr[4][1] + "\t" + organism.getGenesNumber() + 
							"\t" + organism.getRibGenesNumber() + "\t" + organism.getGenomeSize() + 
							"\t" + organism.getGcComposition() + "\t" + organism.gettRnaGenesNumber());
					
					
					//Вывод в консоль индекса, который работает в данном организме
					say = "index = " + (index + 1);
					console.append(say);
					LOG.info(say);
					
					//Вывод результатов в файл
					LOG.debug("Saveing EEI to files...");
					Auxiliary_Class.outEEI(resultFolder1, organism);
					
					//Сортировка списка по увеличению номера гена
					Sorter.sortNup(organism);
					
					//Запись результатов в файлы
					LOG.debug("Saveing results to files...");
					Auxiliary_Class.resultsToFiles(organism, parameters, resultFolder1, files);
					
					//Определение времени вывода результатов
					timer = Auxiliary_Class.time("7th stage - Out", timer);
					//TODO Построение профилей LCIij
					LOG.debug("Counting LCIij profiles...");
					timer = countLCIij(parameters.isCountLCIij(), resultFolder1, organism.getName(), timer);
					//*/
					
					//Конец рассмотрения одной папки (организма), переход к рассмотрению следующей
					timerOrganism = Auxiliary_Class.time(listFolders[folderNumber].getName(), timerOrganism);
					console.append("\n\n");
				} else {
					LOG.warn("There are no ribosomal proteins in organism");
					///Чтоб для вирусов работало.
					Auxiliary_Class.numberEEItoFile(organism, resultFolder1);
					//Определение времени вывода результатов
					timer = Auxiliary_Class.time("7th stage - Out", timer);
					//TODO Построение профилей LCIij
					LOG.debug("Counting LCIij profiles...");
					timer = countLCIij(parameters.isCountLCIij(), resultFolder1, organism.getName(), timer);
					//*/
					
					//Конец рассмотрения одной папки (организма), переход к рассмотрению следующей
					timerOrganism = Auxiliary_Class.time(listFolders[folderNumber].getName(), timerOrganism);
					console.append("\n\n");
					//continue;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (orgInd != null) {
				orgInd.close();
			}
		}
		
		//Поиск ортологов
		if (parameters.isSearchOrthologs()) {
			LOG.debug("Searching orthologs...");
			Orthologs.search(files, resultFolder, parameters, console);
			//Определение времени поиска ортологов
			timer = Auxiliary_Class.time("10th stage - Orthologs", timer);
		}
		
		//Определение времени работы всей программы
		Auxiliary_Class.time("All", timeStart);
		
		//Завершение работы программы
		console.append("Done!");
		LOG.debug("Stop EEIEngine.engine().");
	}

	private static void countAvNumbOfRepeats(Organism organism) {
		double srkolvonorib = 0;
		int numbOfNoRib = 0;
		double srkolvorib = 0;
		int numbOfRib = 0;
		for (int i = 0; i < organism.getGenes().size(); i++) {
			if (!organism.getGenes().get(i).isRibosomal()) {
				numbOfNoRib++;
				srkolvonorib += organism.getGenes().get(i).getLci1()*(organism.getGenes().get(i).getLength()-3);
			}
			else {
				numbOfRib++;
				srkolvorib += organism.getGenes().get(i).getLci1()*(organism.getGenes().get(i).getLength()-3);
			}
		}
		LOG.info("Average number of repeats in all no rib genes: " + srkolvonorib/numbOfNoRib
				+"\tin all rib genes: "+srkolvorib/numbOfRib);
	}

	private static String readRibPosFile(File fileRibPoz) {
		String ribpoz = null;
		BufferedReader readerRibPoz = null;
		try {
			readerRibPoz = new BufferedReader(new FileReader(fileRibPoz));
			ribpoz = readerRibPoz.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (readerRibPoz != null) {
				try {
					readerRibPoz.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		return ribpoz;
	}

	private static void writeRibPosFileHead(File fileRibPoz, String organism_name) {
		PrintStream outRibPoz = null;
		try {				
			outRibPoz = new PrintStream(fileRibPoz);
			outRibPoz.println("Positions of ribosomal genes " + organism_name);
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (outRibPoz != null) {
				outRibPoz.close();
			}
		}
	}

	private static void countCodonFrequencies(Organism organism,
			Parameters parameters, File resultFolder1, int number2,
			JTextArea console) {
		//Подсчет частот кодонов в генах без выбрасывания флангов
		FreqCounter freq = null;
		try {
			freq = new FreqCounter();
		} catch (IOException e) {
			console.append(e.toString());
		} catch (NumberFormatException e) {
			console.append(e.toString());
		}
		freq.countFreq0(organism.getGenes());
		if (parameters.isMakeFileFreq()) {
			freq.makeFileFreq(organism.getGenes(), organism.getName(), resultFolder1);
		}
		//Автоматический подбор флангов (если задано)
		Auxiliary_Class.flangSizeAutocounting(parameters, organism, number2, freq);
		//Cортировка по увеличению порядкового номера гена. Приведение списка генов в первоначальный вид
		Sorter.sortNup(organism);
		//Подсчет частот кодонов в генах без учета кодонов на краях
		freq.countFreq(organism.getGenes(), parameters.getUnUseCodonsL1(), parameters.getUnUseCodonsR1());
		//Определение времени расчета частот кодонов и размеров флангов (если задано)
	}
	
	private static long countLCIij(final boolean countLCIij, final File resultFolder1, final String organismName, final long timer) {
		long timerLocal = timer;
		//Построение профилей LCIij
		if (countLCIij) {
			//Исходный файл с LCIij
			File fileInput1 = new File(resultFolder1, organismName + "_lciij.txt");
			//Исходный файл с EEI
			File fileInput2 = new File(resultFolder1, organismName + "_number_eei.txt");
			//Файл с результирующими профилями
			File fileOutput = new File(resultFolder1, organismName + "_lciij_profile.xlsx");
			
			//Расчеты профилей
			LCIijDataProcessing.readLCIij(fileInput1, fileInput2, fileOutput);
			//Определение времени расчета LCIij_profile
			timerLocal = Auxiliary_Class.time("8th stage - LCIij_profile", timerLocal);
			
			//Определение времени перезаписи InOut
			timerLocal = Auxiliary_Class.time("9th stage - InOut", timerLocal);
		}
		return timerLocal;
	}
		
}
