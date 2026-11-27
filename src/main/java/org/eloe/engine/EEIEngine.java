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
	
	//Filter for listFiles that returns subdirectories of dir/name.
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
		
		//Start time for the entire run
		long timeStart = System.currentTimeMillis();
		//Timer for individual stages
		long timer = System.currentTimeMillis();
		
		//List of files with results for each organism
		List<String> files = new ArrayList<String>();
		
		//TODO Read parameters==================================================================
		LOG.debug("Parameters reading...");
		Parameters parameters = new Parameters(ps);
		
		//Read the list of folders in the "use" folder.
		File listFolders[] = sourceFolder.listFiles(new DirFilter());
		
		//File of organism indices
		File orgIndexes = new File(resultFolder + File.separator + "organisms_index.txt");
		Path orgIndexesPath = orgIndexes.toPath();
		PrintStream orgInd = null;
		try {
			if (Files.notExists(orgIndexesPath)) {
				orgInd = new PrintStream(new FileOutputStream(orgIndexes, true));//true enables append mode
				orgInd.println("Organism\tIndx\tM1\tR1\tM2\tR2\tM3\tR3\tM4\tR4\tM5\tR5\tNumber of all genes"
						+ "\tNumber of ribos. genes\tGenome length, bp\tGC composition, %\ttRNA");
			} else {
				//N.B. Keep file creation inside this branch so the existence check remains meaningful.
				orgInd = new PrintStream(new FileOutputStream(orgIndexes, true));//true enables append mode
			}
			
			//Process each organism subdirectory in "Use".
			for (int folderNumber = 0; folderNumber < listFolders.length; folderNumber++) {
				//Timer to determine the timing of the calculation of one organism
				long timerOrganism = System.currentTimeMillis();
				
				//Output in the console of the sequence number and the name of the folder (organism) for calculations
				String say = (folderNumber + 1) + " of " + listFolders.length + "\n" + listFolders[folderNumber].getName();
				console.append(say);
				LOG.info(say);
	        	
	        	//Create an output subdirectory. Skip the organism if it already exists.
				String curResFolderName = listFolders[folderNumber].getName().replace(" ", "_");
	        	File resultFolder1 = new File(resultFolder, curResFolderName);
	        	if (!resultFolder1.exists()) {
	        		resultFolder1.mkdir();
	        	} else {
	        		continue;//Skip if the output subdirectory already exists.
	        	}
	        	
	        	//Start time for the current stage 
	        	timer = System.currentTimeMillis();
	    		
	    		//TODO Extraction of organism data from directory
	        	LOG.debug("Extraction...");
	        	Organism organism = OrganismBuilder.getOrganism(listFolders[folderNumber], parameters, resultFolder1);
	        	if (organism == null) {
	        		continue;
	        	}
	        	//Determination of extraction time
	        	timer = Auxiliary_Class.time("1.1st stage - Extract", timer);
	        	
	        	//TODO Processing
	        	LOG.debug("Processing...");
	        	//The Processing class removes unwanted genes from the list.
				Processing.processing(organism, parameters);
				//Determination of processing time
				timer = Auxiliary_Class.time("1.2st stage - Processing", timer);
	        	
	        	//Number of genes in the training sample
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
					//If the loop parameters for LCI1 and LCI2 match
					if (parameters.getSmin1() == parameters.getSmin2() 
							&& parameters.getSmax1() == parameters.getSmax2() 
							&& parameters.getLmin1() == parameters.getLmin2() 
							&& parameters.getLmax1() == parameters.getLmax2())
					{
						//Class "LCI12_Counter"
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
						//Average LCI1 value
						m1 = lci12.getM1();
						//Average deviation of LCI1
						sigma1 = lci12.getSigma1();
						//Average value of LCI2
						m2 = lci12.getM2();
						//Average deviation of LCI2
						sigma2 = lci12.getSigma2();
						//Determination of the calculation time of LCI12
						timer = Auxiliary_Class.time("3th stage - LCI12", timer);
					}
					else
					{
						//Class "LCI1_Counter"
						LCI1Counter lci1 = new LCI1Counter(parameters);
						lci1.countLCI(organism.getGenes());
						//Average LCI1 value
						m1 = lci1.getM();
						//Average deviation of LCI1
						sigma1 = lci1.getSigma();
						//Determination of the calculation time of LCI1
						timer = Auxiliary_Class.time("3.1th stage - LCI1", timer);
									
						//Class "LCI2_Counter"
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
						//Average value of LCI2
						m2 = lci2.getM();
						//Average deviation of LCI2
						sigma2 = lci2.getSigma();
						//Determination of the calculation time of LCI2
						timer = Auxiliary_Class.time("3.2th stage - LCI2", timer);
					}
				}
				
				// For counting average number of repeats in all no ribosomal genes
				if (ps.countRepeatsInGenes) {
					countAvNumbOfRepeats(organism);
				}
							
				//TODO LCIij_Counter class
				LOG.debug("LCIij counting...");
				if (parameters.isCountLCIij())
				{
					LCIij_Counter lciij = new LCIij_Counter(parameters, resultFolder1 + File.separator + organism.getName() + "_lciij.txt");
					//Counting LCIij
					lciij.countLCIij(organism);
					timer = Auxiliary_Class.time("3.3th stage - LCIij", timer);
				}
				
				//TODO Class "Te_Counter", calculation Tei
				LOG.debug("Te counting...");
				if (parameters.isUseUserParam1()) {
					m1 = parameters.getM11();
					sigma1 = parameters.getSigma11();
					m2 = parameters.getM22();
					sigma2 = parameters.getSigma22();
				}
				
				//Calculation p(i)1 // excluding energy (false)
				Te_Counter.countTei(organism.getGenes(), m1, sigma1, false);
				//Calculation of p(i)2 // with energy (true)
				Te_Counter.countTei(organism.getGenes(), m2, sigma2, true);
				//Determination of Tei calculation time
				timer = Auxiliary_Class.time("4th stage - Tei", timer);
				
				//Write the header to the file containing ribosomal gene positions.
				File fileRibPoz = new File(resultFolder1, organism.getName() + "_ribpos.txt");
				writeRibPosFileHead(fileRibPoz,	organism.getName());
				
				//TODO Class "EEI_Counter"
				LOG.debug("EEI counting...");
				EEI_Counter eei = new EEI_Counter(number2, parameters, fileRibPoz);
				//Array for indices (5 index types, 2 M and R parameters)
				int[][] mr = eei.countEEI(organism, console);
				//Determination of the EEI calculation time
				timer = Auxiliary_Class.time("6th stage - EEI", timer);
				
				//Reading a file with positions of ribosomal genes
				String ribpoz = readRibPosFile(fileRibPoz);
				if (!ribpoz.contains("//")) {
					if (ps.showRibPos) {
						//Plot ribosomal gene positions.
						Auxiliary_Class.ribPosGraphic(fileRibPoz, organism);
					}
					
					//Determination of the index that works in this organism
					int index = Auxiliary_Class.indexCounter(mr);
					organism.setIndex(index);
					
					//Output in the console of maximum value M and minimum R
					say = "Max index: M" + (index + 1) + "=" + mr[index][0] + "\tR" + (index + 1) + "=" + mr[index][1];
					console.append("\n" + say + "\n");
					LOG.info(say);
					
					//Log file entry
					orgInd.println(organism.getName() + "\t" + (index + 1) + 
							"\t" + mr[0][0] + "\t" + mr[0][1] + "\t" + mr[1][0] + "\t" + mr[1][1] + 
							"\t" + mr[2][0] + "\t" + mr[2][1] + "\t" + mr[3][0] + "\t" + mr[3][1] + 
							"\t" + mr[4][0] + "\t" + mr[4][1] + "\t" + organism.getGenesNumber() + 
							"\t" + organism.getRibGenesNumber() + "\t" + organism.getGenomeSize() + 
							"\t" + organism.getGcComposition() + "\t" + organism.gettRnaGenesNumber());
					
					
					//Output in the console of the index that works in this organism
					say = "index = " + (index + 1);
					console.append(say);
					LOG.info(say);
					
					//Output of results into file
					LOG.debug("Saveing EEI to files...");
					Auxiliary_Class.outEEI(resultFolder1, organism);
					
					//Sort the list by ascending gene number.
					Sorter.sortNup(organism);
					
					//Record the results in files
					LOG.debug("Saveing results to files...");
					Auxiliary_Class.resultsToFiles(organism, parameters, resultFolder1, files);
					
					//Determination of output time
					timer = Auxiliary_Class.time("7th stage - Out", timer);
					//TODO Building LCIij profiles
					LOG.debug("Counting LCIij profiles...");
					timer = countLCIij(parameters.isCountLCIij(), resultFolder1, organism.getName(), timer);
					//*/
					
					//Finish the current organism folder and move to the next one.
					timerOrganism = Auxiliary_Class.time(listFolders[folderNumber].getName(), timerOrganism);
					console.append("\n\n");
				} else {
					LOG.warn("There are no ribosomal proteins in organism");
					///Enable processing for viruses.
					Auxiliary_Class.numberEEItoFile(organism, resultFolder1);
					//Determination of output time
					timer = Auxiliary_Class.time("7th stage - Out", timer);
					//TODO Building LCIij profiles
					LOG.debug("Counting LCIij profiles...");
					timer = countLCIij(parameters.isCountLCIij(), resultFolder1, organism.getName(), timer);
					//*/
					
					//Finish the current organism folder and move to the next one.
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
		
		//Search for orthologs
		if (parameters.isSearchOrthologs()) {
			LOG.debug("Searching orthologs...");
			Orthologs.search(files, resultFolder, parameters, console);
			//Measure the time spent searching for orthologs.
			timer = Auxiliary_Class.time("10th stage - Orthologs", timer);
		}
		
		//Determination of the duration of the entire program
		Auxiliary_Class.time("All", timeStart);
		
		//Completion of the programme
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
		//Count codon frequencies in genes with flanks retained.
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
		//Automatic selection of flanks (if specified)
		Auxiliary_Class.flangSizeAutocounting(parameters, organism, number2, freq);
		//Sort by ascending gene number to restore the original order.
		Sorter.sortNup(organism);
		//Counting the frequencies of codons in genes without taking into account codons at edges
		freq.countFreq(organism.getGenes(), parameters.getUnUseCodonsL1(), parameters.getUnUseCodonsR1());
		//Determination of the time of calculation of codon frequencies and flank sizes (if given)
	}
	
	private static long countLCIij(final boolean countLCIij, final File resultFolder1, final String organismName, final long timer) {
		long timerLocal = timer;
		//Building LCIij profiles
		if (countLCIij) {
			//Source file with LCIij
			File fileInput1 = new File(resultFolder1, organismName + "_lciij.txt");
			//Source file from EEI
			File fileInput2 = new File(resultFolder1, organismName + "_number_eei.txt");
			//File with resultant profiles
			File fileOutput = new File(resultFolder1, organismName + "_lciij_profile.xlsx");
			
			//Profile calculations
			LCIijDataProcessing.readLCIij(fileInput1, fileInput2, fileOutput);
			//Calculation time of LCIij_profile
			timerLocal = Auxiliary_Class.time("8th stage - LCIij_profile", timerLocal);
			
			//Determining the time of rewriting InOut
			timerLocal = Auxiliary_Class.time("9th stage - InOut", timerLocal);
		}
		return timerLocal;
	}
		
}
