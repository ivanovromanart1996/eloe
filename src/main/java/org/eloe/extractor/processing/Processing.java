package org.eloe.extractor.processing;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;
import org.eloe.parameters.PProcessing;

//Класс для удаления из списка "плохих" генов
public abstract class Processing {
	
	private static final Logger LOG = Logger.getLogger(Processing.class);
	
	public static void processing (Organism organism, PProcessing parameters) {
		int badGenes = 0;
		int ribGenes = 0;
		
		//Удаление генов со значками "n"
		for (int i = 0; i < organism.getGenes().size(); i++) {
			if (organism.getGenes().get(i).isFlagUse()) {
				if (organism.getGenes().get(i).getNAString().contains("n")) {
					//Отметка гена, как неиспользуемого
					organism.getGenes().get(i).setFlagUse(false);
					if (organism.getGenes().get(i).isRibosomal()) {
						ribGenes++;
						organism.delRibGene();
					}
					badGenes++;
				}
			}
		}
		
		if (parameters.isCheckStart() || parameters.isDeleteGwBC() || parameters.isCheckStop()) {
			//Заполнение карты стартовыми кодонами
			Set<String> setStart = new HashSet<String>(parameters.getStartCodons().length);
			for (int j = 0; j < parameters.getStartCodons().length; j++) {
				setStart.add(parameters.getStartCodons()[j].trim());
			}
			//Заполнение карты "плохими" кодонами
			Set<String> setMiddle = new HashSet<String>(parameters.getMidlCodons().length);
			for (int j = 0; j < parameters.getMidlCodons().length; j++) {
				setMiddle.add(parameters.getMidlCodons()[j].trim());
			}
			//Заполнение карты стоп-кодонами
			Set<String> setStop = new HashSet<String>(parameters.getStopCodons().length);
			for (int j = 0; j < parameters.getStopCodons().length; j++) {
				setStop.add(parameters.getStopCodons()[j].trim());
			}
			
			for (int i = 0; i < organism.getGenes().size(); i++) {
				if (organism.getGenes().get(i).isFlagUse()) {
					//Длина гена без флангов
					int lenght = organism.getGenes().get(i).getLength();
					//Для последовательности гена без флангов
					char gene[] = organism.getGenes().get(i).getNA();
					
					char c[] = new char[3];
					String codon = "";
					
					if (parameters.isCheckStart()) {
						System.arraycopy(gene, 0, c, 0, 3);
						codon = String.valueOf(c);
						//Проверка наличия в начале гена стартового кодона
						if (!setStart.contains(codon)) {
							//Отметка гена, как неиспользуемого
							organism.getGenes().get(i).setFlagUse(false);
							if (organism.getGenes().get(i).isRibosomal()) {
								ribGenes++;
								organism.delRibGene();
							}
							badGenes++;
						}
					}
					
					if (parameters.isDeleteGwBC()) {
						//Число кодонов в гене без последнего
						int codons = organism.getGenes().get(i).getCodonsNumber()-1;
						for (int j = 1; j < codons; j++) {
							System.arraycopy(gene, 3*j, c, 0, 3);
							codon = String.valueOf(c);
							//Проверка наличия в середине гена плохого кодона
							if (setMiddle.contains(codon)) {
								//Отметка гена, как неиспользуемого
								organism.getGenes().get(i).setFlagUse(false);
								if (organism.getGenes().get(i).isRibosomal()) {
									ribGenes++;
									organism.delRibGene();
								}
								badGenes++;
								break;
							}
						}
					}
					
					if (parameters.isCheckStop()) {
						System.arraycopy(gene, lenght-3, c, 0, 3);
						codon = String.valueOf(c);
						//Проверка наличия на конце гена стоп-кодона
						if (!setStop.contains(codon)) {
							//Отметка гена, как неиспользуемого
							organism.getGenes().get(i).setFlagUse(false);
							if (organism.getGenes().get(i).isRibosomal()) {
								ribGenes++;
								organism.delRibGene();
							}
							badGenes++;
						}
					}
				}
			}
		}
		
		//Удаление из списка всех помеченных генов
		int counterDeletGenes = 0;
		for (int i = 0; i < organism.getGenes().size(); i++) {
			if (!organism.getGenes().get(i).isFlagUse()) {
				organism.getGenes().remove(i);
				counterDeletGenes++;
				i--;
			}
		}
		((ArrayList<Gene>) organism.getGenes()).trimToSize();
		
		int pseudogenes = organism.getPseudoGenesNumber();
		if (parameters.isUsePseudoGenes()) {
			pseudogenes = 0;
		}
		//Сообщение о количестве удаленных рибосомных генов
		//Сообщение о том, что короткие гены были удалены
		//Сообщение о том, что псевдо гены были удалены
		//Сообщение о том, что плохие гены были удалены
		//Сообщение об общем количестве удаленных генов
		//Сообщение о количестве генов, участвующих в расчетах
		LOG.info(ribGenes+" ribosomal genes were deleted!\n"
				+organism.getShortGenesNumber()+" short genes were deleted!\n"
				+pseudogenes+" pseudogenes were deleted!\n"
				+badGenes+" bad genes were deleted!\n"
				+ "In all "+counterDeletGenes+" genes were deleted!\n"
				+ "Number of genes used in counting: "+organism.getGenes().size());
	}
}
