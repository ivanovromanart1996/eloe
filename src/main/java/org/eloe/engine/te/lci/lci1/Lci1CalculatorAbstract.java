package org.eloe.engine.te.lci.lci1;

import org.apache.log4j.Logger;

import org.eloe.engine.te.lci.ILciCalculator;
import org.eloe.organism.IGeneWithFlanks;
import org.eloe.parameters.PLCI1;

public abstract class Lci1CalculatorAbstract implements ILciCalculator {
	
	private static final Logger LOG = Logger.getLogger(Lci1CalculatorAbstract.class);
	
	private PLCI1 parameters;
	private int distance;
	private int rightFlank;
	
	public Lci1CalculatorAbstract(PLCI1 parameters) {
		this.parameters = parameters;
		distance = parameters.getLmax1() - parameters.getLmin1();
		rightFlank = parameters.getSmax1() + parameters.getLmax1();
	}

	@Override
	public double calculate(IGeneWithFlanks gene) {
		
		if (rightFlank > gene.getRightFlankLength()) {
			LOG.info("The right flank length of gene "+gene.getName()+" is not enough for LCI1 calculation.");
			return 0;
		}
		
		//Длина исследуемой части гена (3 - стоп-кодон)
		int lenght = gene.getLength()-3+rightFlank;
		//Граница левого края скользящей рамки
		int edge = lenght-2*parameters.getSmax1()-parameters.getLmax1()+1;
		if (edge < 0) {
			LOG.info("The length of gene "+gene.getName()+" with flanks is not enough for LCI1 calculation.");
			return 0;
		}
		
		char[] na = new char[lenght];
		//Достаем ген без флангов + h нуклеотидов на 3' конце
		System.arraycopy(gene.getNaWithFlanks(), gene.getLeftFlankLength(), na, 0, lenght);
		
		//Счетчик повторов
		int counter = 0;
		for (int j = 0; j < edge; j++) {
			for (int s = parameters.getSmin1(); s <= parameters.getSmax1(); s++) {
				
				//Строка для повтора
				char[] char1 = new char[s];
				//Копирование последовательности
				System.arraycopy(na, j, char1, 0, s);
				//Строка для комплементарной последовательности
				String con1 = String.valueOf(reverseComplement(char1));
				
				//Длина участка поиска
				int lenght1 = s+distance;
				//Участок поиска
				char[] char2 = new char[lenght1];
				//Достаем участок поиска
				System.arraycopy(na, j+s+parameters.getLmin1(), char2, 0, lenght1);
				//Участок поиска
				String con2 = String.valueOf(char2);
				
				int poz = 0;
				//Поиск повторов (комплементарных участков)
				int k = con2.indexOf(con1, poz);
				while (k != -1) {
					//Подсчет повторов
					counter++;
					poz = k + 1;
					k = con2.indexOf(con1, poz);
				}
			}
		}
		
		return (double) counter/(gene.getLength()-3);
	}

	abstract protected char[] reverseComplement(char[] char1);
}
