package org.eloe.engine.te.lci.lci2;

import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import org.eloe.engine.te.lci.ILciCalculator;
import org.eloe.organism.IGeneWithFlanks;
import org.eloe.parameters.PLCI2;

public abstract class Lci2CalculatorAbstract implements ILciCalculator {
	
private static final Logger LOG = Logger.getLogger(Lci2CalculatorAbstract.class);
	
	private PLCI2 parameters;
	private int distance;
	private int rightFlank;
	//Карта для записи энергии динуклеотидов
	private Map<String, Double> map = new HashMap<String, Double>();
	//Массив с энергиями петель
	private double loopenergy[];
	
	public Lci2CalculatorAbstract(PLCI2 parameters, Map<String, Double> map, double[] loopenergy) {
		this.parameters = parameters;
		distance = parameters.getLmax2()-parameters.getLmin2();
		rightFlank = parameters.getSmax2()+parameters.getLmax2();
		this.map = map;
		this.loopenergy = loopenergy;
	}

	@Override
	public double calculate(IGeneWithFlanks gene) {
		
		if (rightFlank > gene.getRightFlankLength()) {
			LOG.info("The right flank length of gene "+gene.getName()+" is not enough for LCI2 calculation.");
			return 0;
		}
		
		//Длина исследуемой части гена (3 - стоп-кодон)
		int lenght = gene.getLength()-3+rightFlank;
		//Граница левого края скользящей рамки
		int edge = lenght-2*parameters.getSmax2()-parameters.getLmax2()+1;
		if (edge < 0) {
			LOG.info("The length of gene "+gene.getName()+" with flanks is not enough for LCI2 calculation.");
			return 0;
		}
		
		char[] na = new char[lenght];
		//Достаем ген без флангов + h нуклеотидов на 3' конце
		System.arraycopy(gene.getNaWithFlanks(), gene.getLeftFlankLength(), na, 0, lenght);
		
		//Полная энергия взаимодействия
		double energyFull = 0;
		for (int j = 0; j < edge; j++) {
			//Энергия шпилек
			double energy = 0;
			for (int s = parameters.getSmin2(); s <= parameters.getSmax2(); s++) {
				
				//Строка для повтора
				char[] char1 = new char[s];
				//Копирование последовательности
				System.arraycopy(na, j, char1, 0, s);
				
				//Энергия выбранного повтора
				double energy_char1 = 0;
				
				for (int m = 0; m < s-1; m++) {
					String dinukleotid = char1[m]+""+char1[m+1];
					
					if (map.containsKey(dinukleotid)) {
						//Энергия взаимодействия повторов
						energy_char1 += map.get(dinukleotid);
					}
				}
				
				//Строка для комплементарной последовательности
				String con1 = String.valueOf(reverseComplement(char1));
				
				//Длина участка поиска
				int lenght1 = s+distance;
				//Участок поиска
				char[] char2 = new char[lenght1];
				//Достаем участок поиска
				System.arraycopy(na, j+s+parameters.getLmin2(), char2, 0, lenght1);
				//Участок поиска
				String con2 = String.valueOf(char2);
				
				int poz = 0;
				//Энергия петель
				double loop1 = 0;
				
				//Поиск повторов (комплементарных участков)
				int k = con2.indexOf(con1, poz);
				while (k != -1) {
					int loopsize = parameters.getLmin2()+k-1;
					
					//Энергия петель
					loop1 = loopenergy[(loopsize < loopenergy.length) ? loopsize : (loopenergy.length-1)];
					double e = energy_char1 - loop1;
					
					//Если энергия больше порога существоавания петли, тогда учитываем эту петлю
					if (e >= parameters.getMinloopEnergy()) {
						energy += e;
					}
					poz = k + 1;
					k = con2.indexOf(con1, poz);
				}
			}
			//Полная энергия
			energyFull += energy;
		}
		
		return energyFull/(gene.getLength()-3);
	}
	
	abstract protected char[] reverseComplement(char[] char1);
}
