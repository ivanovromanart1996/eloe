package org.eloe.engine.te.lci;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eloe.engine.te.lci.lci2.Lci2Calculator;
import org.eloe.organism.Gene;
import org.eloe.parameters.PLCI2;

// Програма для подсчета индекса локальной комплементарности

public class LCI2_Counter {
	
	private static final String ENERGY_TXT = "lci2double.txt";
	
	//Среднее значение LCI
	private double m2;
	//Среднее отклонение LCI
	private double sigma2;
	
	private ILciCalculator lci2Calculator;
	
	//Конструктор
	public LCI2_Counter(PLCI2 parameters) throws IOException, NumberFormatException {
		
		//Карта для записи энергии динуклеотидов
		Map<String, Double> map = new HashMap<String, Double>();
		StringBuilder loopsenergy = new StringBuilder();
		BufferedReader reader = null;
		try {
			String line = new String();
			//Чтение файла с энергиями динуклеотидов и петель
			reader = new BufferedReader(new InputStreamReader(
					LCI2_Counter.class.getClassLoader().getResourceAsStream(ENERGY_TXT)));
			while (!(line = reader.readLine()).startsWith("//")){
				String arr[] = line.split(" - ");
				map.put(arr[0], Double.parseDouble(arr[1]));
			}
			
			while (!((line = reader.readLine()).startsWith("////"))) {
				loopsenergy.append(line);
			}
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		String loopsenergy1[] = loopsenergy.toString().split(";");
		//Массив для записи энергии петель
		double loopenergy[] = new double[loopsenergy1.length];
		for (int i = 0; i < loopsenergy1.length; i++) {
			loopenergy[i] = Double.parseDouble(loopsenergy1[i].trim());
		}
		
		lci2Calculator = new Lci2Calculator(parameters, map, loopenergy);
	}
	
	//TODO Рассчет LCI2
	public void countLCI(List<Gene> genes) {
		
		double sum = 0;
		double sum2 = 0;
		double lci2 = 0;

		for (Gene gene : genes) {
			lci2 = lci2Calculator.calculate(gene);
			gene.setLci2(lci2);
			sum += lci2;
			sum2 += lci2*lci2;
		}
		
		double average = sum/genes.size();
		double disp = Math.sqrt(sum2/genes.size()-average*average);
		
		//фортран
		this.m2 = average+disp;
		this.sigma2 = 2*disp;
		
		//статьи
		//this.m2 = average;
		//this.sigma2 = forSigma2*disp;
	}
	
	//TODO Возвращает среднее значение LCI
	public double getM() {
		return this.m2;
	}
	
	//TODO Возвращает среднее отклонение LCI
	public double getSigma() {
		return this.sigma2;
	}
}
