package org.eloe.engine.te.lci;

import java.util.List;

import org.eloe.engine.te.lci.lci1.Lci1Calculator;
import org.eloe.organism.Gene;
import org.eloe.parameters.PLCI1;

// Програма для подсчета индекса локальной комплементарности
public class LCI1Counter {
	
	//Среднее значение LCI
	private double m1;
	//Среднее отклонение LCI
	private double sigma1;
	
	private ILciCalculator lci1Calculator;
	
	//Конструктор
	public LCI1Counter(PLCI1 parameters) {
		lci1Calculator = new Lci1Calculator(parameters);
	}
	
	//TODO Рассчет LCI1
	public final void countLCI(List<Gene> genes) {
		
		double sum = 0;
		double sum2 = 0;
		double lci1 = 0;

		for (Gene gene : genes) {
			lci1 = lci1Calculator.calculate(gene);
			gene.setLci1(lci1);
			sum += lci1;
			sum2 += lci1*lci1;
		}
		
		double average = sum/genes.size();
		double disp = Math.sqrt(sum2/genes.size()-average*average);
		
		//фортран
		this.m1 = average+2*disp;
		this.sigma1 = 4*disp;
		
		//статьи
		//this.m1 = average;
		//this.sigma1 = forSigma1*disp;
	}
	
	//TODO Возвращает среднее значение LCI
	public final double getM() {
		return this.m1;
	}
	
	//TODO Возвращает среднее отклонение LCI
	public final double getSigma() {
		return this.sigma1;
	}
}
