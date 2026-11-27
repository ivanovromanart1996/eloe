package org.eloe.engine.ta;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;
import org.eloe.other.Sorter;

public class M1_Counter {
	
	private static final Logger LOG = Logger.getLogger(M1_Counter.class);
	
	public static int[] count_m1(Organism organism, int number2, double minKc, double ve, int unUseCodonsL, int unUseCodonsR) {
		//Подсчет Tai (не работает без "Подсчета частот кодонов в гене")
		Ta_Counter a = new Ta_Counter(minKc, ve, number2, unUseCodonsL, unUseCodonsR);
		a.countTaiFirst(organism.getGenes());
		
		//Сортировка списка по Tai
		Sorter.sortTadown(organism);
		
		//Порядок генов в обучающей выборке
		int number = organism.getGenes().size();
		int numberd = number-number2;
		Set<Integer> set1 = new HashSet<Integer>(number2);
		
		//Запись положения генов обучающей выборки
		for (int i = 0; i < number2; i++) {
			set1.add(organism.getGenes().get(numberd+i).getNumber());
		}
		boolean repeatflag = false;
		int repeatcounter = 0;
		
		//Подсчет Tai
		do {
			repeatflag = false;
			//Счетчик количества циклов
			repeatcounter++;
			//Подсчет Tai (не работает без "Подсчета частот кодонов в гене")
			a.countTai(organism.getGenes());
			
			//Сортировка списка по убыванию Tai
			Sorter.sortTadown(organism);
			Set<Integer> set2 = new HashSet<Integer>(number2);
			//Проверка положений контрольных генов
			for (int i = 0; i < number2; i++) {
				set2.add(organism.getGenes().get(numberd+i).getNumber());
			}
			
			//Если обучающие выборки не совпадают, повтор цикла
			if (!set1.containsAll(set2)) {
				set1 = set2;
				repeatflag = true;
			}
		} while (repeatflag && (repeatcounter < 30));
		
		//В случае зацикливания
		if (repeatcounter >= 30) {
			LOG.warn("Warning! Many cicles! RepeatCounter >= 30 (have to be <30)");
		}
		
		//Нахождение и запись в файл позиций рибосомных генов
		int djmass[] = ribpozition1(organism.getGenes(), organism.getRibGenesNumber());
		
		//Среднее и смещение позиций рибосомных генов
		return djmass;
	}
	
	//Позиция рибосомальных генов
	private static int[] ribpozition1(List<Gene> l, int ribgenes) {
		int djmass[] = {0, 0};
		if (ribgenes != 0) {
			int ribPozMass[] = new int[ribgenes];
			int j = 0;
			int ng = 0;
			for (int i = 0; i < l.size(); i++) {
				if (l.get(i).isRibosomal()) {
					ribPozMass[j] = i+1;
					j++;
					ng += (i+1);
				}
			}
			double m = (double)ng/ribgenes;
			double r = 0;
			for (int i = 0; i < ribgenes; i++) {
				r += Math.pow((ribPozMass[i]-m), 2);
			}
			r = Math.sqrt(r/ribgenes);
			int dj = (int)Math.round(100*(2*m-l.size()-1)/(l.size()-ribgenes));
			int djr = (int)Math.round(100*(2*r-ribgenes-1)/(l.size()-ribgenes));
			/*
			int dj = (int)Math.round(100*(2*(m-1)/(l.size()-1)-1));
			int djr = (int)Math.round(200*(r-1)/(l.size()-1));
			*/
			djmass[0] = dj;
			djmass[1] = djr;
		}
		return djmass;
	}
}
