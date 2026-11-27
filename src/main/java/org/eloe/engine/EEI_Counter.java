package org.eloe.engine;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JTextArea;

import org.apache.log4j.Logger;

import org.eloe.engine.ta.Ta_Counter;
import org.eloe.engine.te.Te_Counter;
import org.eloe.organism.Organism;
import org.eloe.other.Auxiliary_Class;
import org.eloe.other.Sorter;
import org.eloe.parameters.PEEI;

public class EEI_Counter {
	
	private static final Logger LOG = Logger.getLogger(EEI_Counter.class);
	
	private int number2; //Количество генов в обучающей выборке
	private PEEI parameters;
	private File file_ribPoz; //Файл для записи позиций рибосомных генов
	
	//Конструктор
	public EEI_Counter(int number2_1, PEEI parameters, File file_ribPoz_1) {
		this.number2 = number2_1;
		this.parameters = parameters;
		this.file_ribPoz = file_ribPoz_1;
	}
	
	public int[][] countEEI(Organism organism, JTextArea console) {
		//Массив для индексов (5 - типов индекса, 2 - параметры M и R)
		int[][] mr = new int[5][2];
		
		//EEI(1) - кодонный состав
		mr[0] = count_EEI1(organism);
		String say = "M1=" + mr[0][0] + "\tR1=" + mr[0][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(2) - шпильки без энергии (false)
		mr[1] = count_EEI23(organism, false);
		say = "M2=" + mr[1][0] + "\tR2=" + mr[1][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(3) - шпильки с энергией (true)
		mr[2] = count_EEI23(organism, true);
		say = "M3=" + mr[2][0] + "\tR3=" + mr[2][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(4) - кодонный состав + шпильки без энергии (false)
		mr[3] = count_EEI45(organism, false);
		say = "M4=" + mr[3][0] + "\tR4=" + mr[3][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(5) - кодонный состав + шпильки с энергией (true)
		mr[4] = count_EEI45(organism, true);
		say = "M5=" + mr[4][0] + "\tR5=" + mr[4][1];
		console.append("\n" + say);
		LOG.info(say);
		return mr;
	}
	
	//TODO Расчет EEI1
	private int[] count_EEI1(Organism organism) {
		//Подсчет Ta(i) (не работает без "Подсчета частот кодонов в гене")
		Ta_Counter a = new Ta_Counter(parameters.getMinKc(), parameters.getVe(), number2, 
				parameters.getUnUseCodonsL1(), parameters.getUnUseCodonsR1());
		a.countTaiFirst(organism.getGenes());
		
		//Сортировка списка по убыванию Tai
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
		
		//Подсчет EEI
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
		int djmass[] = Auxiliary_Class.ribpozition(organism.getGenes(), organism.getRibGenesNumber(), file_ribPoz);
		
		//Запись значений EEI1 (=1/Ta)
		for (int i = 0; i < organism.getGenes().size(); i++) {
			organism.getGenes().get(i).putEEI1(organism.getGenes().get(i).getTai());
		}
		
		//Среднее и смещение позиций рибосомных генов
		return djmass;
	}
	
	//TODO Расчет EEI2 и EEI3
	private int[] count_EEI23(Organism organism, boolean energy) {
		//Вычисление Te(i)
		Te_Counter.tei(organism.getGenes(), 0.1, 1.0, energy, parameters.isUseUserParam2(), 
				parameters.getTmin1(), parameters.getTmax1());
		//Сортировка списка по убыванию LCI
		Sorter.sortTedown(organism);
		
		//Без учета энергии (LCI1)
		if (!energy) {			
			//Запись значений EEI2 (=1/Te(LCI1))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI2(organism.getGenes().get(i).getTei());
			}
		}
		//C учетом энергии (LCI2)
		else {
			//Запись значений EEI3 (=1/Te(LCI2))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI3(organism.getGenes().get(i).getTei());
			}
		}
		//Нахождение и запись в файл позиций рибосомных генов
		int djmass[] = Auxiliary_Class.ribpozition(organism.getGenes(), organism.getRibGenesNumber(), file_ribPoz);
		//Среднее и смещение позиций рибосомных генов
		return djmass;
	}
	
	//TODO Расчет EEI4 и EEI5
	private int[] count_EEI45(Organism organism, boolean energy) {
		//Подсчет Ta(i) (не работает без "Подсчета частот кодонов в гене")
		Ta_Counter a = new Ta_Counter(parameters.getMinKc(), parameters.getVe(), number2, 
				parameters.getUnUseCodonsL1(), parameters.getUnUseCodonsR1());
		a.countTaiFirst(organism.getGenes());
		
		//Получаем значения Tmin и Tmax
		double tmin = a.getTmin();
		double tmax = a.getTmax();
		
		//Вычисление Te(i)
		Te_Counter.tei(organism.getGenes(), tmin, tmax, energy, parameters.isUseUserParam2(), 
				parameters.getTmin1(), parameters.getTmax1());
		//Вычисление EEI
		eei(organism);
		
		//Сортировка по убыванию EEI
		Sorter.sorteeidown(organism);
		
		//Порядок генов в выборке
		int number = organism.getGenes().size();
		int numberd = number-number2;
		Set<Integer> set1 = new HashSet<Integer>(number2);
		//Запись положения генов обучающей выборки
		for (int i = 0; i < number2; i++) {
			set1.add(organism.getGenes().get(numberd+i).getNumber());
		}
		boolean repeatflag = false;
		int repeatcounter = 0;
		
		//Подсчет EEI
		do {
			repeatflag = false;
			//Счетчик повторов цикла
			repeatcounter++;
			//Подсчет Tai (не работает без "Подсчета частот кодонов в гене")
			a.countTai(organism.getGenes());
			tmin = a.getTmin();
			tmax = a.getTmax();
			
			//Вычисление Te(i)
			Te_Counter.tei(organism.getGenes(), tmin, tmax, energy, parameters.isUseUserParam2(), 
					parameters.getTmin1(), parameters.getTmax1());
			//Вычисление EEI
			eei(organism);
			
			//Сортировка по EEI по убыванию
			Sorter.sorteeidown(organism);
			Set<Integer> set2 = new HashSet<Integer>(number2);
			//Проверка положений генов обучающей выборки
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
		int djmass[] = Auxiliary_Class.ribpozition(organism.getGenes(), organism.getRibGenesNumber(), file_ribPoz);
		
		//Запись значений EEI
		if (!energy) {
			//Запись значений EEI4 (=2/(Ta+Te))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI4(organism.getGenes().get(i).getEei());
			}
		}
		else {
			//Запись значений EEI5 (=2/(Ta+Te))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI5(organism.getGenes().get(i).getEei());
			}
		}
		
		//Среднее и смещение позиций рибосомных генов
		return djmass;
	}
	
	//TODO Вычисление EEI = Ta + Te
	private void eei(Organism organism) {
		for (int i = 0; i < organism.getGenes().size(); i++) {
			//Кодонный состав + шпильки
			double eei = organism.getGenes().get(i).getTai()+organism.getGenes().get(i).getTei();
			organism.getGenes().get(i).setEei(eei);
		}
	}
}
