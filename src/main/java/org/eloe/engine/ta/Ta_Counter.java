package org.eloe.engine.ta;

import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;

public class Ta_Counter {
	
	private static final Logger LOG = Logger.getLogger(Ta_Counter.class);
	
	private static final int SUMFREQSIZE = 64; //Размер массива summFreq (количество триплетов)
	private double atrmin = 0.1;
	private double atrmax = 1.0;
	private double atrmn1 = 0.1;
	private double tmin = 0.1;
	private double tmax = 1.0;
	//private double minKc; //Минимальное значение Kc
	//private double ve; //Нормировочный множитель
	private int number2; //Количество генов в обучающей выборке
	private int unUseCodonsL; //Количество неиспользуемых в расчетах кодонов на левом фланге
	private int unUseCodonsR; //Количество неиспользуемых в расчетах кодонов на правом фланге
	
	//Конструктор
	public Ta_Counter(double minKc_1, double ve_1, int number2_1, int unUseCodonsL_1, int unUseCodonsR_1) {
		//this.minKc = minKc_1;
		//this.ve = ve_1;
		this.number2 = number2_1;
		this.unUseCodonsL = unUseCodonsL_1;
		this.unUseCodonsR = unUseCodonsR_1;
	}
	
	//TODO Первый расчет Ta(i)
	public void countTaiFirst(List<Gene> l) {
		//Массив корней из суммы частот кодонов
		double kc[] = new double[SUMFREQSIZE];
		//Создание массива случайных чисел номеров генов для расчета
		Set<Integer> setOfRandomInt = setOfRandomInt(l.size());
		//Суммирование частот одинаковых кодонов по генам из выборки
		for (Integer i : setOfRandomInt) {
			int triplFreq0[] = l.get(i).getCodonComposition0();
			for (int j = 0; j < SUMFREQSIZE; j++) {
				kc[j] += triplFreq0[j];
			}
		}
		//Основной метод, рассчитывающий tai
		tai(l, kc);
	}

	//TODO Расчет Ta(i)
	public void countTai(List<Gene> l) {
		//Массив корней из суммы частот кодонов
		double kc[] = new double[SUMFREQSIZE];
		//Суммирование частот одинаковых кодонов по генам из выборки
		for (int i = l.size()-1, max = l.size()-number2-1; i > max; i--) {
			int triplFreq0[] = l.get(i).getCodonComposition0();
			for (int j = 0; j < SUMFREQSIZE; j++) {
				kc[j] += triplFreq0[j];
			}
		}
		//Основной метод, рассчитывающий tai
		tai(l, kc);
	}
	
	//TODO Вернуть Tmin
	public double getTmin(){
		return this.tmin;
	}
	
	//TODO Вернуть Tmax
	public double getTmax(){
		return this.tmax;
	}
	
	private Set<Integer> setOfRandomInt(int l_size) {
		Random rand = new Random(System.currentTimeMillis());
		Set<Integer> setOfIntegers = new HashSet<Integer>();
		do {
			//Добавление в набор случайного номера гена
			setOfIntegers.add(rand.nextInt(l_size));
		} while (setOfIntegers.size() < number2);
		return setOfIntegers;
	}
	
	//TODO Основной метод, расчитывающий tai
	private void tai(List<Gene> l, double kc[]) {
		double saai = 0;
		double kcmin = Double.MAX_VALUE;
		for (int i = 0; i < SUMFREQSIZE; i++) {
			//Извлечение корней
			kc[i] = (double)Math.sqrt(kc[i]);
			if (kcmin > kc[i]) {
				kcmin = kc[i];
			}
			//Сумма всех элементов массива
			saai += kc[i];
		}
		//Среднее значение суммы корней, деленной на корень
		double srednt = 0;
		for (int i = 0; i < SUMFREQSIZE; i++) {
			kc[i] = (kc[i] == 0) ? 64 : saai/kc[i]; //64-среднее значение для kc (=saai/(saai/64)))
			srednt += kc[i];
		}
		//Среднее (/64)
		srednt /= SUMFREQSIZE;
		
		double aaimin = Double.MAX_VALUE;
		for (int i = 0; i < SUMFREQSIZE; i++) {
			kc[i] /= srednt;
			if (aaimin > kc[i]) {
				aaimin = kc[i];
			}
		}
		for (int i = 0; i < SUMFREQSIZE; i++) {
			kc[i] = atrmin*kc[i]/aaimin;
			if (kc[i] > atrmax) {
				kc[i] = atrmax;
			}
			if (kc[i] < atrmn1) {
				kc[i] = atrmn1;
			}
		}
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		for (int i = 0; i < l.size(); i++) {
			double ta = 0;
			int triplFreq[] = l.get(i).getCodonComposition();
			int codons = l.get(i).getCodonsNumber()-(unUseCodonsL+unUseCodonsR);
			if (codons > 0) {
				for (int j = 0; j < SUMFREQSIZE; j++) {
					ta += triplFreq[j]*kc[j];
				}
				//Вычисление Tai
				ta /= codons;
			}
			l.get(i).setTai(ta);
			if (min > ta) {
				min = ta;
			}
			if (max < ta) {
				max = ta;
			}
		}
		double xaix1 = 0;
		double xaix2 = 0;
		if (max > min) {
			xaix1 = (this.tmax - this.tmin)/(max - min);
			xaix2 = (this.tmin*max-this.tmax*min)/(max-min);
		}
		else {
			LOG.error("Warning! max <= min. unUseCodonsL = "+unUseCodonsL+", unUseCodonsR = "+unUseCodonsR);
		}
		//Вычисление и запись Tai
		for (int i = 0; i < l.size(); i++) {
			l.get(i).setTai(l.get(i).getTai()*xaix1+xaix2);
		}
	}
}
