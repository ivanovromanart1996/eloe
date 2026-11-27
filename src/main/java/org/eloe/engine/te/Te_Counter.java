package org.eloe.engine.te;

import java.util.List;

import org.apache.log4j.Logger;

import org.eloe.math.Function;
import org.eloe.math.GammaFunction;
import org.eloe.math.GaussChristoffel7;
import org.eloe.math.Integral;
import org.eloe.math.SigmaFunction2;
import org.eloe.organism.Gene;


public abstract class Te_Counter {
	
	private static final Logger LOG = Logger.getLogger(Te_Counter.class);
	
	//Вычисление P_i
	public static void countTei(List<Gene> l, double m, double sigma, boolean x) {
		//m - среднее значение LCI, sigma - среднее отклонение LCI от среднего
		double k = m/sigma/sigma; //Параметр k
		double n = m*k; //Параметр n
		
		if (n > 168) {
			System.out.println("Danger! Gamma function may lead to overflow of type Double! (x32 bit pricessors)");
			LOG.fatal("Warning! Gamma function may lead to overflow of type Double! (x32 bit pricessors)");
			return;
		}
		//Гамма функция G(n+1)
		double g = GammaFunction.calculate(n+1);
		
		//When change function, change integration limits
		Function f = new SigmaFunction2(n);
		Integral integral = new GaussChristoffel7();
		//Tочность вычисления интеграла
		double epsilon = 1e-8;
		
		//Если без учета энергии шпильки
		if (!x) {
			//Перебор массива генов и подсчет интеграла для каждого
			for (int i = 0; i < l.size(); i++){
				double pi1 = integral.calculate(0, k*l.get(i).getLci1(), f, epsilon)/g;
				l.get(i).setPi1(pi1);
			}
		}
		//Если с учетом энергии шпильки
		else {
			//Перебор массива генов и подсчет интеграла для каждого
			for (int i = 0; i < l.size(); i++){
				double pi2 = integral.calculate(0, k*l.get(i).getLci2(), f, epsilon)/g;
				l.get(i).setPi2(pi2);
			}
		}
	}

	//Вычисление Te(i)
	public static void tei(List<Gene> l, double tmin, double tmax, boolean x, boolean useUserParam2, double tmin1, double tmax1) {
		if (useUserParam2) {
			tmin = tmin1;
			tmax = tmax1;
		}
		//Без энергии
		if (!x) {
			for (int i = 0; i < l.size(); i++){
				//Вероятность, вычисленная из среднего и отклонения LCI1
				double p = l.get(i).getPi1();
				//Вычисление Te(i)
				double tei = tmin*(1-p)+tmax*p;
				l.get(i).setTei(tei);
			}
		}
		//С энергией
		else {
			for (int i = 0; i < l.size(); i++){
				//Вероятность, вычисленная из среднего и отклонения LCI2
				double p = l.get(i).getPi2();
				//Вычисление Te(i)
				double tei = tmin*(1-p)+tmax*p;
				l.get(i).setTei(tei);
			}
		}
	}
}
