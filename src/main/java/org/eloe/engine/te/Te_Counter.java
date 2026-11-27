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
	
	//Computation of P_i
	public static void countTei(List<Gene> l, double m, double sigma, boolean x) {
		//m is the mean value of LCI, sigma is the mean deviation of LCI from the mean
		double k = m/sigma/sigma; //Parameter k
		double n = m*k; //Parameter n
		
		if (n > 168) {
			System.out.println("Danger! Gamma function may lead to overflow of type Double! (x32 bit pricessors)");
			LOG.fatal("Warning! Gamma function may lead to overflow of type Double! (x32 bit pricessors)");
			return;
		}
		//Gamma function G(n+1)
		double g = GammaFunction.calculate(n+1);
		
		//When change function, change integration limits
		Function f = new SigmaFunction2(n);
		Integral integral = new GaussChristoffel7();
		//Accuracy of calculation of integral
		double epsilon = 1e-8;
		
		//If hairpin energy is excluded
		if (!x) {
			//Searching the array of genes and counting the integral for each
			for (int i = 0; i < l.size(); i++){
				double pi1 = integral.calculate(0, k*l.get(i).getLci1(), f, epsilon)/g;
				l.get(i).setPi1(pi1);
			}
		}
		//If hairpin energy is included
		else {
			//Searching the array of genes and counting the integral for each
			for (int i = 0; i < l.size(); i++){
				double pi2 = integral.calculate(0, k*l.get(i).getLci2(), f, epsilon)/g;
				l.get(i).setPi2(pi2);
			}
		}
	}

	//Calculation of Te(i)
	public static void tei(List<Gene> l, double tmin, double tmax, boolean x, boolean useUserParam2, double tmin1, double tmax1) {
		if (useUserParam2) {
			tmin = tmin1;
			tmax = tmax1;
		}
		//Without energy
		if (!x) {
			for (int i = 0; i < l.size(); i++){
				//Probability calculated from the mean and deviation of LCI1
				double p = l.get(i).getPi1();
				//Calculation of Te(i)
				double tei = tmin*(1-p)+tmax*p;
				l.get(i).setTei(tei);
			}
		}
		//With energy
		else {
			for (int i = 0; i < l.size(); i++){
				//Probability calculated from the mean and deviation of LCI2
				double p = l.get(i).getPi2();
				//Calculation of Te(i)
				double tei = tmin*(1-p)+tmax*p;
				l.get(i).setTei(tei);
			}
		}
	}
}
