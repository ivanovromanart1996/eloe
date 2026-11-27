package org.eloe.engine.ta;

import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.log4j.Logger;

import org.eloe.organism.Gene;

public class Ta_Counter {
	
	private static final Logger LOG = Logger.getLogger(Ta_Counter.class);
	
	private static final int SUMFREQSIZE = 64; //Size of the summFreq array (number of triplets)
	private double atrmin = 0.1;
	private double atrmax = 1.0;
	private double atrmn1 = 0.1;
	private double tmin = 0.1;
	private double tmax = 1.0;
	//private double minKc; //Minimum value of Kc
	//private double ve; //Normalization factor
	private int number2; //Number of genes in the training sample
	private int unUseCodonsL; //Number of unused codons on the left flank
	private int unUseCodonsR; //Number of unused codons on the right flank
	
	//Constructor
	public Ta_Counter(double minKc_1, double ve_1, int number2_1, int unUseCodonsL_1, int unUseCodonsR_1) {
		//this.minKc = minKc_1;
		//this.ve = ve_1;
		this.number2 = number2_1;
		this.unUseCodonsL = unUseCodonsL_1;
		this.unUseCodonsR = unUseCodonsR_1;
	}
	
	//TODO First calculation Ta(i)
	public void countTaiFirst(List<Gene> l) {
		//Array of square roots of summed codon frequencies
		double kc[] = new double[SUMFREQSIZE];
		//Build an array of randomly selected gene numbers for the calculation.
		Set<Integer> setOfRandomInt = setOfRandomInt(l.size());
		//Summarizing the frequencies of the same codons by genes from the sample
		for (Integer i : setOfRandomInt) {
			int triplFreq0[] = l.get(i).getCodonComposition0();
			for (int j = 0; j < SUMFREQSIZE; j++) {
				kc[j] += triplFreq0[j];
			}
		}
		//The basic method of calculating tai
		tai(l, kc);
	}

	//TODO Calculation of Ta(i)
	public void countTai(List<Gene> l) {
		//Array of square roots of summed codon frequencies
		double kc[] = new double[SUMFREQSIZE];
		//Summarizing the frequencies of the same codons by genes from the sample
		for (int i = l.size()-1, max = l.size()-number2-1; i > max; i--) {
			int triplFreq0[] = l.get(i).getCodonComposition0();
			for (int j = 0; j < SUMFREQSIZE; j++) {
				kc[j] += triplFreq0[j];
			}
		}
		//The basic method of calculating tai
		tai(l, kc);
	}
	
	//TODO Return Tmin.
	public double getTmin(){
		return this.tmin;
	}
	
	//TODO Return Tmax.
	public double getTmax(){
		return this.tmax;
	}
	
	private Set<Integer> setOfRandomInt(int l_size) {
		Random rand = new Random(System.currentTimeMillis());
		Set<Integer> setOfIntegers = new HashSet<Integer>();
		do {
			//Adding a random gene number to the set
			setOfIntegers.add(rand.nextInt(l_size));
		} while (setOfIntegers.size() < number2);
		return setOfIntegers;
	}
	
	//TODO The basic method for calculating tai
	private void tai(List<Gene> l, double kc[]) {
		double saai = 0;
		double kcmin = Double.MAX_VALUE;
		for (int i = 0; i < SUMFREQSIZE; i++) {
			//Root extraction
			kc[i] = (double)Math.sqrt(kc[i]);
			if (kcmin > kc[i]) {
				kcmin = kc[i];
			}
			//Sum of all elements of the array
			saai += kc[i];
		}
		//Average value of the sum of roots divided by the root
		double srednt = 0;
		for (int i = 0; i < SUMFREQSIZE; i++) {
			kc[i] = (kc[i] == 0) ? 64 : saai/kc[i]; //64-mean value for kc (=saai/(saai/64))
			srednt += kc[i];
		}
		//Average (/64)
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
				//Tai computation
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
		//Calculation and recording of Tai
		for (int i = 0; i < l.size(); i++) {
			l.get(i).setTai(l.get(i).getTai()*xaix1+xaix2);
		}
	}
}
