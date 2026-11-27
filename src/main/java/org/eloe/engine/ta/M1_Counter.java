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
		//Tai counting (does not work without "Counting codon frequencies in a gene")
		Ta_Counter a = new Ta_Counter(minKc, ve, number2, unUseCodonsL, unUseCodonsR);
		a.countTaiFirst(organism.getGenes());
		
		//Sorting the list by Tai
		Sorter.sortTadown(organism);
		
		//The order of genes in the training sample
		int number = organism.getGenes().size();
		int numberd = number-number2;
		Set<Integer> set1 = new HashSet<Integer>(number2);
		
		//Recording the position of the genes of the training sample
		for (int i = 0; i < number2; i++) {
			set1.add(organism.getGenes().get(numberd+i).getNumber());
		}
		boolean repeatflag = false;
		int repeatcounter = 0;
		
		//Tai counting
		do {
			repeatflag = false;
			//Cycle counter
			repeatcounter++;
			//Tai counting (does not work without "Counting codon frequencies in a gene")
			a.countTai(organism.getGenes());
			
			//Sorting the list by descending Tai
			Sorter.sortTadown(organism);
			Set<Integer> set2 = new HashSet<Integer>(number2);
			//Checking the positions of control genes
			for (int i = 0; i < number2; i++) {
				set2.add(organism.getGenes().get(numberd+i).getNumber());
			}
			
			//If the training samples do not match, repeat the cycle
			if (!set1.containsAll(set2)) {
				set1 = set2;
				repeatflag = true;
			}
		} while (repeatflag && (repeatcounter < 30));
		
		//If the iteration cycles
		if (repeatcounter >= 30) {
			LOG.warn("Warning! Many cicles! RepeatCounter >= 30 (have to be <30)");
		}
		
		//Finding and writing positions of ribosomal genes
		int djmass[] = ribpozition1(organism.getGenes(), organism.getRibGenesNumber());
		
		//Mean and displacement of ribosomal gene positions
		return djmass;
	}
	
	//Position of ribosomal genes
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
