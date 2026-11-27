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
	
	private int number2; //Number of genes in the training sample
	private PEEI parameters;
	private File file_ribPoz; //File for recording positions of ribosomal genes
	
	//Constructor
	public EEI_Counter(int number2_1, PEEI parameters, File file_ribPoz_1) {
		this.number2 = number2_1;
		this.parameters = parameters;
		this.file_ribPoz = file_ribPoz_1;
	}
	
	public int[][] countEEI(Organism organism, JTextArea console) {
		//Array for indices (5 index types, 2 M and R parameters)
		int[][] mr = new int[5][2];
		
		//EEI(1) - codon composition
		mr[0] = count_EEI1(organism);
		String say = "M1=" + mr[0][0] + "\tR1=" + mr[0][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(2) - hairpins without energy (false)
		mr[1] = count_EEI23(organism, false);
		say = "M2=" + mr[1][0] + "\tR2=" + mr[1][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(3) - hairpins with energy (true)
		mr[2] = count_EEI23(organism, true);
		say = "M3=" + mr[2][0] + "\tR3=" + mr[2][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(4) - codon composition + hairpins without energy (false)
		mr[3] = count_EEI45(organism, false);
		say = "M4=" + mr[3][0] + "\tR4=" + mr[3][1];
		console.append("\n" + say);
		LOG.info(say);
		
		//EEI(5) - codon composition + hairpins with energy (true)
		mr[4] = count_EEI45(organism, true);
		say = "M5=" + mr[4][0] + "\tR5=" + mr[4][1];
		console.append("\n" + say);
		LOG.info(say);
		return mr;
	}
	
	//TODO Calculation EEI1
	private int[] count_EEI1(Organism organism) {
		//Ta(i) counting (does not work without "Counting codon frequencies in a gene")
		Ta_Counter a = new Ta_Counter(parameters.getMinKc(), parameters.getVe(), number2, 
				parameters.getUnUseCodonsL1(), parameters.getUnUseCodonsR1());
		a.countTaiFirst(organism.getGenes());
		
		//Sorting the list by descending Tai
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
		
		//Estimate EEI
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
		int djmass[] = Auxiliary_Class.ribpozition(organism.getGenes(), organism.getRibGenesNumber(), file_ribPoz);
		
		//Record the values of EEI1 (=1/Ta)
		for (int i = 0; i < organism.getGenes().size(); i++) {
			organism.getGenes().get(i).putEEI1(organism.getGenes().get(i).getTai());
		}
		
		//Mean and displacement of ribosomal gene positions
		return djmass;
	}
	
	//TODO Calculation of EEI2 and EEI3
	private int[] count_EEI23(Organism organism, boolean energy) {
		//Calculation of Te(i)
		Te_Counter.tei(organism.getGenes(), 0.1, 1.0, energy, parameters.isUseUserParam2(), 
				parameters.getTmin1(), parameters.getTmax1());
		//Sorting the list by descending LCI
		Sorter.sortTedown(organism);
		
		//Exclusive of energy (LCI1)
		if (!energy) {			
			//Record the values of EEI2 (=1/Te(LCI1))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI2(organism.getGenes().get(i).getTei());
			}
		}
		//Includes energy (LCI2)
		else {
			//Record the values of EEI3 (=1/Te(LCI2))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI3(organism.getGenes().get(i).getTei());
			}
		}
		//Finding and writing positions of ribosomal genes
		int djmass[] = Auxiliary_Class.ribpozition(organism.getGenes(), organism.getRibGenesNumber(), file_ribPoz);
		//Mean and displacement of ribosomal gene positions
		return djmass;
	}
	
	//TODO Calculation of EEI4 and EEI5
	private int[] count_EEI45(Organism organism, boolean energy) {
		//Ta(i) counting (does not work without "Counting codon frequencies in a gene")
		Ta_Counter a = new Ta_Counter(parameters.getMinKc(), parameters.getVe(), number2, 
				parameters.getUnUseCodonsL1(), parameters.getUnUseCodonsR1());
		a.countTaiFirst(organism.getGenes());
		
		//Get the values of Tmin and Tmax
		double tmin = a.getTmin();
		double tmax = a.getTmax();
		
		//Calculation of Te(i)
		Te_Counter.tei(organism.getGenes(), tmin, tmax, energy, parameters.isUseUserParam2(), 
				parameters.getTmin1(), parameters.getTmax1());
		//Calculation of EEI
		eei(organism);
		
		//Sorting by descending EEI
		Sorter.sorteeidown(organism);
		
		//The order of genes in the sample
		int number = organism.getGenes().size();
		int numberd = number-number2;
		Set<Integer> set1 = new HashSet<Integer>(number2);
		//Recording the position of the genes of the training sample
		for (int i = 0; i < number2; i++) {
			set1.add(organism.getGenes().get(numberd+i).getNumber());
		}
		boolean repeatflag = false;
		int repeatcounter = 0;
		
		//Estimate EEI
		do {
			repeatflag = false;
			//Cycle repeater
			repeatcounter++;
			//Tai counting (does not work without "Counting codon frequencies in a gene")
			a.countTai(organism.getGenes());
			tmin = a.getTmin();
			tmax = a.getTmax();
			
			//Calculation of Te(i)
			Te_Counter.tei(organism.getGenes(), tmin, tmax, energy, parameters.isUseUserParam2(), 
					parameters.getTmin1(), parameters.getTmax1());
			//Calculation of EEI
			eei(organism);
			
			//Sorting by EEI by descending
			Sorter.sorteeidown(organism);
			Set<Integer> set2 = new HashSet<Integer>(number2);
			//Checking the positions of the genes of the training sample
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
		int djmass[] = Auxiliary_Class.ribpozition(organism.getGenes(), organism.getRibGenesNumber(), file_ribPoz);
		
		//Recording of EEI values
		if (!energy) {
			//Record the EEI4 values (=2/(Ta+Te))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI4(organism.getGenes().get(i).getEei());
			}
		}
		else {
			//Record EEI5 values (=2/(Ta+Te))
			for (int i = 0; i < organism.getGenes().size(); i++) {
				organism.getGenes().get(i).putEEI5(organism.getGenes().get(i).getEei());
			}
		}
		
		//Mean and displacement of ribosomal gene positions
		return djmass;
	}
	
	//TODO Calculation of EEI = Ta + Te
	private void eei(Organism organism) {
		for (int i = 0; i < organism.getGenes().size(); i++) {
			//Codon composition + hairpins
			double eei = organism.getGenes().get(i).getTai()+organism.getGenes().get(i).getTei();
			organism.getGenes().get(i).setEei(eei);
		}
	}
}
