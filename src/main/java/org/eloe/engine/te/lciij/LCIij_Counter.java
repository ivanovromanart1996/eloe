package org.eloe.engine.te.lciij;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import org.eloe.engine.ta.FreqCounter;
import org.eloe.organism.Organism;
import org.eloe.other.Sorter;
import org.eloe.parameters.PLCI2;

// Програма для подсчета индекса локальной комплементарности

public class LCIij_Counter {
	private static final String ENERGY = "lci2double.txt";
	
	private PLCI2 parameters;
	//private int forSigma2; //Множитель для сигмы
	private double loopenergy[] = null;
	private Map<String, Double> map = new HashMap<String, Double>();
	private String resultFileLCIij; //Файл для записи результатов рассчетов LCIij
	
	//Конструктор
	public LCIij_Counter(PLCI2 parameters, String resultFileLCIij) {
		this.parameters = parameters;
		//this.forSigma2 = parameters.forSigma2;
		this.resultFileLCIij = resultFileLCIij;
		String loopsenergy = "";
		BufferedReader reader = null;
		String line = new String();
		try {
			//Чтение файла с энергиями динуклеотидов и петель
			reader = new BufferedReader(new InputStreamReader(FreqCounter.class.getClassLoader().getResourceAsStream(ENERGY)));
			while (!(line = reader.readLine()).startsWith("//")){
				String arr[] = line.split(" - ");
				this.map.put(arr[0], Double.parseDouble(arr[1]));
			}
			while (!((line = reader.readLine()).startsWith("////"))) {
				loopsenergy += line;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		String loopsenergy1[] = loopsenergy.split(";");
		//Число длин петель
		int massivesize = loopsenergy1.length;
		//Массив для записи энергии петель
		loopenergy = new double[massivesize];
		for (int i = 0; i < massivesize; i++) {
			loopenergy[i] = Double.parseDouble(loopsenergy1[i].trim());
		}
	}
	
	//TODO Расчет LCIij
	public void countLCIij(Organism organism) {
		PrintStream outLCIij = null;
		try {
			outLCIij = new PrintStream(resultFileLCIij);
			if (organism.isEukaryota()) {
				outLCIij.println("#\tLCIij\t(eucaryota)");
			}
			else {
				outLCIij.println("#\tfirst&last\tfirst&!last\tmiddle\tlast&!first\tLCIij\t(procaryota)");
			}
			
			for (int i = 0; i < organism.getGenes().size(); i++) {
				int d = parameters.getLmax2()-parameters.getLmin2();
				int lenght = organism.getGenes().get(i).getLenghtWithFlanks(); //Длина исследуемой части гена 
				char gene[] = organism.getGenes().get(i).getNaWithFlanks(); //Достаем ген с флангами
				int edge = lenght-2*parameters.getSmax2()-parameters.getLmax2()+1; //Граница левого края скользящей рамки
				if (edge < 0) continue;
				double lciij[] = new double[lenght]; // LCI индивидуального нуклеотида
				for (int j = 0; j < edge; j++){
					double energy = 0; //Энергия кодонного взаимодействия
					for (int s = parameters.getSmin2(); s <= parameters.getSmax2(); s++){
						char con1[] = new char[s]; //Строка для повтора
						System.arraycopy(gene, j, con1, 0, s); //Копирование последовательности
						double energy1 = 0; //Энергия кодонного взаимодействия
						for (int m = 0; m < s-1; m++) {
							String dinukleotid = con1[m]+""+con1[m+1];
							if (map.containsKey(dinukleotid)) {
								energy1 += map.get(dinukleotid); //Энергия взаимодействия повторов
							}
						}
						StringBuilder con2 = new StringBuilder(""); //Строка для комплементарной последовательности
						for (int k = s-1; k > -1 ; k--){ //Построение комплементарной последовательности
							if (con1[k] == 'c') {
								con2.append("g");
							}
							else {
								if (con1[k] == 'a') {
									con2.append("t");
								}
								else {
									if (con1[k] == 'g') {
										con2.append("c");
									}
									else con2.append("a");
								}
							}
						}
						String con22 = con2.toString();
						int lenght1 = s+d; //Длина участка поиска
						char con3[] = new char[lenght1]; //Участок поиска
						System.arraycopy(gene, j+s+parameters.getLmin2(), con3, 0, lenght1); //Достаем участок поиска
						String gene1 = String.valueOf(con3); //Участок поиска
						int poz = 0;
						double loop1 = 0; //Энергия петель
						int k = gene1.indexOf(con22, poz);
						while (k != -1) {
							int loopsize = parameters.getLmin2()+k-1;
							//Энергия петель
							loop1 = loopenergy[(loopsize < loopenergy.length) ? loopsize : (loopenergy.length-1)];
							double e = energy1 - loop1;
							//Если энергия больше порога существования петли, тогда учитываем эту петлю
							if (e >= parameters.getMinloopEnergy()) {
								energy += e;
							}
							poz = k+1;
							k = gene1.indexOf(con22, poz);
						}
						for (int jj = j; jj < j+s; jj++) {
							lciij[jj] += 10*energy;
						}
					}
				}
				//Запись индекса локальной комплементарности индивидуального нуклеотида
				StringBuilder s = new StringBuilder("");
				s.append(organism.getGenes().get(i).getNumber());
				if (!organism.isEukaryota()) {
					if (organism.getGenes().get(i).isFirstCistron() && organism.getGenes().get(i).isLastCistron()) {
						s.append("\t1\t0\t0\t0");
					}
					else {
						if (organism.getGenes().get(i).isMiddleCistron()) {
							s.append("\t0\t0\t1\t0");
						}
						else {
							if (organism.getGenes().get(i).isFirstCistron() && !organism.getGenes().get(i).isLastCistron()) {
								s.append("\t0\t1\t0\t0");
							}
							else {
								s.append("\t0\t0\t0\t1");
							}
						}
					}
				}
				for (int j = 0; j < lciij.length; j++) {
					s.append("\t");
					s.append(lciij[j]);
				}
				outLCIij.println(s.toString());
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (outLCIij != null) {
				outLCIij.close();
			}
		}
	}	
}
