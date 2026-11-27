package org.eloe.engine.te.lci;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import org.eloe.engine.ta.FreqCounter;
import org.eloe.organism.Organism;
import org.eloe.parameters.PLCI2;

// Програма для подсчета индекса локальной комплементарности

public class LCI12_Counter {
	
	private static final Logger LOG = Logger.getLogger(LCI12_Counter.class);

	private static final String ENERGY = "lci2double.txt";
	private PLCI2 parameters;
	//private int forSigma1; //Множитель для сигмы
	//private int forSigma2; //Множитель для сигмы
	private double m1 = 0; //Среднее значение LCI
	private double sigma1 = 0; //Среднее отклонение LCI
	private double m2 = 0; //Среднее значение LCI
	private double sigma2 = 0; //Среднее отклонение LCI
	private double sred1 = 0;
	private double disp1 = 0;
	private double sred2 = 0;
	private double disp2 = 0;
	private double loopenergy[] = null;
	private Map<String, Double> map = new HashMap<String, Double>(); //Карта для записи энергии динуклеотидов
	
	//Конструктор
	public LCI12_Counter(PLCI2 parameters) throws IOException, NumberFormatException {
		this.parameters = parameters;
		//Вспомогательная переменная
		//this.forSigma1 = parameters.forSigma1;
		//this.forSigma2 = parameters.forSigma2;
		String loopsenergy = "";
		
		BufferedReader reader = null;
		String line = new String();
		try {
			//Чтение файла с энергиями динуклеотидов и петель
			reader = new BufferedReader(new InputStreamReader(FreqCounter.class.getClassLoader().getResourceAsStream(ENERGY)));
			while (!(line = reader.readLine()).startsWith("//")){
				String arr[] = line.split(" - ");
				double a = Double.parseDouble(arr[1]);
				map.put(arr[0], a);
			}
			while (!((line = reader.readLine()).startsWith("////"))) {
				loopsenergy += line;
			}
		} catch (NumberFormatException e) {
			e.printStackTrace();
			LOG.error("Sorry, something wrong!", e);
			throw e;
		} catch (IOException e) {
			e.printStackTrace();
			LOG.error("Sorry, something wrong!", e);
			throw e;
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
		try {
			for (int i = 0; i < massivesize; i++) {
				loopenergy[i] = Double.parseDouble(loopsenergy1[i].trim());
			}
		} catch (NumberFormatException e) {
			e.printStackTrace();
			LOG.error("Sorry, something wrong!", e);
			throw e;
		}
	}
	
	//TODO Рассчет LCI12
	public void countLCI(Organism organism) {
		//LCI1
		double lci1 = 0;
		//LCI2
		double lci2 = 0;
		
		for (int i = 0; i < organism.getGenes().size(); i++) {
			//Заступ справа за край гена
			int h = parameters.getSmax2()+parameters.getLmax2();
			if (h > organism.getGenes().get(i).getRightFlankLength()) {
				//h = rightFl;
				continue;
			}
			
			//Длина исследуемой части гена (3 - стоп-кодон)
			int lenght = organism.getGenes().get(i).getLength()-3+h;
			
			//Строка для записи исследуемого гена
			char gene[] = new char[lenght];
			
			//Длина левого фланга
			int leftFl = organism.getGenes().get(i).getLeftFlankLength();
			
			//Достаем ген без флангов + h нуклеотидов на 3' конце
			System.arraycopy(organism.getGenes().get(i).getNaWithFlanks(), leftFl, gene, 0, lenght);
			
			//Граница левого края скользящей рамки
			int edge = lenght-2*parameters.getSmax2()-parameters.getLmax2()+1;
			if (edge < 0) continue;
			
			//Полная энергия взаимодействия
			double energyFull = 0;
			
			//Счетчик повторов
			int counter = 0;
			
			for (int j = 0; j < edge; j++){
				//Энергия шпилек
				double energy = 0;
				
				for (int s = parameters.getSmin2(); s <= parameters.getSmax2(); s++){
					//Строка для повтора
					char con1[] = new char[s];
					//Копирование последовательности
					System.arraycopy(gene, j, con1, 0, s);
					//Энергия выбранного повтора
					double energy_con1 = 0;
					
					for (int m = 0; m < s-1; m++) {
						String dinukleotid = con1[m]+""+con1[m+1];
						if (map.containsKey(dinukleotid)) {
							//Энергия взаимодействия повторов
							energy_con1 += map.get(dinukleotid);
						}
					}
					
					//Строка для комплементарной последовательности
					StringBuilder con2 = new StringBuilder("");
					
					//Построение комплементарной последовательности
					for (int k = s-1; k > -1 ; k--){
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
					
					//Длина участка поиска
					int lenght1 = s+this.parameters.getLmax2()-this.parameters.getLmin2();
					//Массив для участка поиска
					char con3[] = new char[lenght1];
					//Достаем участок поиска
					System.arraycopy(gene, j+s+parameters.getLmin2(), con3, 0, lenght1);
					//Участок поиска
					String gene1 = String.valueOf(con3);
					
					int poz = 0;
					//Энергия петель
					double loop1 = 0;
					
					int k = gene1.indexOf(con22, poz);
					while (k != -1) {
						//Подсчет повторов
						counter++;
						
						int loopsize = parameters.getLmin2()+k-1;
						
						//Энергия петель
						loop1 = loopenergy[(loopsize < loopenergy.length) ? loopsize : (loopenergy.length-1)];
						double e = energy_con1 - loop1;
						//Если энергия больше порога существоавания петли, тогда учитываем эту петлю
						if (e >= parameters.getMinloopEnergy()) {
							energy += e;
						}
						poz = k+1;
						k = gene1.indexOf(con22, poz);
					}
				}
				
				//Полная энергия
				energyFull += energy;
			}
			
			//Число повторов на один нуклеотид (3 - стоп-кодон)
			lci1 = (double)counter/(organism.getGenes().get(i).getLength()-3);
			//Запись LCI1 в карточку гена
			organism.getGenes().get(i).setLci1(lci1);
			this.sred1 += lci1;
			this.disp1 += lci1*lci1;
			
			//Энергия на один нуклеотид (3 - стоп-кодон)
			lci2 = energyFull/(organism.getGenes().get(i).getLength()-3);
			//Запись lci2
			organism.getGenes().get(i).setLci2(lci2);
			this.sred2 += lci2;
			this.disp2 += lci2*lci2;
		}
		
		this.sred1 /= organism.getGenes().size();
		this.disp1 = this.disp1/organism.getGenes().size() - this.sred1*this.sred1;
		this.disp1 = Math.pow(this.disp1, 0.5);
		
		//фортран
		this.m1 = this.sred1 + 2*this.disp1;
		this.sigma1 = 4*this.disp1;
		
		//статьи
		//this.m1 = this.sred1;
		//this.sigma1 = forSigma1*this.disp1;
		
		this.sred2 /= organism.getGenes().size();
		this.disp2 = this.disp2/organism.getGenes().size() - this.sred2*this.sred2;
		this.disp2 = Math.pow(this.disp2, 0.5);
		
		//фортран
		this.m2 = this.sred2 + this.disp2;
		this.sigma2 = 2*this.disp2;
		
		//статьи
		//this.m2 = this.sred2;
		//this.sigma2 = forSigma2*this.disp2;
	}
	
	//TODO Возвращает среднее значение LCI
	public double getM1() {
		return this.m1;
	}
	
	//TODO Возвращает среднее отклонение LCI
	public double getSigma1() {
		return this.sigma1;
	}
	
	//TODO Возвращает среднее значение LCI
	public double getM2() {
		return this.m2;
	}
	
	//TODO Возвращает среднее отклонение LCI
	public double getSigma2() {
		return this.sigma2;
	}
}
