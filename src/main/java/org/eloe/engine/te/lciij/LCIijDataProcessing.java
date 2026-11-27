package org.eloe.engine.te.lciij;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.log4j.Logger;
import org.apache.poi.ss.usermodel.Chart;
import org.apache.poi.ss.usermodel.ClientAnchor;
import org.apache.poi.ss.usermodel.Drawing;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.charts.AxisCrosses;
import org.apache.poi.ss.usermodel.charts.AxisPosition;
import org.apache.poi.ss.usermodel.charts.ChartAxis;
import org.apache.poi.ss.usermodel.charts.ChartDataSource;
import org.apache.poi.ss.usermodel.charts.DataSources;
import org.apache.poi.ss.usermodel.charts.LineChartData;
import org.apache.poi.ss.usermodel.charts.ValueAxis;
import org.apache.poi.ss.util.CellRangeAddress;
import org.apache.poi.ss.util.WorkbookUtil;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import org.eloe.other.Auxiliary_Class;

//TODO Reading of LCIij values and calculation of correlation
public class LCIijDataProcessing {
	
	private static final Logger LOG = Logger.getLogger(LCIijDataProcessing.class);
	
	//Additional data for prokaryotes
	private static final int FIRSTPRO = 5;
	//Additional data for eukaryotes
	private static final int FIRSTEU = 1;
	//Minimum gene length
	private static final int BETWEEN = 3;
	//Size of flanks
	private static final int FLANGSIZE = 600;
	
	public static void readLCIij(File fileInput1, File fileInput2, File fileOutput) {
		boolean procaryota = false;
		int counter = 0;
		int first = 0;
		long time = System.currentTimeMillis();
		int razmer = 2*FLANGSIZE+BETWEEN;
		List<double[]> matrixL = new ArrayList<double[]>();
		List<double[]> matrixR = new ArrayList<double[]>();
		List<double[]> matrixEEI = new ArrayList<double[]>();
		
		BufferedReader readerLCIij = null;
		BufferedReader readerEEI = null;
		try {
			//First file.
			readerLCIij = new BufferedReader(new FileReader(fileInput1));
			String line = readerLCIij.readLine();
			if (line.contains("eucaryota")) {
				first = FIRSTEU;
				procaryota = false;
			}
			else {
				first = FIRSTPRO;
				procaryota = true;
			}
			while ((line = readerLCIij.readLine()) != null) {
				String massString[] = line.split("\t");
				int pozition = massString.length-(first+razmer);
				if (pozition >= 0) {
					double massDoubleL[] = new double[first+razmer];
					double massDoubleR[] = new double[first+razmer];
					for (int i = 0; i < first; i++) {
						massDoubleL[i] = Double.parseDouble(massString[i]);
						massDoubleR[i] = Double.parseDouble(massString[i]);
					}
					for (int i = first; i < first+razmer; i++) {
						massDoubleL[i] = Double.parseDouble(massString[i]);
						massDoubleR[i] = Double.parseDouble(massString[pozition+i]);
					}
					matrixL.add(massDoubleL);
					matrixR.add(massDoubleR);
				}
				else {
					counter++;
				}
			}
			
			//Second file.
			readerEEI = new BufferedReader(new FileReader(fileInput2));
			line = readerEEI.readLine();
			while ((line = readerEEI.readLine()) != null) {
				String massString[] = line.split("\t");
				double massDoubleEEI[] = new double[first+1];
				for (int i = 0; i < massString.length; i++) {
					massDoubleEEI[i] = Double.parseDouble(massString[i]);
				}
				matrixEEI.add(massDoubleEEI);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (readerLCIij != null) {
				try {
					readerLCIij.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			if (readerEEI != null) {
				try {
					readerEEI.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		time = Auxiliary_Class.time("reading LCIij data", time);
		
		//TODO Calculation of correlation
		//Average EEI
		double sredv = 0;
		if (matrixEEI.size() == matrixL.size()) {
			for (int i = 0; i < matrixEEI.size(); i++) {
				sredv += matrixEEI.get(i)[first];
			}
			sredv /= matrixEEI.size();
		}
		else {
			counter++;
			LOG.error("matrixEEI.size() != matrixL.size()\n"
					+ "matrixEEI.size(): "+matrixEEI.size()+" matrixL.size(): "+matrixL.size());
		}
		//Average LCIij values of first and last nucleotides
		double sredmL[] = new double[razmer];
		double sredmLFandL[] = new double[razmer];
		double sredmLFandNoL[] = new double[razmer];
		double sredmLMid[] = new double[razmer];
		double sredmLLandNoF[] = new double[razmer];
		double sredmR[] = new double[razmer];
		double sredmRFandL[] = new double[razmer];
		double sredmRFandNoL[] = new double[razmer];
		double sredmRMid[] = new double[razmer];
		double sredmRLandNoF[] = new double[razmer];
		int _FandL = 0;
		int _FandNoL = 0;
		int _Mid = 0;
		int _LandNoF = 0;
		if (matrixL.size() == matrixR.size()) {
			for (int j = 0; j < razmer; j++) {
				for (int i = 0; i < matrixL.size(); i++) {
					sredmL[j] += matrixL.get(i)[first+j];
					sredmR[j] += matrixR.get(i)[first+j];
					if (procaryota) {
						if (matrixL.get(i)[1] == 1.) { //If the gene is first&last
							if (j == 0) _FandL++;
							sredmLFandL[j] += matrixL.get(i)[first+j];
							sredmRFandL[j] += matrixR.get(i)[first+j];
						}
						else {
							if (matrixL.get(i)[3] == 1.) { //If the gene is middle
								if (j == 0) _Mid++;
								sredmLMid[j] += matrixL.get(i)[first+j];
								sredmRMid[j] += matrixR.get(i)[first+j];
							}
							else {
								if (matrixL.get(i)[2] == 1.) { //If the gene is first&!last
									if (j == 0) _FandNoL++;
									sredmLFandNoL[j] += matrixL.get(i)[first+j];
									sredmRFandNoL[j] += matrixR.get(i)[first+j];
								}
								else { //If the gene is last&!first
									if (j == 0) _LandNoF++;
									sredmLLandNoF[j] += matrixL.get(i)[first+j];
									sredmRLandNoF[j] += matrixR.get(i)[first+j];
								}
							}
						}
					}
				}
				sredmL[j] /= matrixL.size();
				sredmR[j] /= matrixR.size();
				if (procaryota) {
					sredmLFandL[j] /= matrixL.size();
					sredmRFandL[j] /= matrixR.size();
					sredmLMid[j] /= matrixL.size();
					sredmRMid[j] /= matrixR.size();
					sredmLFandNoL[j] /= matrixL.size();
					sredmRFandNoL[j] /= matrixR.size();
					sredmLLandNoF[j] /= matrixL.size();
					sredmRLandNoF[j] /= matrixR.size();
				}
			}
		}
		else {
			counter++;
			System.out.println("matrixL.size() != matrixR.size(), || sredm1.length != sredm2.length");
			System.out.println("matrixL.size(): "+matrixL.size()+" matrixR.size(): "+matrixR.size());
			LOG.error("matrixL.size() != matrixR.size(), || sredm1.length != sredm2.length\nmatrixL.size(): "+matrixL.size()+" matrixR.size(): "+matrixR.size());
		}
		//Delta vector
		double deltv = 0;
		for (int i = 0; i < matrixEEI.size(); i++) {
			deltv += Math.pow((matrixEEI.get(i)[first]-sredv), 2);
		}
		deltv = Math.sqrt(deltv/matrixEEI.size());
		//Delta of the first matrix
		double deltm1[] = new double[razmer];
		Arrays.fill(deltm1, 0);
		for (int j = 0; j < razmer; j++) {
			for (int i = 0; i < matrixL.size(); i++) {
				deltm1[j] += Math.pow((matrixL.get(i)[first+j]-sredmL[j]), 2);
			}
			deltm1[j] = Math.sqrt(deltm1[j]/matrixL.size());
		}
		//Delta of the second matrix
		double deltm2[] = new double[razmer];
		Arrays.fill(deltm2, 0);
		for (int j = 0; j < razmer; j++) {
			for (int i = 0; i < matrixR.size(); i++) {
				deltm2[j] += Math.pow((matrixR.get(i)[first+j]-sredmR[j]), 2);
			}
			deltm2[j] = Math.sqrt(deltm2[j]/matrixR.size());
		}
		//r1
		double r1[] = new double[razmer];
		Arrays.fill(r1, 0);
		for (int j = 0; j < razmer; j++) {
			for (int i = 0; i < matrixL.size(); i++) {
				r1[j] += (matrixL.get(i)[first+j]-sredmL[j])*(matrixEEI.get(i)[first]-sredv);
			}
			r1[j] = r1[j]/(deltm1[j]*deltv*matrixL.size());
		}
		//r2
		double r2[] = new double[razmer];
		Arrays.fill(r2, 0);
		for (int j = 0; j < razmer; j++) {
			for (int i = 0; i < matrixR.size(); i++) {
				r2[j] += (matrixR.get(i)[first+j]-sredmR[j])*(matrixEEI.get(i)[first]-sredv);
			}
			r2[j] = r2[j]/(deltm2[j]*deltv*matrixR.size());
		}
		//t1
		double t1[] = new double[razmer];
		for (int j = 0; j < t1.length; j++) {
			t1[j] = Math.abs(r1[j])/Math.sqrt((1-r1[j]*r1[j])/(matrixL.size()-2));
		}
		//t2
		double t2[] = new double[razmer];
		for (int j = 0; j < t2.length; j++) {
			t2[j] = Math.abs(r2[j])/Math.sqrt((1-r2[j]*r2[j])/(matrixR.size()-2));
		}
		time = Auxiliary_Class.time("counting LCIij profiles", time);
		
		//=============================================================================================================
		//=============================================================================================================
		//=============================================================================================================
		//TODO Output of calculations into file
		
		int shift = 4;
		
		FileOutputStream fileOut = null;
		Workbook wb = new XSSFWorkbook();
		String safeName = WorkbookUtil.createSafeSheetName("new sheet");
		Sheet sheet = wb.createSheet(safeName);
		
		Row row0 = sheet.createRow(0);
		row0.createCell(0).setCellValue("Number of genes:");
		row0.createCell(1).setCellValue(matrixEEI.size());
		row0.createCell(2).setCellValue("Dostovernost:");
		row0.createCell(3).setCellValue(0.05);
		row0.createCell(4).setCellValue("STUDRASPOBR:");
		row0.createCell(5).setCellFormula("TINV(D1,B1-2)");
		
		Row row1 = sheet.createRow(1);
		row1.createCell(0).setCellValue("Left flank");
		
		Row row2 = sheet.createRow(2);
		row2.createCell(0).setCellValue("All genes:");
		row2.createCell(1).setCellValue(matrixL.size());
		for (int i = 0; i < sredmL.length-2; i++) {
			row2.createCell(i+2).setCellValue(-600+i);
		}
		
		Row row3 = sheet.createRow(3);
		row3.createCell(0).setCellValue("all");
		row3.createCell(1).setCellValue(matrixL.size());
		for (int i = 0; i < sredmL.length; i++) {
			row3.createCell(i+2).setCellValue(sredmL[i]);
		}
		
		if (procaryota) {
			shift = 0;
			Row row4 = sheet.createRow(4);
			row4.createCell(0).setCellValue("first&last");
			row4.createCell(1).setCellValue(_FandL);
			for (int i = 0; i < sredmLFandL.length; i++) {
				row4.createCell(i+2).setCellValue(sredmLFandL[i]);
			}
			
			Row row5 = sheet.createRow(5);
			row5.createCell(0).setCellValue("first&!last");
			row5.createCell(1).setCellValue(_FandNoL);
			for (int i = 0; i < sredmLFandNoL.length; i++) {
				row5.createCell(i+2).setCellValue(sredmLFandNoL[i]);
			}
			
			Row row6 = sheet.createRow(6);
			row6.createCell(0).setCellValue("middle");
			row6.createCell(1).setCellValue(_Mid);
			for (int i = 0; i < sredmLMid.length; i++) {
				row6.createCell(i+2).setCellValue(sredmLMid[i]);
			}
			
			Row row7 = sheet.createRow(7);
			row7.createCell(0).setCellValue("last&!first");
			row7.createCell(1).setCellValue(_LandNoF);
			for (int i = 0; i < sredmLLandNoF.length; i++) {
				row7.createCell(i+2).setCellValue(sredmLLandNoF[i]);
			}
		}		
		
		Row row8 = sheet.createRow(8-shift);
		row8.createCell(0).setCellValue("r1");
		for (int i = 0; i < r1.length; i++) {
			row8.createCell(i+2).setCellValue(r1[i]);
		}
		
		Row row9 = sheet.createRow(9-shift);
		row9.createCell(0).setCellValue("t1");
		for (int i = 0; i < t1.length; i++) {
			row9.createCell(i+2).setCellValue(t1[i]);
		}
		
		Row row10 = sheet.createRow(10-shift);
		row10.createCell(0).setCellValue("dostovernost");
		if (procaryota) {
			row10.createCell(2).setCellValue("=ЕСЛИ(C10>$F$1;1;0)");
		}
		else {
			row10.createCell(2).setCellValue("=ЕСЛИ(C6>$F$1;1;0)");
		}
		
		Row row11 = sheet.createRow(11-shift);
		row11.createCell(0).setCellValue("Right flank");
		
		Row row12 = sheet.createRow(12-shift);
		row12.createCell(0).setCellValue("All genes:");
		row12.createCell(1).setCellValue(matrixR.size());
		
		Row row13 = sheet.createRow(13-shift);
		row13.createCell(0).setCellValue("all");
		row13.createCell(1).setCellValue(matrixR.size());
		for (int i = 0; i < sredmR.length; i++) {
			row13.createCell(i+2).setCellValue(sredmR[i]);
		}
		
		shift = 8;
		if (procaryota) {
			shift = 0;
			Row row14 = sheet.createRow(14);
			row14.createCell(0).setCellValue("first&last");
			row14.createCell(1).setCellValue(_FandL);
			for (int i = 0; i < sredmRFandL.length; i++) {
				row14.createCell(i+2).setCellValue(sredmRFandL[i]);
			}
			
			Row row15 = sheet.createRow(15);
			row15.createCell(0).setCellValue("first&!last");
			row15.createCell(1).setCellValue(_FandNoL);
			for (int i = 0; i < sredmRFandNoL.length; i++) {
				row15.createCell(i+2).setCellValue(sredmRFandNoL[i]);
			}
			
			Row row16 = sheet.createRow(16);
			row16.createCell(0).setCellValue("middle");
			row16.createCell(1).setCellValue(_Mid);
			for (int i = 0; i < sredmRMid.length; i++) {
				row16.createCell(i+2).setCellValue(sredmRMid[i]);
			}
			
			Row row17 = sheet.createRow(17);
			row17.createCell(0).setCellValue("last&!first");
			row17.createCell(1).setCellValue(_LandNoF);
			for (int i = 0; i < sredmRLandNoF.length; i++) {
				row17.createCell(i+2).setCellValue(sredmRLandNoF[i]);
			}
		}		
		
		Row row18 = sheet.createRow(18-shift);
		row18.createCell(0).setCellValue("r2");
		for (int i = 0; i < r2.length; i++) {
			row18.createCell(i+2).setCellValue(r2[i]);
		}
		
		Row row19 = sheet.createRow(19-shift);
		row19.createCell(0).setCellValue("t2");
		for (int i = 0; i < t2.length; i++) {
			row19.createCell(i+2).setCellValue(t2[i]);
		}
		
		Row row20 = sheet.createRow(20-shift);
		row20.createCell(0).setCellValue("dostovernost");
		if (procaryota) {
			row20.createCell(2).setCellValue("=ЕСЛИ(C20>$F$1;1;0)");
		}
		else {
			row20.createCell(2).setCellValue("=ЕСЛИ(C12>$F$1;1;0)");
		}
		
		Drawing drawing = sheet.createDrawingPatriarch();
		
		ClientAnchor anchor1 = drawing.createAnchor(0, 0, 0, 0, 2, 22-shift, 10, 37-shift);
		Chart chart1 = drawing.createChart(anchor1);
		LineChartData data1 = chart1.getChartDataFactory().createLineChartData();
		ChartAxis bottomAxis1 = chart1.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
		ValueAxis leftAxis1 = chart1.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
		leftAxis1.setCrosses(AxisCrosses.AUTO_ZERO);
		ChartDataSource<Number> xs = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(2, 2, 2, sredmL.length-1));
		ChartDataSource<Number> ys1 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(3, 3, 2, sredmL.length-1));
		data1.addSerie(xs, ys1);
		chart1.plot(data1, bottomAxis1, leftAxis1);
		
		Row row22 = sheet.createRow(22-shift);
		row22.createCell(10).setCellValue("all");
		
		ClientAnchor anchor2 = drawing.createAnchor(0, 0, 0, 0, 11, 22-shift, 19, 37-shift);
		Chart chart2 = drawing.createChart(anchor2);
		LineChartData data2 = chart2.getChartDataFactory().createLineChartData();
		ChartAxis bottomAxis2 = chart2.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
		ValueAxis leftAxis2 = chart2.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
		leftAxis2.setCrosses(AxisCrosses.AUTO_ZERO);
		
		shift = 4;
		if (procaryota) shift = 0;
		
		ChartDataSource<Number> ys2 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(13-shift, 13-shift, 2, sredmR.length-1));
		data2.addSerie(xs, ys2);
		chart2.plot(data2, bottomAxis2, leftAxis2);
		
		if (procaryota) {
			ClientAnchor anchor3 = drawing.createAnchor(0, 0, 0, 0, 2, 38, 10, 53);
			Chart chart3 = drawing.createChart(anchor3);
			LineChartData data3 = chart3.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis3 = chart3.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis3 = chart3.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis3.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys3 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(4, 4, 2, sredmR.length-1));
			data3.addSerie(xs, ys3);
			chart3.plot(data3, bottomAxis3, leftAxis3);
			
			Row row38 = sheet.createRow(38);
			row38.createCell(10).setCellValue("first&last");
			
			ClientAnchor anchor4 = drawing.createAnchor(0, 0, 0, 0, 11, 38, 19, 53);
			Chart chart4 = drawing.createChart(anchor4);
			LineChartData data4 = chart4.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis4 = chart4.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis4 = chart4.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis4.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys4 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(14, 14, 2, sredmR.length-1));
			data4.addSerie(xs, ys4);
			chart4.plot(data4, bottomAxis4, leftAxis4);
			
			ClientAnchor anchor5 = drawing.createAnchor(0, 0, 0, 0, 2, 54, 10, 69);
			Chart chart5 = drawing.createChart(anchor5);
			LineChartData data5 = chart5.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis5 = chart5.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis5 = chart5.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis5.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys5 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(5, 5, 2, sredmR.length-1));
			data5.addSerie(xs, ys5);
			chart5.plot(data5, bottomAxis5, leftAxis5);
			
			Row row54 = sheet.createRow(54);
			row54.createCell(10).setCellValue("first&!last");
			
			ClientAnchor anchor6 = drawing.createAnchor(0, 0, 0, 0, 11, 54, 19, 69);
			Chart chart6 = drawing.createChart(anchor6);
			LineChartData data6 = chart6.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis6 = chart6.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis6 = chart6.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis6.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys6 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(15, 15, 2, sredmR.length-1));
			data6.addSerie(xs, ys6);
			chart6.plot(data6, bottomAxis6, leftAxis6);
			
			ClientAnchor anchor7 = drawing.createAnchor(0, 0, 0, 0, 2, 70, 10, 85);
			Chart chart7 = drawing.createChart(anchor7);
			LineChartData data7 = chart7.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis7 = chart7.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis7 = chart7.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis7.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys7 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(6, 6, 2, sredmR.length-1));
			data7.addSerie(xs, ys7);
			chart7.plot(data7, bottomAxis7, leftAxis7);
			
			Row row70 = sheet.createRow(70);
			row70.createCell(10).setCellValue("middle");
			
			ClientAnchor anchor8 = drawing.createAnchor(0, 0, 0, 0, 11, 70, 19, 85);
			Chart chart8 = drawing.createChart(anchor8);
			LineChartData data8 = chart8.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis8 = chart8.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis8 = chart8.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis8.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys8 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(16, 16, 2, sredmR.length-1));
			data8.addSerie(xs, ys8);
			chart8.plot(data8, bottomAxis8, leftAxis8);
			
			ClientAnchor anchor9 = drawing.createAnchor(0, 0, 0, 0, 2, 86, 10, 101);
			Chart chart9 = drawing.createChart(anchor9);
			LineChartData data9 = chart9.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis9 = chart9.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis9 = chart9.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis9.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys9 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(7, 7, 2, sredmR.length-1));
			data9.addSerie(xs, ys9);
			chart9.plot(data9, bottomAxis9, leftAxis9);
			
			Row row86 = sheet.createRow(86);
			row86.createCell(10).setCellValue("last&!first");
			
			ClientAnchor anchor10 = drawing.createAnchor(0, 0, 0, 0, 11, 86, 19, 101);
			Chart chart10 = drawing.createChart(anchor10);
			LineChartData data10 = chart10.getChartDataFactory().createLineChartData();
			ChartAxis bottomAxis10 = chart10.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
			ValueAxis leftAxis10 = chart10.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
			leftAxis10.setCrosses(AxisCrosses.AUTO_ZERO);
			ChartDataSource<Number> ys10 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(17, 17, 2, sredmR.length-1));
			data10.addSerie(xs, ys10);
			chart10.plot(data10, bottomAxis10, leftAxis10);
		}
		
		int shift2 = 72;
		if (procaryota) shift2 = 0;
		
		ClientAnchor anchor11 = drawing.createAnchor(0, 0, 0, 0, 2, 102-shift2, 10, 117-shift2);
		Chart chart11 = drawing.createChart(anchor11);
		LineChartData data11 = chart11.getChartDataFactory().createLineChartData();
		ChartAxis bottomAxis11 = chart11.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
		ValueAxis leftAxis11 = chart11.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
		leftAxis11.setCrosses(AxisCrosses.AUTO_ZERO);
		ChartDataSource<Number> ys11 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(8-shift, 8-shift, 2, sredmR.length-1));
		data11.addSerie(xs, ys11);
		chart11.plot(data11, bottomAxis11, leftAxis11);
		
		Row row102 = sheet.createRow(102-shift2);
		row102.createCell(10).setCellValue("correlation");
		
		shift = 8;
		if (procaryota) shift = 0;
		
		ClientAnchor anchor10 = drawing.createAnchor(0, 0, 0, 0, 11, 102-shift2, 19, 117-shift2);
		Chart chart10 = drawing.createChart(anchor10);
		LineChartData data10 = chart10.getChartDataFactory().createLineChartData();
		ChartAxis bottomAxis10 = chart10.getChartAxisFactory().createCategoryAxis(AxisPosition.BOTTOM);
		ValueAxis leftAxis10 = chart10.getChartAxisFactory().createValueAxis(AxisPosition.LEFT);
		leftAxis10.setCrosses(AxisCrosses.AUTO_ZERO);
		ChartDataSource<Number> ys10 = DataSources.fromNumericCellRange(sheet, new CellRangeAddress(18-shift, 18-shift, 2, sredmR.length-1));
		data10.addSerie(xs, ys10);
		chart10.plot(data10, bottomAxis10, leftAxis10);
		
		try {
			fileOut = new FileOutputStream(fileOutput.toString().replace(".txt", ".xlsx"));
			wb.write(fileOut);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (fileOut != null) {
				try {
					fileOut.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		//=============================================================================================================
		//=============================================================================================================
		//=============================================================================================================

		Auxiliary_Class.time("out LCIij profiles", time);
		String say = counter == 0 ? "OK" : "NO";
		LOG.info("LCIij profile "+say);
	}
}
