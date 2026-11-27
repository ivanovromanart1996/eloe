package org.eloe.other;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JTextArea;

import org.apache.log4j.Logger;

import org.eloe.parameters.POrthologs;

public class Orthologs {
	
	private static final Logger LOG = Logger.getLogger(Orthologs.class);
	
	public static void search(List<String> files, File resultFolder, POrthologs parameters, JTextArea console) {
		if (files.size() > 1) {
			if (files.size() > 10) {
				console.append("\nToo much organisms for orthologs search. It should be <= 10");
				LOG.warn("Too much organisms for orthologs search. It should be <= 10");
				return;
			}
			File result_file = new File(resultFolder, "orthologs.xls");
			PrintStream out = null;
			try {
				out = new PrintStream(result_file);
				//Delete existing records
				out.print("");
			} catch (IOException e) {
				e.printStackTrace();
			} finally {
				if (out != null) {
					out.close();
				}
			}
			
			for (int i = 0; i < files.size()-1; i++) {
				String first = files.get(i);
				//List of names of genes of the first organism
				List<String> first_organism = readGeneNames(first);
				//Name of the first organism
				String first_organism_name = first.substring(first.lastIndexOf(File.separator)+1);
				for (int j = i+1; j < files.size(); j++) {
					String second = files.get(j);
					//List of names of genes of the second organism
					List<String> second_organism = readGeneNames(second);
					//Name of the second organism
					String second_organism_name = second.substring(second.lastIndexOf(File.separator)+1);
					orthologs(first_organism, second_organism, first_organism_name, second_organism_name, result_file, parameters.getNumber_orthologs());
				}
			}
		}
	}
	
	private static void orthologs(List<String> first_organism, List<String> second_organism, 
			String first_organism_name, String second_organism_name, File result_file, int n) {
		PrintStream out = null;
		try {
			out = new PrintStream(new BufferedOutputStream(new FileOutputStream(result_file, true)));
			out.println(first_organism_name+": "+first_organism.size());
			out.println(second_organism_name+": "+second_organism.size());
			//Find orthologs of the first organism in the second.
			searchOrth(first_organism, second_organism, first_organism_name, second_organism_name, n, out);
			//Find orthologs of the second organism in the first.
			searchOrth(second_organism, first_organism, second_organism_name, first_organism_name, n, out);
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (out != null) {
				out.close();
			}
		}
	}

	private static void searchOrth(List<String> first_organism, List<String> second_organism, 
			String first_organism_name, String second_organism_name, int n, PrintStream out) {
		out.println("Top "+n+" from \""+first_organism_name+"\" in \""+second_organism_name+"\"");
		out.println("\t"+first_organism_name+"\t"+second_organism_name);
		
		int sum = 0;
		int sumsq = 0;
		int counter = 0;
		
		for (int i = 0; i < n; i++) {
			//If the gene name is nonempty and the second organism has a gene with that name
			if (!first_organism.get(i).equals("-") && second_organism.contains(first_organism.get(i))) {
				//Record the position of the gene in the second organism
				int poz = second_organism.indexOf(first_organism.get(i))+1;
				counter++;
				sum += poz;
				sumsq += poz*poz;
				//Gene positions in the first and second organism
				out.println((i+1)+"\t"+(i+1)+"\t"+poz+"\t"+first_organism.get(i));
			}
		}
		
		double sred = 0;
		double disp = 0;
		
		if (counter != 0) {
			sred = sum/counter;
			disp = Math.pow((sumsq/counter-sred*sred), 0.5);
			//M
			int m = (int)Math.round(100*(2*sred-second_organism.size()-1)/(second_organism.size()-counter));
			//R
			int r = (int)Math.round(100*(2*disp-counter-1)/(second_organism.size()-counter));
			out.println("M = "+m+"\tR = "+r);
			out.print("\n");
		}
		else {
			LOG.info("There are not ortholog genes!");
		}
	}

	private static List<String> readGeneNames(String org) {		
		//Read the file containing organism data.
		List<String> second_organism = new ArrayList<String>();
		
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(org));
			reader.readLine();
			String line;
			while ((line = reader.readLine()) != null) {
				String[] mass = line.split("\t");
				//Reading a gene name (sixth column in a row)
				second_organism.add(mass[5].toLowerCase());
			}
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
		return second_organism;
	}	
}
