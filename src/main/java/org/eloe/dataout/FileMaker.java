package org.eloe.dataout;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

import org.eloe.organism.Gene;
import org.eloe.organism.Organism;

public class FileMaker {
	
	Organism organism;
	File folder;
	
	public FileMaker(Organism organism, File folder) {
		this.organism = organism;
		this.folder = folder;
	}

	public void makeFile(String fileName, ILineMaker lm) {
		
		PrintStream out = null;
		try {
			out = new PrintStream(new File(folder, fileName));
			for (Gene gene : organism.getGenes()) {
				out.println(lm.makeLine(gene));
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} finally {
			if (out != null)
				out.close();
		}
	}

}
