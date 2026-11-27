package org.eloe.extractor;

import java.io.File;
import java.util.regex.Pattern;

import org.eloe.parameters.PExtraction;

public abstract class ExtractorFactory {
	private static final Pattern GENBANK = Pattern.compile(".+\\.gb[kf]*$");//old: ".+\\.gb(k)?$"
	private static final Pattern EMBL = Pattern.compile(".+\\.emb$");
	
	public static Extractor getExtractor(File orgFile, PExtraction parameters, File resultFolder) {
		if (GENBANK.matcher(orgFile.toString()).matches()) {
			return new GBKFileExtractor(orgFile, parameters, resultFolder);
		}
		else if (EMBL.matcher(orgFile.toString()).matches()) {
			return new EMBFileExtractor(orgFile, parameters, resultFolder);
		}
		else return null;
	}
}
