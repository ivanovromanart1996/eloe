package org.eloe.parameters;

public interface PProcessing {

	boolean isCheckStart();

	String[] getStartCodons();

	boolean isDeleteGwBC();

	String[] getMidlCodons();

	boolean isCheckStop();

	String[] getStopCodons();
		
	boolean isUsePseudoGenes();

}