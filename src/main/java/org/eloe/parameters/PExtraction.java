package org.eloe.parameters;

public interface PExtraction {

	String getContext();

	int getlFlang();

	int getrFlang();

	int getDistBOper();

	int getMinLenGen();
	
	int getSmax1();
	
	int getLmax1();
	
	int getSmax2();
	
	int getLmax2();

	boolean isUsePseudoGenes();

	boolean isMakeFileWholeGenome();

	int getMaxNumbCod();

	boolean isUseRibGenesAsHighExpressed1();

	boolean isUseOtherGenesAsHighExpressed1();

}