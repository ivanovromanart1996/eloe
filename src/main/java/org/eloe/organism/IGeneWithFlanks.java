package org.eloe.organism;

public interface IGeneWithFlanks extends IGene {
	
	char[] getNaWithFlanks();
	
	int getLeftFlankLength();

	int getRightFlankLength();
	
}
