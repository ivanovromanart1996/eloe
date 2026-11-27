package org.eloe.engine.ta;

import org.eloe.organism.IGene;

public interface ICodonCompositionCalculator {
	
	int[] calculate(IGene gene, int unUsedCodonsL, int unUsedCodonsR);

}
