package org.eloe.engine.te.lci.lci1;

import org.eloe.other.Auxiliary_Class;
import org.eloe.parameters.PLCI1;

public class Lci1Calculator extends Lci1CalculatorAbstract {

	public Lci1Calculator(PLCI1 parameters) {
		super(parameters);
	}

	@Override
	protected char[] reverseComplement(char[] char1) {
		return Auxiliary_Class.reversComplement(char1);
	}
}
