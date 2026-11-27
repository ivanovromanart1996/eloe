package org.eloe.engine.te.lci.lci2;

import java.util.Map;

import org.eloe.other.Auxiliary_Class;
import org.eloe.parameters.PLCI2;

public class Lci2Calculator extends Lci2CalculatorAbstract {

	public Lci2Calculator(PLCI2 parameters, Map<String, Double> map, double[] loopenergy) {
		super(parameters, map, loopenergy);
	}

	@Override
	protected char[] reverseComplement(char[] char1) {
		return Auxiliary_Class.reversComplement(char1);
	}

}
