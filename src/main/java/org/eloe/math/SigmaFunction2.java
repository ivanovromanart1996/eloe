package org.eloe.math;

public class SigmaFunction2 implements Function {
	private double n;
	
	public SigmaFunction2(double n) {
		this.n = n;
	}

	@Override
	public double f(double x) {
		return Math.pow(x, n)/Math.exp(x);
	}

}
