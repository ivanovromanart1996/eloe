package org.eloe.math;

public abstract class Integral {
	public final double calculate(double a, double b, Function f, double epsilon) {
		long numberOfIntervals = 1;
		double y0;
		double y1 = integral(a, b, f);
		do {
			if (numberOfIntervals >= Long.MAX_VALUE/2) {
				System.out.println("Error! Calculation cann't be finished. Number of intervals is bigger than Long.MAX_VALUE.");
				return 0;
			}
			numberOfIntervals *= 2;
			y0 = y1;
			y1 = 0;
			double h = (b-a)/numberOfIntervals;
			for (int i = 0; i < numberOfIntervals; i++) {
				y1 += integral(a+i*h, a+(i+1)*h, f);
			}
		} while (Math.abs(y1-y0) > epsilon);
		return y1;
	}
	
	protected abstract double integral(double a, double b, Function f);
}
