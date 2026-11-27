package org.eloe.math;

//Integral calculation by Gauss-Christoffel method with roots of Legendre polynomials in power 7
public class GaussChristoffel7 extends Integral {
	//Корни полинома Лежандра степени 7
	public static final double X[] = {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0000000000000000, 
		0.4058451513773972, 0.7415311855993945, 0.9491079123427585};
	//Весовые коэффициенты 7s
	public static final double A[] = {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 
		0.3818300505051189, 0.2797053914892766, 0.1294849661688697};
		
	@Override
	protected final double integral(double a, double b, Function f) {
		double apb = (a+b)/2;
		double bma = (b-a)/2;
		double sum = 0;
		for (int i = 0; i < A.length; i++) {
			sum += A[i]*f.f(apb+bma*X[i]);
		}
		return bma*sum;
	}
}
