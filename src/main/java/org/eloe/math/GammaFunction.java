package org.eloe.math;

public abstract class GammaFunction {
	//Coefficients of gamma-function approximation polynomials
	private static final double P[]={-1.71618513886549492533811e+0, 2.47656508055759199108314e+1, -3.79804256470945635097577e+2, 6.29331155312818442661052e+2,
            8.66966202790413211295064e+2, -3.14512729688483675254357e+4, -3.61444134186911729807069e+4, 6.64561438202405440627855e+4};
	private static final double Q[]={-3.08402300119738975254353e+1, 3.15350626979604161529144e+2, -1.01515636749021914166146e+3, -3.10777167157231109440444e+3,
            2.25381184209801510330112e+4, 4.75584627752788110767815e+3, -1.34659959864969306392456e+5, -1.15132259675553483497211e+5};
		
	//A method that calculates the approximate value of the gamma function
	public static double calculate(double z){
	    if ((z >= 1.0) && (z <= 2.0)) return gammaapprox(z); // use an approximation for 1<=z<=2
	    if ((z > 0.0) && (z < 1.0))	return calculate(z+1.0)/z; // Recurrence relation for 0<z<1
	    if (z > 2) return (z-1)*calculate(z-1);	// Recurrence relation for z>2
	    else { // if z<=0
	    	System.out.println("Gamma(z); z <= 0!");
	    	return 0;
	    }
	}
	
	//Calculation of the ratio of gamma function approximation polynomials
	private static double gammaapprox(double z) {
		double x = z-1.0;
	    double a = 0.0;
	    double b = 1.0;
	    for (int i = 0; i < P.length; i++){
	         a =(a+P[i])*x;
	         b = b*x+Q[i];
	    }
	    return (a/b+1.0);
	}
}
