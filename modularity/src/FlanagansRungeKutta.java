/*
*   Class RungeKutta
*       requires interfaces DerivFunction and DerivnFunction
*
*   Contains the methods for the Runge-Kutta procedures for solving
*   single or solving sets of ordinary differential equations (ODEs)
*   [draws heavily on the approach adopted in Numerical Recipes
*   (C language version)http://www.nr.com]
*
*   A single ODE is supplied by means of an interface,
*       DerivFunction
*   A set of ODEs is supplied by means of an interface,
*       DerivnFunction
*
*   WRITTEN BY: Dr Michael Thomas Flanagan
*
*   DATE:	    February 2002
*   UPDATES:    22 June 2003,  April 2004,
*               15 September 2006 (to incorporate improvements suggested by Klaus Benary [Klaus.Benary@gede.de])
*
*   DOCUMENTATION:
*   See Michael Thomas Flanagan's Java library on-line web page:
*   RungeKutta.html
*
*   Copyright (c) April 2004, September 2006
*
*   PERMISSION TO COPY:
*   Permission to use, copy and modify this software and its documentation for
*   NON-COMMERCIAL purposes is granted, without fee, provided that an acknowledgement
*   to the author, Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies.
*
*   Dr Michael Thomas Flanagan makes no representations about the suitability
*   or fitness of the software for any or for a particular purpose.
*   Michael Thomas Flanagan shall not be liable for any damages suffered
*   as a result of using, modifying or distributing this software or its derivatives.
*
***************************************************************************************/


// Class for Runge-Kutta solution of ordinary differential equations
public class FlanagansRungeKutta{
        public static double H; //Martin's addition for testing
        public FlanagansRungeKutta(){
        }

        private static double SAFETY=0.9D; // safety scaling factor for Runge Kutta Fehlberg tolerance check

    	// Fourth order Runge-Kutta for n (nequ) ordinary differential equations (ODE)
	    public static double[] fourthOrder(DerivnFunction g, double x0, double[] y0, double xn, double h){
      		int nequ = y0.length;
        	double[] k1 =new double[nequ];
        	double[] k2 =new double[nequ];
        	double[] k3 =new double[nequ];
        	double[] k4 =new double[nequ];
        	double[] y =new double[nequ];
        	double[] yd =new double[nequ];
        	double[] dydx =new double[nequ];
        	double x = 0.0D;

        	// Calculate nsteps
        	double ns = (xn - x0)/h;
        	ns = Math.rint(ns);
        	int nsteps = (int) ns;
        	h = (xn - x0)/ns;

        	// initialise
        	for(int i=0; i<nequ; i++)y[i] = y0[i];

            // iteration over allowed steps
	        for(int j=0; j<nsteps; j++){
        	    	x  = x0 + j*h;
            		dydx = g.derivn(x, y);
            		for(int i=0; i<nequ; i++)k1[i] = h*dydx[i];

	            	for(int i=0; i<nequ; i++)yd[i] = y[i] + k1[i]/2;
        		    dydx = g.derivn(x + h/2, yd);
            		for(int i=0; i<nequ; i++)k2[i] = h*dydx[i];

	            	for(int i=0; i<nequ; i++)yd[i] = y[i] + k2[i]/2;
        	    	dydx = g.derivn(x + h/2, yd);
            		for(int i=0; i<nequ; i++)k3[i] = h*dydx[i];

	            	for(int i=0; i<nequ; i++)yd[i] = y[i] + k3[i];
        	    	dydx = g.derivn(x + h, yd);
            		for(int i=0; i<nequ; i++)k4[i] = h*dydx[i];

	            	for(int i=0; i<nequ; i++)y[i] += k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6;

        	}
        	return y;
    	}

    	// Runge-Kutta-Cash-Karp for n (nequ) ordinary differential equations (ODEs
    	public static double[] cashKarp(DerivnFunction g, double x0, double[] y0, double xn, double h, double abstol, double reltol, int maxiter){
        	int nequ = y0.length;
        	double[] k1 =new double[nequ];
        	double[] k2 =new double[nequ];
        	double[] k3 =new double[nequ];
        	double[] k4 =new double[nequ];
        	double[] k5 =new double[nequ];
        	double[] k6 =new double[nequ];
        	double[] y =new double[nequ];
        	double[] y6 =new double[nequ];
        	double[] y5 =new double[nequ];
        	double[] yd =new double[nequ];
        	double[] dydx =new double[nequ];

        	double x = 0.0D, err = 0.0D, maxerr = 0.0D, delta = 0.0D, tol = 1.0D;
        	int ii = 0;

        	// initialise
        	for(int i=0; i<nequ; i++)y[i] = y0[i];
        	x = x0;

        	while(x<xn){
            		ii++;
            		if(ii>maxiter)throw new ArithmeticException("Maximum number of iterations exceeded");

            		dydx = g.derivn(x, y);
            		for(int i=0; i<nequ; i++)k1[i] = h*dydx[i];

		            for(int i=0; i<nequ; i++)yd[i] = y[i] + k1[i]/5.0;
            		dydx = g.derivn(x + h/5.0, yd);
            		for(int i=0; i<nequ; i++)k2[i] = h*dydx[i];

            		for(int i=0; i<nequ; i++)yd[i] = y[i] + (3.0*k1[i] + 9.0*k2[i])/40.0;
            		dydx = g.derivn(x + 3.0*h/10.0, yd);
            		for(int i=0; i<nequ; i++)k3[i] = h*dydx[i];

            		for(int i=0; i<nequ; i++)yd[i] = y[i] + (3.0*k1[i] - 9.0*k2[i] + 12.0*k3[i])/10.0;
            		dydx = g.derivn(x + 3.0*h/5.0, yd);
            		for(int i=0; i<nequ; i++)k4[i] = h*dydx[i];

	            	for(int i=0; i<nequ; i++)yd[i] = y[i] -11.0*k1[i]/55.0 + 5.0*k2[i]/2.0 - 70.0*k3[i]/27.0 + 35.0*k4[i]/27.0;
            		dydx = g.derivn(x + h, yd);
            		for(int i=0; i<nequ; i++)k5[i] = h*dydx[i];

            		for(int i=0; i<nequ; i++)yd[i] = y[i] + 1631.0*k1[i]/55296.0 + 175.0*k2[i]/512.0 + 575.0*k3[i]/13824.0 + 44275.0*k4[i]/110592.0 + 253.0*k5[i]/4096.0;
            		dydx = g.derivn(x + 7.0*h/8.0, yd);
            		for(int i=0; i<nequ; i++)k6[i] = h*dydx[i];

            		maxerr=0.0D;
            		for(int i=0; i<nequ; i++){
                		y5[i] = y[i] + 2825.0*k1[i]/27648.0 + 18575.0*k3[i]/48384.0 + 13525.0*k4[i]/55296.0 + 277.0*k5[i]/14336.0 + k6[i]/4.0;
                		y6[i] = y[i] + 37*k1[i]/378.0 + 250.0*k3[i]/621.0 + 125.0*k4[i]/594.0  + 512.0*k6[i]/1771.0;
                		err = Math.abs(y6[i] - y5[i]);
                		tol= Math.abs(y5[i])*reltol +  abstol;
                		maxerr = Math.max(maxerr,err/tol);
            		}
            		if(maxerr<=1.0D){
                		x += h;
                		delta = SAFETY*Math.pow(maxerr, -0.2);
                		if(delta>4.0){
                		    h*= 4.0;
                	    }
                    	else if (delta > 1.0){
                    	    h*=delta;
                    	}
                        //Martin's addition for testing - normally the following line should not be commented out
                		//if(x+h > xn)h = xn-x;
                        //Martin's addition for testing (line below)
                        H = h;
                		y = y5;
            		}
             		else{
                		delta = SAFETY*Math.pow(maxerr,-0.25);
                		if(delta < 0.1D) h *= 0.1;
                    		else h *= delta;
           		    }
        	}
        	return y;
    	}
}
