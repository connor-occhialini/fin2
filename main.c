#include <stdio.h>
#include <gsl/gsl_math.h>

#include <gsl/gsl_monte.h>		//Relevant GSL functions for VEGAS numerical integration
#include <gsl/gsl_monte_vegas.h>

#include "timer.c"			//Timing Functions
#include "plotfun.h"			//Automatic Plotting Function

void timer_start (void);
double timer_stop (void);

void pfun(double x[], double y[], double u[], double v[], int points);

extern double g (double *t, size_t dim, void *params);
struct my_params {double r;};

double f (double x1, double x2, double y1, double y2, double z1, double z2, double r);

double dipole_approx (double r);

int main (void)
{
	int np = 20;
	double distance[20];
	double vegas_energy[20];
	double res, err;
        double error1[20];
	double error2[20];
	

	size_t dim = 6;
	double xl[] = {0., 0., 0., 0., 0., 0.};
	double xu[] = {1., 1., 1., 1., 1., 1.};
	
	gsl_rng *r = gsl_rng_alloc (gsl_rng_taus2);
	unsigned long seed = 1UL;
	
	gsl_rng_set (r, seed);
	
	size_t calls = 1000000;

	double dd = (4.-1.001)/np;
	

	timer_start();

	for(int q = 0; q < np; q++)
	{
		distance[q] = 1.001 + q*dd;
	}
	
	//VEGAS Integration
	
	for(int j = 0; j < np; j++)
	{
		struct my_params params = {distance[j]};
		gsl_monte_function G = { &g, dim , &params};

		gsl_monte_vegas_state *sv = gsl_monte_vegas_alloc (dim);
	
		gsl_monte_vegas_init (sv);

		gsl_monte_vegas_integrate (&G, xl, xu, dim, calls / 10, r, sv, &res, &err);

		do
		{
			gsl_monte_vegas_integrate (&G, xl, xu, dim, calls, r, sv, &res, &err);
			fflush (stdout);
		}
		while (fabs (gsl_monte_vegas_chisq (sv) - 1.0) > 0.2);
		vegas_energy[j] = fabs(res);
		error1[j] = fabs(err);

		gsl_monte_vegas_free (sv);
	}

	double t1 = timer_stop();

	gsl_rng_free (r);

	int iterations = 1000000;
	int tot = 20;
	double x1, x2, y1, y2, z1, z2;
	double handmade_energy[20];

	gsl_rng *zed = gsl_rng_alloc (gsl_rng_taus2);

	gsl_rng_set (zed, seed);

	timer_start();

	//Handmade Monte Carlo Integrator

	for(int j = 0; j < 20.; j++)
	{
	    double d = 1.001 + j * dd;
	    double ival = 0;
	    double i2val = 0;   
	    for(int z = 0; z < 16; z++)
	    {
	      double value = 0;
       		for (int q = 0; q < iterations; q++)
        	{
            		x1 = gsl_rng_uniform(zed);
            		x2 = gsl_rng_uniform(zed);
            		y1 = gsl_rng_uniform(zed);
            		y2 = gsl_rng_uniform(zed);
            		z1 = gsl_rng_uniform(zed);
            		z2 = gsl_rng_uniform(zed);
            		value += f(x1, x2, y1, y2, z1, z2, d);
       		}
		ival += value;
		i2val += value*value;
	    }

	    double result = ival/iterations/16;
	    double result2 = i2val/iterations/iterations/16;

	    error2[j] = sqrt(fabs((result*result)-result2));
	    
	    handmade_energy[j] = fabs(result);
	}

	double t2 = timer_stop();

	//Dipole Approximation

	double dipole_energy[20];
	for(int u = 0; u < tot; u++)
	{
		dipole_energy[u] = dipole_approx(distance[u]);
	}

	printf("Distance   VEGAS      Error      Handmade   Error      Dipole\n");
		
	for(int w = 0; w < tot; w++)
	{
	  printf("%.6f   %.6f   %.6f   %.6f   %.6f   %.6f\n", distance[w], vegas_energy[w], error1[w], handmade_energy[w], error2[w], dipole_energy[w]);
	}

	printf("Time for VEGAS:  %.6f\nTime for HANDMADE:  %.6f\n", t1, t2);	

	//Plots the 3 methods vs. distance automatically (function in plotfun.h)

	pfun(distance, vegas_energy, handmade_energy, dipole_energy, 20);
	
	return 0;
}

double g (double *t, size_t dim, void *params)
{
    	double dist2, delta;
    	double rho1, rho2;
    	int sdim = ((int) dim)/2;
    	double r = *((double *) params);

    	dist2 = 0.;
    	for (int i = 0; i < sdim; i++)
    	{
        	delta = t[i] - t[i + sdim]; 
        	if (i == 0)
        	{
           		delta += r;
        	}
        	dist2 += delta*delta; 
    	}

    	double norm = 512./(9.*(M_PI - 2.));
    	rho1 = norm * atan(2*(t[0]-.5))*pow(sin(M_PI*t[1]),4.)*pow(sin(M_PI*t[2]),4.);
    	rho2 = norm * atan(2*(t[3]-.5))*pow(sin(M_PI*t[4]),4.)*pow(sin(M_PI*t[5]),4.);

    	return rho1*rho2/sqrt(dist2);
}

double f (double x1, double x2, double y1, double y2, double z1, double z2, double r)
{
	double dist2;
    	double rho1, rho2;
    
    	dist2 = (x1 - x2 + r)*(x1 - x2 + r) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2);

    	double norm = 512./(9.*(M_PI - 2.));
    	rho1 = norm * atan(2*(x1-.5))*pow(sin(M_PI*y1),4.)*pow(sin(M_PI*z1),4.);
    	rho2 = norm * atan(2*(x2-.5))*pow(sin(M_PI*y2),4.)*pow(sin(M_PI*z2),4.);

    	return rho1*rho2/sqrt(dist2);
}

double dipole_approx (double r)
{
	return 2/(r*r*r);
}

