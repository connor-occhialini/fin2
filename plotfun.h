#include <stdlib.h>
#include <stdio.h>
#define COMMANDS 2

//takes arguments from array structures and outputs a GNUPLOT of the data with connecting lines.
void pfun(double x[], double y[], double u[], double v[], int points)
{
    char * commandsForGnuplot[] = {"set terminal png","set output 'energy.png'","set title \"Integral Values as a Function of Distance\"","set logscale x 2", "set logscale y 10", "plot 'data.temp' u 1:2 title 'VEGAS' w l, 'data.temp' u 1:3 title 'Handmade Integrator' w l, 'data.temp' u 1:4 title 'Dipole Approximation' w l"};
    FILE * temp = fopen("data.temp", "w"); // write coordinates here.
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;
    int q;
	for(q = 0; q < points ;q++)
	{
		fprintf(temp, "%lf %lf %lf %lf\n", x[q], y[q], u[q], v[q]);
	}

 
    for (i=0; i < 6; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
}

