#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include "softplus.h"

/*------------------------  MODEL PARAMETERS  -------------------------*/

const struct Params default_params = {0.5, 0.1, 0.01, 20., .001, 0.25, 2.*M_PI/360., .9, 0.00025, 1.0};

const int dim = 5;

/*------------------------  MODEL FUNCTIONS  -------------------------*/

double seasonality(double t, Params *params){
	return 1.+params->s*cos(params->w*t);
}

double foi(double H, double I, Params *params){
	double slope = params->Mmax/params->Hthres;
	return (slope*params->epsilon*log(exp(1/params->epsilon*(params->Hthres-H))+1)+1.-params->Mmax)*I*params->beta;
}

int sirODE(double t, const double y[], double f[], void *params){
	for (int i = 0; i<dim; i++)
		f[i] = 0;
	
	struct Params *my_params = params;

	
	f[0] += -y[0]*foi(y[4],y[1],my_params)*seasonality(t,my_params) + my_params->nu*y[2];
	f[1] += y[0]*foi(y[4],y[1],my_params)*seasonality(t,my_params) - my_params->gamma*y[1];
	f[2] += my_params->gamma*y[1] -my_params->nu*y[2];
	f[3] += (my_params->g*y[1]-y[3])/my_params->tau;
	f[4] += (y[3]-y[4])/my_params->tau;
	
	return GSL_SUCCESS;
}
