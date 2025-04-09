#ifndef SOFTPLUS_REPARAM_H
#define SOFTPLUS_REPARAM_H

#include <stdio.h>
#include <assert.h>

typedef struct Params
{
	double beta, gamma, nu, tau, s, w, Mmax, epsilon, g, foi_slope;
} Params;

/*------------------------  MODEL PARAMETERS  -------------------------*/

extern const struct Params default_params;

extern const int dim;


double seasonality(double t, Params *params);

double foi(double H, double I, Params *params);

/* Numerical Integration */
int sirODE(double t, const double y[], double f[], void *params);

#endif
