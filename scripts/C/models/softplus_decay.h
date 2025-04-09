#ifndef SOFTPLUS_DECAY_H
#define SOFTPLUS_DECAY_H

#include <stdio.h>
#include <assert.h>

typedef struct Params
{
	double beta, gamma, nu, tau, Hthres, s, w, Mmax, epsilon, delta;
} Params;

/*------------------------  MODEL PARAMETERS  -------------------------*/

extern const struct Params default_params;

extern const int dim;


double seasonality(double t, Params *params);

double foi(double H, double I, Params *params);

/* Numerical Integration */
int sirODE(double t, const double y[], double f[], void *params);

#endif
