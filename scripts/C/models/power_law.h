#ifndef POWER_LAW_H
#define POWER_LAW_H

#include <stdio.h>
#include <assert.h>

typedef struct Params
{
	double beta, gamma, nu, tau, alpha, s, w;
} Params;

/*------------------------  MODEL PARAMETERS  -------------------------*/

extern const struct Params default_params;

extern const int dim;


double seasonality(double t, Params *params);

double foi(double H, double I, Params *params);

int sirODE(double t, const double y[], double f[], void *params);

#endif