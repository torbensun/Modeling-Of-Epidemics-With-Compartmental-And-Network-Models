#ifndef MYNUMERICS_H
#define MYNUMERICS_H


// Function that returns array filled with commuters from a cell
double *commutersFrom(double **commuters, int dimension, int l);
// Function that returns array filled with commuters to a cell
double *commutersTo(double **commuters, int dimension, int l);
//
double *effective_population(double **commuters, double *population, int dimension);

#endif