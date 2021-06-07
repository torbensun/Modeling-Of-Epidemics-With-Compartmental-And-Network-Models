#ifndef MYNUMERICS_H
#define MYNUMERICS_H


// Function that returns array filled with commuters from a cell
double *commutersFrom(double **commuters, int dimension, int l);
// Function that returns array filled with commuters to a cell
double *commutersTo(double **commuters, int dimension, int l);
// Funtion to calculate the effective population in every cell
double *effective_population(double **commuters, double *population, int dimension);
// Function to calculate the effective infected in every cell
double *effective_infected(double **commuters, double *population, int dimension, double *infected);

// typedef for function to be integrated
typedef
int ode_func(double, double[], double[], double**, double*);
// Runge Kutta 4 for this system
void rk4_step(double t, double delta_t, double y[], ode_func func, int dim, double **commuters, double *population);



#endif