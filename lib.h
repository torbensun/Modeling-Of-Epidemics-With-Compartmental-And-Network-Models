#ifndef MYNUMERICS_H
#define MYNUMERICS_H


// Function that returns array filled with commuters from a cell
void commutersFrom(double **commuters, int dimension, int l, double *a);
// Function that returns array filled with commuters to a cell
void commutersTo(double **commuters, int dimension, int l, double *a);
// Funtion to calculate the effective population in every cell
void effective_population(double **commuters, double *population, int dimension, double *array);
// Function to calculate the effective infected in every cell
void effective_infected(double **commuters, double *population, int dimension, double *infected, double *Ieff);

// typedef for function to be integrated
typedef
int ode_func(double, double[], double[], double**, double*);
// Runge Kutta 4 for this system
void rk4_step(double t, double delta_t, double y[], ode_func func, int dim, double **commuters, double *population);



#endif