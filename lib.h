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
int ode_func(double, const double[], double[], void*);
// Runge Kutta 4 for this system
void rk4_step(double t, double delta_t, double y[], ode_func func, int dim, double **commuters, double *population);


//----------------------- LINEAR ALGEBRA --------------------------/*

/* Gauss elimination algorithm that transforms a nxm matrix into a triangular matrix */
int gauss(int n, int m, double A[n*m], int pivoting);
/* Returns the determinant of a nxn matrix A */
double det(int n, double A[n*n], int pivoting);
/* Calculates the solution x of linear equation Ax=b */
void lgs_solve(int n, double A[n*n],double b[n], double x[n], int pivoting);
/* Calculates the inverse of nxn matrix A */
void inverse(int n, double A[n*n], double A_inverse[n*n]);
/* Calculates the matrix product of nxm matrix A and mxp matrix B */
void matrix_product(int n, int m, int p, double A[n*m], double B[m*p], double AB[n*p]);
/* Calculates the transpose of nxm matrix A */
void transpose(int n, int m, double A[n*m], double A_transpose[m*n]);


#endif