#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<tgmath.h>


#include"lib.h"

// FUNTION INITIALIZATION


int function_of_system(double t, double y[], double f[], double **commuters, double *population);

void read_commuters(double **commu, int N);

void fill_pop_array();



// initializing parameters. For more explanation of these see the Latex/PDF 

const double alpha = 0.2; // rate of infected
const double beta = 1.0/14.0; // recovery rate
const int dimension = 12; // dimension of the system
const double p = 0.0264; // death rate
// TODO find a source or a better estimate for the value (in case it isn't good enough already)
const double t0 = 10.0/24.0; // ratio of "commuters not home" (in hours) to hours a day.


// ####################################################################
// #                            MAIN                                  #
// ####################################################################

int main(){
    // size of matrix
    int N = 38;
    // creating matrix for commuters. Double pointer 
    double **commu = (double**)malloc(sizeof(double) * N);
    // set the double pointer on mallocs
    for (int i = 0; i < N; i++){
        commu[i] = (double*) malloc(sizeof(double) * N);
    }


    // filling the population vector
    double population[dimension] = {2.758170000000000000e+05, 7.045800000000000000e+04, 1.362920000000000000e+05, 1.402510000000000000e+05, 1.322850000000000000e+05, 3.260410000000000000e+05, 2.133100000000000000e+05, 2.367640000000000000e+05, 2.021370000000000000e+05, 1.006290000000000000e+05, 1.000060000000000000e+05, 8.341600000000000000e+04};// array with population of system
    //fill_pop_array(population);

    

    // filling the commuter matrix
    read_commuters(commu, N);


    // solving the differential equation using rk4
    // the time period is from July 24. to November 1. 

    // time of starting. Note that the unit is days
    double t_start = 0.0;
    // The simmulation spans 100 days
    double t_end = 100.0;

    // choosing a delta_t
    // TODO maybe adjust this for better resolution/results
    double delta_t = 0.001; // this time step corresponds to about 1.4 minutes (someone please check the math on this)

    // dimension of differential equation (as there are 4 differential equations per cell)
    int dim_deq = 4 * dimension;

    // array with initial conditions
    // TODO import initial conditions of the system
    double y[dim_deq] = {9.984119905589575739e-01, 6.526066196064781778e-05, 1.490118448101458628e-03, 3.263033098032390889e-05, 9.984955576371739028e-01, 0.000000000000000000e+00, 1.405092395469641556e-03, 9.934996735643929645e-05, 9.980703196079007133e-01, 2.201156340797699005e-05, 1.702227570216887330e-03, 2.054412584744519116e-04, 9.973119621250472466e-01, 3.565036969433373228e-05, 2.516916100419961527e-03, 1.354714048384681867e-04, 9.988585251540235133e-01, 1.511887213213894254e-05, 1.088558793514003958e-03, 3.779718033034735889e-05, 9.960127713999159527e-01, 3.373808815455725004e-05, 3.711189697001297545e-03, 2.423008149281838885e-04, 9.987764286718859852e-01, 0.000000000000000000e+00, 1.200131264357038978e-03, 2.344006375697341754e-05, 9.981922927472082208e-01, 6.757784122586203419e-05, 1.609197344190839820e-03, 1.309320673751076896e-04, 9.980211440755527574e-01, 6.431281754453662964e-05, 1.870018848602680297e-03, 4.452425830006381531e-05, 9.976746265986942142e-01, 9.937493167973447138e-06, 2.156436017450238159e-03, 1.589998906875751542e-04, 9.983900965942043015e-01, -9.999400035997840970e-06, 1.499910005399676064e-03, 1.199928004319740781e-04, 9.991848086697995290e-01, 2.397621559413062403e-05, 7.672388990121799691e-04, 2.397621559413062403e-05, 9.975464841437541308e-01, 4.407513514214230984e-05, 2.300549210752604427e-03, 1.088915103511751232e-04};


    // file for saving the solution of calculations
    FILE *sol = fopen("DiffEq.txt", "w");

    // loop to solve the differential equation
    while (t_start < t_end){
        // applying the rk4 step function to calculate the next step in y
        rk4_step(t_start, delta_t, y, function_of_system, dim_deq, commu, population);

        // adding time to the file
        fprintf(sol, "&lf,", t_start);

        // adding all the y values to file
        for (int i = 0; i < dim_deq; i++){
            fprintf(sol, "%lf,", y[i]);
        }

        // new line in file
        fprintf(sol, "\n");


        // increasing time by delta_t
        t_start += delta_t;

    }


    // free the mallocs
    for (int i = 0; i < N; i++){
        free(commu[i]);
    }
    free(commu);

    // closing file
    fclose(sol);

    return 0;

}


// ####################################################################
// #                            FUNCTIONS                             #
// ####################################################################


// ################################### COMMUTERS ###################################

/**
 * @brief function to read the txt file with the commuters
 * 
 * @param com Array to fill with the commuters
 * @param N The size of the matrix that will be read. 
 */
void read_commuters(double **commu, int N){

    // Commuters
    // loading txt with commuters
    FILE *com = fopen("Pendler.txt", "r");

    // filling array
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            // Break condition if it can't read in value
            if(!fscanf(com, "%lf\t", &commu[i][j])){
                printf("The matrix didn’t read properly");
                break;
            } 
        }
    }

    fclose(com);
    return;
}

// ##################################### FUNCTION OF SYSTEM ######################################################


/**
 * @brief The function of the system. The system is in this case is modeled by the approach of constant coefficients, i.e. the weighted mean of the standard
 *        and the dynamic model
 * 
 * @param t time at which it is to be calculated
 * @param y array contaning the system quantities [S,...S, I,...I, R,...R, D,...D]
 * @param f array which is to be filled with derivatives
 * @param commuters matrix with commuters
 * @param population array with population
 * @return int 
 */
int function_of_system(double t, double y[], double f[], double **commuters, double *population){
    // arrays for different quantities
    // effective infected
    double *Ieff = (double*)malloc(sizeof(double) * dimension);
    // S
    double *So = (double*)malloc(sizeof(double) * dimension);
    // I
    double *Io = (double*)malloc(sizeof(double) * dimension);
    // R
    double *Ro = (double*)malloc(sizeof(double) * dimension);
    // D
    double *Do = (double*)malloc(sizeof(double) * dimension);
    // effective population
    double *Neff = (double*)malloc(sizeof(double) * dimension);
    
    // filling the arrays
    for (int i = 0; i < dimension; i++){
        So[i] = y[i];
        Io[i] = y[i + dimension];
        Ro[i] = y[i + 2*dimension];
        Do[i] = y[i + 3*dimension];
        Neff[i] = effective_population(commuters, population, dimension)[i];
    }

    // fill the effective infected
    for (int i = 0; i < dimension; i++){
        Ieff[i] = effective_infected(commuters, population, dimension, Io)[i];
    }

    //TODO finish the function
    for (int i = 0; i < dimension; i++){
        // sum for last term in derivative. For detail see PDF/LaTeX
        double sum = 0;
        for (int k = 0; k < dimension; k++){
            sum += commutersFrom(commuters, dimension, i)[k] * Ieff[k];
        }

        // dSdt
        f[i] = - (1 - t0) * alpha * So[i] * Io[i] - t0/population[i] * So[i] * (Neff[i] * Ieff[i] - sum);

        // dIdt
        f[i + dimension] = - f[i] - beta*Io[i];

        // dRdt
        f[i + 2 * dimension] = beta * (1 - p) * Io[i];

        // dDdt
        f[i + 3 * dimension] = beta * p * Io[i];
    }


    // freeing the malloc
    free(Ieff);
    free(So);
    free(Io);
    free(Ro);
    free(Do);
    free(Neff);
    return 0;
}

// ######################################### POPULATION ###################################################
/**
 * @brief Function to read in the population from txt file and save in the array "population"
 * 
 * @param population Array to save the population in
 */
void fill_pop_array(double* population)
{
    FILE* popdata = fopen("Internal Data/popdata38.txt", "r");
    for (int j = 0; j < dimension; ++j)
        fscanf(popdata, "%lf", &population[j]);
    fclose(popdata);
}

