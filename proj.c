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
    double population[dimension]; // array with population of system
    fill_pop_array(population);

    for (int j = 0; j < dimension; ++j)
        printf("%f\n",population[j]);

    // filling the commuter matrix
    read_commuters(commu, N);


    //void *params = commu;



    // free the mallocs
    for (int i = 0; i < N; i++){
        free(commu[i]);
    }
    free(commu);
    return 0;

}


// ####################################################################
// #                            FUNCTIONS                             #
// ####################################################################



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

    return;
}

// ###########################################################################################


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

// ###########################################################################################

void fill_pop_array(double* population)
{
    FILE* popdata = fopen("Internal Data/popdata38.txt", "r");
    for (int j = 0; j < dimension; ++j)
        fscanf(popdata, "%lf", &population[j]);
    fclose(popdata);
}

// TODO: (either here or in lib.c) write a rk4 solver or adjust the one from my_numerics


