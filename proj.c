#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<tgmath.h>


#include"lib.h"

// FUNTION INITIALIZATION


int function_of_system(double t, double y[], double f[], double **commuters, double *population);

void read_commuters(double **commu, int N);



// initializing parameters. For more explanation of these see the Latex/PDF 

const double alpha = 0.2; // rate of infected
const double beta = 1.0/14.0; // recovery rate
const int dimension = 12; // dimension of the system
const double p = 0.0264; // death rate
// TODO fill this or make code read it.
const double population[dimension] = {}; // array with population of system
    




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


    // freeing the malloc
    free(Ieff);
    free(So);
    free(Io);
    free(Ro);
    free(Do);
    free(Neff);
    return 0;
}
 
// TODO: (either here or in lib.c) write a rk4 solver or adjust the one from my_numerics


