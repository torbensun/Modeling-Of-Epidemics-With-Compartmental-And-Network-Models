#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<tgmath.h>


#include"lib.h"

// FUNTION INITIALIZATION


int function_of_system(double t, double y[], double f[], double *params);

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


int function_of_system(double t, double y[], double f[], double *params){
    




    return 0;
}
 



