



/**
 * @brief Function to calculate the commuters going out of cell l
 * 
 * @param commuters The matrix with all the commuters
 * @param dimension Integer telling how big the system is
 * @param l         Integer saying from what cell the commuters come from
 * @return double*  Array with the commuters from cell l
 */
double *commutersFrom(double **commuters, int dimension, int l){
    
    // array to be returned
    double a[dimension];

    // fill
    for (int i = 0; i < dimension; i++){
        a[i] = commuters[i][l];
    }

    return a;
}





/**
 * @brief Function to calculate the commuters going to cell l
 * 
 * @param commuters The matrix with all the commuters
 * @param dimension Integer telling how big the system is
 * @param l         Integer saying to what cell the commuters are going
 * @return double*  Array with the commuters going to cell l
 */
double *commutersTo(double **commuters, int dimension, int l){

     // array to be returned
    double a[dimension];

    // fill
    for (int i = 0; i < dimension; i++){
        a[i] = commuters[l][i];
    }

    return a;

}



/**
 * @brief function to calculate the effective population in a cell. It is defined as initial population subtracted with the commuters commuting out of the cell.
 * 
 * @param commuters Matrix (38 * 38) filled with commuters
 * @param population Array with population
 * @return double*  Array with effective population
 */

double *effective_population(double **commuters, double *population, int dimension){

    // initialize return array
    double array[dimension];
    
    // iterating over every cell and calculating effective population for every cell
    for (int i = 0; i < dimension; i++){

        // setting the value of return array to the population of the respective cell
        array[i] = population[i];

        // iterating over commuters and subtracting them from the population
        for (int k = 0; k < dimension; k++){
            array[i] -= commuters[k][i];
        }

    }

    // returning array
    return array;
}


/**
 * @brief Function to calculate the effective infected for every cell in the system
 * 
 * @param commuters   Matrix with commuters
 * @param population  Array with population
 * @param dimension   Dimension of the system
 * @param infected    Array with the infected
 * @return double*    Array with effective infected
 */

double *effective_infected(double **commuters, double *population, int dimension, double *infected){

    // initialize return array
    double Ieff[dimension];

    // iterating over every cell to fill return array and calculating the effective infected
    for (int i = 0; i < dimension; i++){

        // setting array values
        Ieff[i] = infected[i];

        // setting two help variables 
        double cfrom = 0;
        double cto = 0;

        // loop for adjusting parameters
        for (int k = 0; k < dimension; k++){
            cfrom -= commuters[k][i];
            cto += commuters[i][k]*infected[k];
        }

        // calculating the effective infected
        Ieff[i] += 1.0/population[i]*(cfrom*infected[i] + cto);
    }

    return Ieff;
    
}





