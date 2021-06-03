



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