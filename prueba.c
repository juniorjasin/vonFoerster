#include <stdio.h>
#include <math.h>
#include <stdlib.h> // porque usa malloc
/*
"""
lo que recibe como parametro es un float y un arreglo de tamanio 3 con tres valores float
"""
def BriereI(tmpsi, p):
# BRIEREI - BriereI - 3 par�metros
	r  = 0.0;
	if tmpsi <= p[1]:
	   r = 0.0;
	elif (tmpsi > p[1]) and (tmpsi <= p[2]):
	   r = p[0] * tmpsi * (tmpsi - p[1]) * sqrt(p[2] - tmpsi);
	else:
	   r = 0.0;
	return r;
*/

float brierei(float tmpsi, float p[]);

// calcular el producto punto

/*
params: 
*/

//double dot(double *v, double *u, int n);


/* pongo un double como matriz, porque lo que viene como parametro es una matriz de
n filas y 1 sola columna.

nota: habria que saber de antemano el valor
*/


/*
def cumtrapz(hs): # hs es una lista/vector
    ct = np.zeros(len(hs));
    for i in range(1,len(hs)):
        ct[i] = ct[i-1] + (hs[i-1] + hs[i]) / 2.0;
    return ct;
*/

double * cumtrapz(double *hs, int sz){
    double *ct = malloc(sz * sizeof(double *));
    for(int i = 0; i < sz; i++) ct[i] = 0;

    for(int i = 1; i < sz; i++)
        ct[i] = ct[i - 1] + ( hs[i - 1] + hs[i] ) / 2.0; 

    return ct;
}


double ** transposef(double **m, int r, int c)
{
    //printf("viene matriz de %i filas y %i columnas", r, c);
    double **transpose = malloc(c * sizeof(double *));
	for(int i = 0; i < c; i++)
		transpose[i] = malloc(r * sizeof(double));

    //printf("se creo matriz de %i filas y %i columnas", c, r);

    
    printf("arranca filas:%i - columnas:%i \n", c, r);
    for(int i=0; i < r; i++){
        //printf("vuelta %i", i);
        for(int j=0; j < c; j++){
            //printf("vuelta %i", j);
            //printf("%f ", m[i][j]);
            //printf("%lf ", transpose[j][i]);
            transpose[j][i] = m[i][j];
        }
        printf("\n");
    }

    return transpose;
}


double ** dot(double **v, double *u, int n)
{
    double **res = malloc(n * sizeof(double *));
	for(int i = 0; i < n; i++)
		res[i] = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            res[i][j] = v[i][0]*u[i];
            //printf(" %f ,", res[i][j]);
        }
        //printf("\n");
    }

    return res;
}

void main(){ 

    /* ------------------------------------------------------------------------------------------ 
    //test brierei
    float tmpsi = 10.1;
    float p[] = {1,2,30};
    printf("array: %f - %f - %f ", p[0], p[1], p[2]);
    printf("float: %f", tmpsi);
    float r = brierei(tmpsi, p);
    printf("result: %f", r);
    //*/

    //* ------------------------------------------------------------------------------------------
    // test dot    
    int nrows = 5;
    int ncolumns = 1;

    // creo array y lo lleno con algun valor
    double *array = malloc(nrows * sizeof(double *));
    
    for(int i = 0; i < nrows; i++) array[i] = i;

    // creo matriz y lo lleno con algun valor
    double **array1 = malloc(nrows * sizeof(double *));
	for(int i = 0; i < nrows; i++)
		array1[i] = malloc(ncolumns * sizeof(double));

    for(int i = 0; i < nrows; i++) array1[i][0] = i;

    for(int i = 0; i < nrows; i++){
        printf("val:%f \n", array1[i][0]);
    }

    printf("\n se llama a dot() \n");
    dot(array1, array, nrows);

    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
         //   printf( "val:%f /n", res[i][j] );
        }
    }

    //*/

    /*------------------------------------------------------------------------------------------
    // test cumtrapz
    int nrows = 5;
    double *array = malloc(nrows * sizeof(double *));
    for(int i = 0; i < nrows; i++) array[i] = i;

    double *res = cumtrapz(array, nrows);
    printf("res:%f", *res);
    */

    /* ------------------------------------------------------------------------------------------
    // test transpose
    printf("test transpose");
    int nrows = 5;
    int ncolumns = 3;
    double **array1 = malloc(nrows * sizeof(double *));
	for(int i = 0; i < nrows; i++)
		array1[i] = malloc(ncolumns * sizeof(double));
    printf("se creo matriz \n");

    for(int i=0; i < ncolumns; ++i){
        for(int j=0; j < nrows; ++j){
            printf("%f", array1[j][i] );
        }
        printf("\n");
    }

    
    double ** t = transposef(array1, nrows, ncolumns);

    printf("se llama a transpose() \n");
    for(int i=0; i < nrows; ++i){
        for(int j=0; j < ncolumns; ++j){
            printf("%f", t[j][i] );
        }
        printf("\n");
    }

    */


}


float brierei(float tmpsi, float p[]){
    float r = 0.0;
    if( tmpsi <= p[1] )
        r = 0.0;
    if( tmpsi > p[1] && tmpsi <= p[2] )
        r = p[0] * tmpsi * (tmpsi - p[1]) * sqrt(p[2] - tmpsi );
    else
        r = 0.0;
    
    return r;
}