#include <stdio.h>
#include <math.h>
#include <stdlib.h> // porque usa malloc

// flags compilacion: -lm (porque usamos math.h)

// probar si anda el metodo ones, pienso que necesita & en el param X

#define PI 3.14159265359
#define MAX_FILAS_COLUMNAS 400

void linspace (double a, double b, int c,double *d);
void meshgrid(double *x, double *y, double X[MAX_FILAS_COLUMNAS][MAX_FILAS_COLUMNAS], double Y[MAX_FILAS_COLUMNAS][MAX_FILAS_COLUMNAS], int nt);
void ones(double X[][1], int filas); // creo que ones deberia recibir las columnas para que quede mas generico y X ser [][] o ** o tipo void
void matriz_por_escalar(double k, double m[][1], int filas);
void TempSim(double *t, int T0, int T1, int T365, int nt, double *tmps, double *tmeds);
double sum(double *m, int size);

// jr
float brierei(double tmps, double *p);
double ** dot(double **v, double *u, int n);
double * cumtrapz(double *hs, int sz);
double ** transposef(double **m, int r, int c);
double ** ones2(int filas, int columnas);

void vonFoerster(double dt, double *t, double *tau, int nt, double *tmps, double *hmrs, double *pnu, double *pdes, double *pinput) {
    float tol = 0.0001;
    double T[nt][nt];
    double Tau[nt][nt];
    meshgrid(t,tau,T,Tau,nt); 

    double Pttau[nt][nt];
    double wts[nt][nt];
    double rates[nt]; // esto es un array, estaba como rates[][]
    double hmrsnul = sum(hmrs,nt);    

    if( hmrsnul == 0){
        for(int i = 0; i < nt; i++)
            rates[i] = brierei(tmps[i], pdes);
    }

    /* verifique los ultimos 9 y conciden.
    mis rates:
    rates:0.004572 - rates:0.000000 - rates:0.004204 - rates:0.076955 - rates:0.003840 
    rates:0.000000 - rates:0.003481 - rates:0.075761 - rates:0.003126
    rates del profe:
    4.57214372e-03   0.00000000e+00   4.20356871e-03   7.69550317e-02   3.83976782e-03
    0.00000000e+00   3.48073087e-03   7.57607866e-02   3.12644666e-03
    */

    /* solo para printear rates
    for(int i = 0; i < nt; i++) {
        if(i%5 == 0){
            printf("rates:%f \n",rates[i]);
        }else
            printf("rates:%f - ",rates[i]);
    }
    //*/

    printf( "cumtrapz /n");
    double * d = cumtrapz(rates, nt);

    /* verifique los ultimos 10 y coinciden (solo en python tienen mas precision pero creo que es porque printf redondea el numero)
    mis ultimos cumtrapz
    26.709788,    26.751158,   26.753444,    26.755546,    26.796126,    26.836523,   26.838443,    26.840183,    26.879804,    26.919248

    cumptraz del profe:
    26.70978843,  26.75115837, 26.75344444,  26.75554622,  26.79612552,  26.83652292, 26.83844281,  26.84018317,  26.87980393,  26.91924755

    */

    // para printear cumtrapz
    // for(int i = 0; i < nt; i++) printf("cumtrapz:%f\n", d[i]);



    double **X;
    printf( "ones2 /n");
    double ** x = ones2(nt, 1);
    for(int i = 0; i < nt; i++) printf(":%f\n", x[i][0]);
    for(int i = 0; i < nt; i++) x[i][0] = x[i][0] * dt;
    
    //printf("ones * dt \n");
    //for(int i = 0; i < nt; i++){
        //printf(":%f\n", x[i][0]); // muestra lo mismo que en python
    //}

    printf( "dot \n");
    double ** res = dot(x, d, nt);

    /* verifique 1ros y ultimos 3 y coinciden.
      mis dots:
      0.000000   0.019033    0.057958 , ...  6.710046 6.719951 6.729812
      dots del profe:
     0.          0.01903326  0.05795842...   6.71004579  6.71995098

     // para printear los datos de la funcion dots()
     for (int i = 0; i < nt; i++){
        for (int j = 0; j < nt; j++){
            printf( "%f", res[i][j] );
        }
        printf("\n");
    }    
    */

    double ** RTau = transposef(res, nt, nt);

    /* test transposef
    mis resultados:
    0.0000000, 0.0000000 ...
    0.019033,  0.019033, ...
    0.057958,  0.0579580, ...
    ...
    6.710046, 6.710046 ...
    6.719951, 6.719951 ...
    6.729812, 6.729812 ...

    RTau del profe
[[ 0.          0.          0.         ...,  0.          0.          0.        ]
 [ 0.01903326  0.01903326  0.01903326 ...,  0.01903326  0.01903326 0.01903326]
 [ 0.05795842  0.05795842  0.05795842 ...,  0.05795842  0.05795842 0.05795842]
 ..., 
 [ 6.71004579  6.71004579  6.71004579 ...,  6.71004579  6.71004579 6.71004579]
 [ 6.71995098  6.71995098  6.71995098 ...,  6.71995098  6.71995098 6.71995098]
 [ 6.72981189  6.72981189  6.72981189 ...,  6.72981189  6.72981189 6.72981189]]

    //printear transpose
    for (int i = 0; i < nt; i++){
        for (int j = 0; j < nt; j++){
            printf( "%f", RTau[i][j] );
        }
        printf("\n");
    }    

    */

    

}


void main(){    
	int tmin = 30;
	int tmax = tmin + 100; //130
	int nd = tmax - tmin; // 100
	int td = 4;
	int nt = nd * td; // 400
    double dt = nd/(double)nt; // 0.25	
	double rate = nd / nt;
    double t[nt];
    double tau[nt];
    linspace(tmin,tmax-dt,nt,t);
    linspace(tmin,tmax-dt,nt,tau);
    int T0   = 15;
    int T1   = 15;
    int T365 = 15; 
    double tmps[nt];
    double tmeds[nt];
    TempSim(t,T0,T1,T365,nt,tmps,tmeds);    
    
    double hmrs[nt];
    for(int i=0; i<nt;i++)
        hmrs[i] = 0;

    //fhrates = "BriereI";
    double pdes[3] ;
    pdes[0] =  0.000131734;
    pdes[1] = 10.25740308;  
    pdes[2] = 36.65400490;
    
    //fhnu    = "Varianza";
    double pnu[3];
    pnu[0]  = 0.0;
    pnu[1]  = 0.0;  
    pnu[2]  = 0.000223;

    double pinput[nt];
    for(int i=0; i<4;i++)
        pinput[i] = 25;
    for(int i=4; i<nt;i++)
        pinput[i] = 0;
    //idCorrida = "vFPy";
    double codfig = 8;
    printf("\n\n");
    vonFoerster(dt,t,tau,nt,tmps,hmrs,pnu,pdes,pinput);

    
}


void linspace (double a, double b, int c,double *d){ // funciona bien
    //printf("%f - %f - %i\n",a,b,c);    
    double aa = (double)b-a;
    //printf("%f\n\n",aa);
    double bb = (double)c-1;
    //printf("%f\n\n",bb);    
    double delta =aa/bb;
    //printf("%f\n\n",delta);
    for (int i=0; i<c; ++i){            
            d[i]=a + (i*delta);
            //printf("%lf\n",d[i]);
    }    
}

/*
Devuelve las coordenadas 2-D de la cuadrícula basadas en las coordenadas contenidas en los vectores x, y. 
X es una matriz donde cada fila es una copia de x, y Y es una matriz donde cada columna es una copia de y. 
La cuadrícula representada por las coordenadas X e Y tiene length(y)filas y length(x)columnas.
*/
void meshgrid(double *x, double *y, double X[MAX_FILAS_COLUMNAS][MAX_FILAS_COLUMNAS], double Y[MAX_FILAS_COLUMNAS][MAX_FILAS_COLUMNAS], int nt){
    for(int i=0; i<nt;i++){ // fila i
        for(int j=0;j<nt;j++){ // columna j
            X[i][j] = x[j]; 
        }
    }

    for(int i=0; i<nt;i++){ // columna i
        for(int j=0;j<nt;j++){ // fila j
            Y[j][i] = y[j];
        }
    }
}

/*
    primer parametro es una matriz de 'n' filas y una columna
    segundo parametro es la cantidad de filas que tiene la matriz
    
    simplemente llena de unos dicha matriz

    Aclaración: en el programa de python siempre q se llama a este metodo de la libreria NumPy lo hace con una matriz de una 
    columnas y 'nt' filas, es por esto q decidi q tome como parametro una matriz con esa cantidad de columnas, ademas no encontre
    como hacer q tome una matriz sin importar su tamaño, en todos los ejemplos que vi siempre se pone un numero al menos en la cantidad
    de columnas de la matriz.
*/
void ones(double X[][1], int filas){
    for(int i=0; i<filas;i++){ // fila i
        for(int j=0;j<1;j++){ // columna j
            X[i][j] = 1; 
        }
    }
}

double ** ones2( int filas, int columnas){
    
    double ** res = malloc(filas * sizeof(double *));
	for(int i = 0; i < filas; i++)
		res[i] = malloc(columnas * sizeof(double));

    for (int i = 0; i < filas; i++){
        for (int j = 0; j < columnas; j++){
            res[i][j] = 1;
            //printf(" %f ,", res[i][j]);
        }
        //printf("\n");
    }
    return res;
}


void matriz_por_escalar(double k, double m[][1], int filas){
     for(int i=0; i<filas;i++){ // fila i
        for(int j=0;j<1;j++){ // columna j
            m[i][j] = m[i][j] *  k;             
        }
    }
}

/*
    #TEMPSIM simula tabla de Temperaturas para hemisferio sur
    # t    : vector de tiempos en dias.fraccion desde t[0] instante inicial
    # T0   : temperatura media
    # T1   : amplitud termica anual
    # T365 : amplitud termica diaria
    # nt : tamaño de los vectores
    # tmps : vector a retornar
    # tmeds : otro vector a retornar
    */
void TempSim(double *t, int T0, int T1, int T365, int nt, double *tmps, double *tmeds){    
    T1 /= 2; // 7
    T365 /=2; // 7 
    for(int i=0; i<nt;i++){
        tmeds[i] = T0+T1*cos(2*PI/365*t[i]);
        if(i%5)
            printf("%f - ",tmeds[i]);
        else
            printf("\n%f - ",tmeds[i]);
    }
    printf("\n\n");      
     for(int i=0; i<nt;i++){
        tmps[i] = tmeds[i] - T365*cos(2*PI*(t[i]-(int)t[i]));
        if(i%5)
            printf("%f - ",tmps[i]);
        else
            printf("\n%f - ",tmps[i]);            
    }        
    printf("\n\n");    
}

/*
    para vector: suma cada uno de los valores de las casillas del vector que es enviado como parametro 
    junto con el tamaño del mismo y retorna el valor de dicha suma
*/
double sum(double *m, int size){
    double ret = 0;
    for(int i=0; i<size;i++){
        ret += (int) m[i];        
    }
    return ret;
}

// jr
float brierei(double tmpsi, double *p){
    float r = 0.0;
    if( tmpsi <= p[1] )
        r = 0.0;
    if( tmpsi > p[1] && tmpsi <= p[2] )
        r = p[0] * tmpsi * (tmpsi - p[1]) * sqrt(p[2] - tmpsi );
    else
        r = 0.0;
    
    return r;
}

/* pongo un double como matriz, porque lo que viene como parametro es una matriz de
n filas y 1 sola columna.
*/
double ** dot(double **v, double *u, int n)
{
    double **res = malloc(n * sizeof(double *));
	for(int i = 0; i < n; i++)
		res[i] = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            res[j][i] = v[i][0]*u[i];
            //printf(" %f ,", res[i][j]);
        }
        //printf("\n");
        //if(i > 5) break;
    }

    return res;
}


double * cumtrapz(double *hs, int sz)
{
    double *ct = malloc(sz * sizeof(double *));
    for(int i = 0; i < sz; i++) ct[i] = 0;

    for(int i = 1; i < sz; i++)
        ct[i] = ct[i - 1] + ( hs[i - 1] + hs[i] ) / 2.0; 

    return ct;
}

double ** transposef(double **m, int r, int c)
{
    double **transpose = malloc(c * sizeof(double *));
	for(int i = 0; i < c; i++)
		transpose[i] = malloc(r * sizeof(double));

    for(int i=0; i < r; i++)
        for(int j=0; j < c; j++)
            transpose[j][i] = m[i][j];
    

    return transpose;
}