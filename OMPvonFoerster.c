#include <stdio.h>
#include <math.h>
#include <stdlib.h> // porque usa malloc
#include <sys/time.h>
#include <omp.h>

/* PREGUNTAS
    1) Cada seccion en parallelSections se apodera de un solo hilo?
    2) en los for que hacen inicializaciones vale la pena paralelizar?
    3) en el caso de matrices usamos solo un #pragma omp for,  o hay una manera de usar uno por cada for y que sea eficiente
    4) en el metodo trapz, el ultimo for decidimos no paralelizarlo xq haciendo pruebas por separado, obtuvimos q secuencialemte era mas rapido, 
        esta de acuerdo el profe?. Intentamos paralelizar el for de afuera y luego el for interno, los resultados fueron q el for externo paralelizado
        es mas lento que el interno paralelizado, pero secuencialemte termina siendo muy similar y en casos menor a el for interno paralelizado
    5) preguntar sobre los tres for anidados (teniendo en cuenta lo q sucesio en el punto anterior) en el metodo multiplyMatrices que realiza prod. de matrices
*/

/*
    1) cambiar paralell de transpose y dot
*/

struct timeval start, end;
int chunk = 1;
// flags compilacion: -lm (porque usamos math.h)

// probar si anda el metodo ones, pienso que necesita & en el param X

#define PI 3.14159265359
#define MAX_FILAS_COLUMNAS 400

// andy
void linspace (double a, double b, int c,double *d);
void meshgrid(double *x, double *y, double **X, double **Y, int nt);
void ones(double X[][1], int filas); // creo que ones deberia recibir las columnas para que quede mas generico y X ser [][] o ** o tipo void
void matriz_por_escalar(double k, double m[][1], int filas);
void TempSim(double *t, int T0, int T1, int T365, int nt, double *tmps, double *tmeds);
double sum(double *m, int size);
double * trapz(double **v, double filas, double columnas);

// jr
float brierei(double tmps, double *p);
double ** dot(double **v, double *u, int n);
double * cumtrapz(double *hs, int sz);
double ** transposef(double **m, int r, int c);
double ** ones2(int filas, int columnas);
double ** diff(double ** a, double **b, int f, int c); 
void absMatrix(double **m, int f, int c);
void sumEscalarToMatrix(double **m, int f, int c, double val);
void escalarMatrixMultiplication(double **m, int f, int c, double val);
double ** multiplyMatrices(double ** firstMatrix, double ** secondMatrix, int rowFirst, int columnFirst, int rowSecond, int columnSecond);


void vonFoerster(double dt, double *t, double *tau, int nt, double *tmps, double *hmrs, double *pnu, double *pdes, double *pinput) {
    float tol = 0.0001;
    double **T;
    double **Tau;
    double **Pttau;
    for (int i=0; i<nt; ++i){ 
            //printf("\n soy hilo: %3d (vuelta:%3d) -> Calculado:  d[i]=a + (i*delta) =  %lf  \n", omp_get_thread_num(), i, d[i]);
            //printf("2 -> %lf\n",tmps[i]);
    }

    T = malloc(nt * sizeof(double *));
	for(int i = 0; i < nt; i++)
		T[i] = malloc(nt * sizeof(double));

    Tau = malloc(nt * sizeof(double *));
	for(int i = 0; i < nt; i++)
		Tau[i] = malloc(nt * sizeof(double));

    Pttau = malloc(nt * sizeof(double *));
	for(int i = 0; i < nt; i++)
		Pttau[i] = malloc(nt * sizeof(double));

    meshgrid(t,tau,T,Tau,nt); 
    
    double rates[nt]; // esto es un array, estaba como rates[][]
    double hmrsnul = sum(hmrs,nt);
    // for (int i=0; i<nt; ++i){            
    //         //printf("\n soy hilo: %3d (vuelta:%3d) -> Calculado:  d[i]=a + (i*delta) =  %lf  \n", omp_get_thread_num(), i, d[i]);
    //         printf("%lf\n",tmps[i]);
    // }

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

    // printf( "cumtrapz \n");
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
    // printf( "ones2 \n");
    double ** x = ones2(nt, 1);
    //for(int i = 0; i < nt; i++) printf(":%f\n", x[i][0]);
    for(int i = 0; i < nt; i++) x[i][0] = x[i][0] * dt;
    
    //printf("ones * dt \n");
    //for(int i = 0; i < nt; i++){
        //printf(":%f\n", x[i][0]); // muestra lo mismo que en python
    //}

    // printf( "dot \n");
    double ** RT = dot(x, d, nt);

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

    double ** RTau = transposef(RT, nt, nt);
    if(pnu[0]==0 && pnu[1]==0){
        double nu = pnu[2];

        // voy a hacer RT-RTau        
         for (int i = 0; i < nt; i++){
            for (int j = 0; j < nt; j++){
                if(T[i][j] > Tau[i][j]){    
                    Pttau[i][j] = (exp(-pow((1-(RT[i][j] - RTau[i][j])),2)/(4*nu*(fabs(T[i][j]-Tau[i][j])+tol))))/sqrt(4*PI*nu*(pow(fabs(T[i][j]-Tau[i][j]),3)+tol));
                    //printf( "%e | ",Pttau[i][j]);
                }else
                {
                    Pttau[i][j] = 0;
                    //printf( "%e | ",Pttau[i][j]);
                }
            }                        
            //printf("\n\n");
        }

        /*
        for (int i = 0; i < 3; i++){
                for (int j = 0; j < nt; j++){
                    printf( "%e | ",Pttau[i][j]);
                }
                printf("\n\n");
        }
        */
    }


    //ints = dt*np.trapz(Pttau,axis=1);

    double * ints = trapz(Pttau, nt, nt);
    for (int i = 0; i < nt; i++){
        ints[i] = dt * ints[i];
    }

    /*
    // ints correcto
    printf("ints \n\n");
    for (int i = 0; i < nt; i++){
        printf( " %e ",ints[i]);
        if(i%4 == 0) printf("\n");
    }
    //*/


    // wts  = (ints>tol)*ints +(ints<=tol);
    double * wts = malloc(nt * sizeof(double *));

    // esta bien la primera parte con esto (ints>tol)*ints
    for (int i = 0; i < nt; i++){
        if( ints[i] > tol ){
            wts[i] = ints[i];
        }else{
            wts[i] = 0;
        }
    }

    for (int i = 0; i < nt; i++){
        if( wts[i] <= tol ){
            wts[i] = 1;
        }
    }

    /*
    printf("\n\n wts \n");
    for (int i = 0; i < nt; i++){
        printf( " %e ", wts[i] );
        if(i%4 == 0) printf("\n");
    }
    printf("\n\n fin wts \n");
    */


    // pout = dt*np.dot(pinput/np.transpose(wts), Pttau);

    // np.transpose(wts)
    double **wtsTrans = malloc(1 * sizeof(double *));
	for(int i = 0; i < nt; i++)
		wtsTrans[i] = malloc(nt * sizeof(double));
    // printf("\npase\n");
    
    for(int i = 0; i < nt; i++){
        wtsTrans[0][i] = wts[i];
    }
    // printf("\npase\n");

    // pinput/np.transpose(wts)
    for(int i = 0; i < nt; i++){
        wtsTrans[0][i] = pinput[i] / wtsTrans[0][i];
    }
    // printf("\npase\n");
    // de aca para arriba esta bien

    /*
    printf("\n\n wtsTrans \n");
    for (int i = 0; i < nt; i++){
        printf( " %e ", wtsTrans[0][i] );
        if(i%4 == 0) printf("\n");
    }
    printf("\n\n fin wtsTrans \n");
    //*/
    
    //np.dot(pinput/np.transpose(wts), Pttau);

    /*
    wtsTrans es una matriz de ntxnt, Pttau tambien.
    => resultado deberia ser una matrz de ntxnt pero da otra cosa...

    */
    double ** pout = multiplyMatrices(wtsTrans, Pttau, 1, nt, nt, nt);
    // printf("\npase\n");

    escalarMatrixMultiplication(pout, 1, nt, dt);
    /*
    printf("\n\n FINALLL \n\n");
    for (int i = 0; i < 1; i++){
        for (int j = 0; j < nt; j++){
            printf( "%e  ", pout[i][j] );
            if(j%4 == 0) printf("\n");
        }
        //if(i%10 == 0) break;
    }//*/

    // printf("\n\n FINAAL \n");


    

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

    /*    
    for (int i = 0; i < nt; i++){
        for (int j = 0; j < nt; j++){
            printf( " %f ", TTau[i][j] );
        }
        printf("\n");
    } 
    */
}


double ** multiplyMatrices(double ** firstMatrix, double ** secondMatrix, int rowFirst, int columnFirst, int rowSecond, int columnSecond)
{
	// Initializing elements of matrix mult to 0.
    double **matrix = malloc(1 * sizeof(double *));
	for(int i = 0; i < columnFirst; i++)
		matrix[i] = malloc(columnFirst * sizeof(double));

    for(int i = 0; i < columnFirst; i++)	
        matrix[0][i] = 0;
    
	// Multiplying matrix firstMatrix and secondMatrix and storing in array mult.
	for(int i = 0; i < rowFirst; i++)
	{
		for(int j = 0; j < columnSecond; j++)
		{   
			for(int k = 0; k < columnFirst; k++)
			{
				matrix[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
			}
		}
	}

    return matrix;
}


void main(){
    gettimeofday(&start,NULL);
	int tmin = 30;
	int tmax = tmin + 100; //130
	int nd = tmax - tmin; // 100
	int td = 4;
	int nt = nd * td; // 400
    printf("%i",nt);
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
    for (int i=0; i<nt; ++i){ 
            //printf("\n soy hilo: %3d (vuelta:%3d) -> Calculado:  d[i]=a + (i*delta) =  %lf  \n", omp_get_thread_num(), i, d[i]);
             //printf("1 -> %lf\n",tmps[i]);
    }    
    
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
    // printf("\n\n");
    vonFoerster(dt,t,tau,nt,tmps,hmrs,pnu,pdes,pinput);    
    gettimeofday(&end,NULL);
    double delta = ((end.tv_sec - start.tv_sec)*1000000u+end.tv_usec-start.tv_usec)/1.e6;
    printf("\n\nDemora: %f\n", delta);
    
}


void linspace (double a, double b, int c,double *d){ // funciona bien
    //printf("%f - %f - %i\n",a,b,c);    
    double aa = (double)b-a;
    //printf("%f\n\n",aa);
    double bb = (double)c-1;
    //printf("%f\n\n",bb);    
    double delta =aa/bb;
    //printf("%f\n\n",delta);

    int i;

     #pragma omp parallel shared(a, d, c, delta, chunk) private(i)
    {
        #pragma omp for schedule(static, chunk) nowait
            for (i=0; i<c; i++){            
                    d[i]=a + (i*delta);
            }    
    }

}

/*
Devuelve las coordenadas 2-D de la cuadrícula basadas en las coordenadas contenidas en los vectores x, y. 
X es una matriz donde cada fila es una copia de x, y Y es una matriz donde cada columna es una copia de y. 
La cuadrícula representada por las coordenadas X e Y tiene length(y)filas y length(x)columnas.
*/
void meshgrid(double *x, double *y, double **X, double **Y, int nt){

    int i, j;
    #pragma omp parallel shared(x,X,y,Y, nt) private(i, j)
    {
        //tid = omp_get_thread_num();
        // sections -> vamos a crear secciones paralelas
        // no wait -> cuando termina el hilo no espera al otro
        // cada seccion es para un hilo diferente
        // openmp crea un hilo para cada seccion definida
        #pragma omp sections nowait
        {

        #pragma omp section // primera seccion
        
                for(int i=0; i<nt;i++){ // fila i
                    for(int j=0;j<nt;j++){ // columna j
                        X[i][j] = x[j]; 
                    }
                }

            
        #pragma omp section 

                for(int i=0; i<nt;i++){ // columna i
                    for(int j=0;j<nt;j++){ // fila j
                        Y[j][i] = y[j];
                    }
                }

        }

  }
    // for(int i=0; i<nt;i++){ // fila i
    //     for(int j=0;j<nt;j++){ // columna j
    //         X[i][j] = x[j]; 
    //     }
    // }

    // for(int i=0; i<nt;i++){ // columna i
    //     for(int j=0;j<nt;j++){ // fila j
    //         Y[j][i] = y[j];
    //     }
    // }
 

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
    int i, j;
    double ** res = malloc(filas * sizeof(double *));
	for(int i = 0; i < filas; i++)
		res[i] = malloc(columnas * sizeof(double));

    #pragma omp parallel shared(filas, columnas, res, chunk) private(i,j)
    {
        #pragma omp for schedule(guided, chunk) nowait
            for (i = 0; i < filas; i++){
                for (j = 0; j < columnas; j++){
                    res[i][j] = 1;
                }
            }
    }

    return res;

    // double ** res = malloc(filas * sizeof(double *));
	// for(int i = 0; i < filas; i++)
	// 	res[i] = malloc(columnas * sizeof(double));

    // for (int i = 0; i < filas; i++){
    //     for (int j = 0; j < columnas; j++){
    //         res[i][j] = 1;
    //         //printf(" %f ,", res[i][j]);
    //     }
    //     //printf("\n");
    // }
    // return res;
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
    int i;
    #pragma omp parallel shared(tmeds, t, T0, T1, chunk, nt) private(i)
    {
        #pragma omp for schedule(guided, chunk) nowait

            for(i=0; i<nt;i++){
                tmeds[i] = T0+T1*cos(2*PI/365*t[i]);
            }
    }
    

    // printf("\n\n tmps");

    #pragma omp parallel shared(tmeds, t, tmps, T365,  chunk, nt) private(i)
    {
            #pragma omp for schedule(guided, chunk) nowait      
                for(int i=0; i<nt; i++){
                    tmps[i] = tmeds[i] - T365*cos(2*PI*(t[i]-(int)t[i]));    
                }        

    }
    // printf("\n\n");   


    //  for(int i=0; i<nt;i++){
    //     if(i%5) 
    //         printf("%f - ",tmps[i]);
    //             else
    //         printf("\n%f - ",tmps[i]);
    // }

    // T1 /= 2; // 7
    // T365 /=2; // 7 
    // for(int i=0; i<nt;i++){
    //     tmeds[i] = T0+T1*cos(2*PI/365*t[i]);
    //     /*
    //     if(i%5)
    //         printf("%f - ",tmeds[i]);
    //     else
    //         printf("\n%f - ",tmeds[i]);
    //         */
    // }
    // printf("\n\n");      
    //  for(int i=0; i<nt;i++){
    //     tmps[i] = tmeds[i] - T365*cos(2*PI*(t[i]-(int)t[i]));
    //     /*
    //     if(i%5)
    //         printf("%f - ",tmps[i]);
    //     else
    //         printf("\n%f - ",tmps[i]);
    //         */
    // }        

}

/*
    para vector: suma cada uno de los valores de las casillas del vector que es enviado como parametro 
    junto con el tamaño del mismo y retorna el valor de dicha suma
*/
double sum(double *m, int size){

    double ret = 0;
    int i;
    #pragma omp parallel for default(shared) private(i) schedule(guided,chunk) reduction(+:ret)

        for(i=0; i<size;i++){
            ret += (int) m[i];        
        }
    return ret;

    //  double ret = 0;
    // for(int i=0; i<size;i++){
    //     ret += (int) m[i];        
    // }
    // return ret;
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
    int i;
    double **res = malloc(n * sizeof(double *));
	for(i = 0; i < n; i++)
		res[i] = malloc(n * sizeof(double));
    

    #pragma omp parallel shared(v, u, res, n, chunk) private(i)
    {
        #pragma omp for schedule(guided, chunk) nowait

            for (i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    res[j][i] = v[i][0]*u[i];
                }
            }
    }

    
    return res;

    // double **res = malloc(n * sizeof(double *));
	// for(int i = 0; i < n; i++)
	// 	res[i] = malloc(n * sizeof(double));

    // for (int i = 0; i < n; i++){
    //     for (int j = 0; j < n; j++){
    //         res[j][i] = v[i][0]*u[i];
    //         //printf(" %f ,", res[i][j]);
    //     }
    //     //printf("\n");
    //     //if(i > 5) break;
    // }

    // return res;
}


double * cumtrapz(double *hs, int sz)
{    
    double *ct = malloc(sz * sizeof(double *));
    for(int i = 0; i < sz; i++) ct[i] = 0;    
    int i;
    // shcdule -> Esta variable de entorno establece el tipo de programación y el tamaño de segmento para todos estos bucles. 
    // El tamaño del trozo puede proporcionarse como un número entero, siendo el valor predeterminado 1.
    // static -> las iteraciones se dividen en trozos de tamaño y se asignan estáticamente a cada uno de los hilos de una manera de round robin
    #pragma omp parallel shared(hs,sz, ct, chunk) private(i)
    {
        #pragma omp parallel for schedule(guided,chunk)
         for(i = 1; i < sz; i++)
            ct[i] = ct[i - 1] + ( hs[i - 1] + hs[i] ) / 2.0; 
    }
    // for(int i = 0; i < sz; i++) ct[i] = 0;

    // for(int i = 1; i < sz; i++)
    //     ct[i] = ct[i - 1] + ( hs[i - 1] + hs[i] ) / 2.0; 
    return ct;
}

double ** transposef(double **m, int r, int c)
{
    int i, j;
    double **transpose = malloc(c * sizeof(double *));
	for(i = 0; i < c; i++)
		transpose[i] = malloc(r * sizeof(double));
   
    #pragma omp parallel shared(m, r, c, transpose, chunk) private(j, i)
    {
        #pragma omp for schedule(guided, chunk) nowait

            for(i=0; i < r; i++)
                for( j=0; j < c; j++)
                    transpose[j][i] = m[i][j];
     }           

    return transpose;
}

double * trapz(double **v, double filas, double columnas)
{
    int f = filas;
    int c = columnas;
    double mat1[f][c-1];
	double mat2[f][c-1];    
	// for(int i = 0; i<filas;i++){
	// 	for(int j = 0; j<(columnas-1);j++){
	// 		mat1[i][j] = v[i][j+1];
	// 	}
	// }

	// for(int i = 0; i<filas;i++){
	// 	for(int j = 0; j<(columnas-1);j++){
	// 		mat2[i][j] = v[i][j];
	// 	}
	// }

    int i, j;
    #pragma omp parallel shared(mat1, mat2, v, f, c) private(i, j)
    {    
        #pragma omp sections nowait
        {

        #pragma omp section // primera seccion
            for( i = 0; i<f;i++){
                for( j = 0; j<(c-1);j++){
                    mat1[i][j] = v[i][j+1];
                }
            }
                            
        #pragma omp section 
            for( i = 0; i<f;i++){
                for( j = 0; j<(c-1);j++){
                    mat2[i][j] = v[i][j];
                }
            }           
        }

  }
    /*
	// ----------------------- PRINT DE LO OBTENIDO -----------------------------
	for(int i = 0; i<filas;i++){
		for(int j=0;j<(columnas-1);j++){
			if(!j%(int)columnas)
				printf("\n%lf - ",mat1[i][j]);
			else
			printf("%lf - ",mat1[i][j]);
		}
	}
	printf("\n");
	for(int i = 0; i<filas;i++){
		for(int j=0;j<(columnas-1);j++){
			if(!j%(int)columnas)
				printf("\n%lf - ",mat2[i][j]);
			else
			printf("%lf - ",mat2[i][j]);
		}
	}
	// ----------------------------------------------------------------------------------
*/

        // double mat3[(int)filas][(int)columnas-1];
        // for(int i = 0; i<filas;i++){
        //     for(int j=0;j<(columnas-1);j++){
        //         mat3[i][j] = (mat1[i][j] + mat2[i][j])/2;
        //     }
        // }

double mat3[f][c];
#pragma omp parallel shared(mat1, mat2, mat3,f,c, chunk) private(i,j)
    {
        #pragma omp for schedule(guided, chunk) nowait        
        for(i = 0; i<f;i++){
            for(j=0;j<(c-1);j++){
                mat3[i][j] = (mat1[i][j] + mat2[i][j])/2;
            }
        }
    }
/*
	// ----------------------- PRINT DE LO OBTENIDO -----------------------------
	printf("\n\n");
	for(int i = 0; i<filas;i++){
		for(int j=0;j<(columnas-1);j++){
			if(!j%(int)columnas)
				printf("\n%lf - ",mat3[i][j]);
			else
			printf("%lf - ",mat3[i][j]);
		}
	}
	// ----------------------------------------------------------------------------------*/
	// printf("\n\n");
	double ret[f];
	double sum;
	for(int i = 0; i<f;i++){
		sum = 0;
		for(int j=0;j<(c-1);j++){
			sum += mat3[i][j];
		}
		ret[i]=sum;
	}
/*
	for(int i = 0; i<filas;i++){
		printf("%lf - ",ret[i]);
	}
	printf("\n\n");
*/

    double * result = ret;
}


double ** diff(double ** a, double **b, int f, int c){
    double **arraydiff = malloc(f * sizeof(double *));
	for(int i = 0; i < f; i++)
		arraydiff[i] = malloc(c * sizeof(double));

    for (int i = 0; i < f; i++){
        for (int j = 0; j < c; j++){
           arraydiff[i][j] = a[i][j] - b[i][j];
        }
    }

    return arraydiff;
}


void absMatrix(double **m, int f, int c){
    for (int i = 0; i < f; i++){
        for (int j = 0; j < c; j++){
            if( m[i][j] < 0 ){
                m[i][j] = m[i][j] * (-1);
            }
        }
    }
}


void sumEscalarToMatrix(double **m, int f, int c, double val){
    for (int i = 0; i < f; i++){
        for (int j = 0; j < c; j++){
            m[i][j] = m[i][j] + val;
        }
    }
}


void escalarMatrixMultiplication(double **m, int f, int c, double val){
    int i,j;
    double aux;
    #pragma omp parallel for default(shared) private(i,j) schedule(guided,chunk) reduction(*:aux)
    for ( i = 0; i < f; i++){
        for (j = 0; j < c; j++){
            aux = m[i][j] * val;
            m[i][j] = aux;
            //m[i][j] = m[i][j] * val;
        }
    }
}