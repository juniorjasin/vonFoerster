# -*- coding: cp1252 -*-
import numpy as np
from math import * 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

   
def feval(funcName, *args):      # first argument must be a string
    return eval(funcName)(*args)

"""
zeros: devuelve una matriz/vector llena de ceros, de tamanño indicado
len: devuelve el numero de elementos de la lista
for i in range(1,len(hs)): va desde 1 a len(hs) - 1 =>  [ 1, len(hs) )

retorna una matriz/vector
"""
def cumtrapz(hs): # hs es una lista/vector
    ct = np.zeros(len(hs));
    for i in range(1,len(hs)):
        ct[i] = ct[i-1] + (hs[i-1] + hs[i]) / 2.0;
    return ct;
      
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

"""
recibe (vector,int,int,int)
return: tuplas de int
"""	
def  TempSim(t, T0, T1, T365):
#TEMPSIM simula tabla de Temperaturas para hemisferio sur
# t    : vector de tiempos en dias.fraccion desde t[0] instante inicial
# T0   : temperatura media
# T1   : amplitud termica anual
# T365 : amplitud termica diaria
    vcos  = np.vectorize(cos);
    print "vcos:"
    print vcos
    vint  = np.vectorize(int);
    T1    = T1/2;
    T365  = T365/2;
    tmeds = T0+T1*vcos(2*pi/365*t);
    tmps  = tmeds - T365*vcos(2*pi*(t - vint(t)));
    return (tmps, tmeds);
    
def Varianza(tmps, p):
#VARIANZA calcula la varianza de la distribuci�n del desarrollo 
#         en funci�n de la temperatura
    numin = 0.000223;
    nu = p[0] * tmps^2 + p[1] * tmps + p[2];
    #   nu = numin*ones(1,length(tmps));
    nu = (nu > numin) * nu + (nu <= numin) * numin;

"""
para metodo np.linespace() -> https://stackoverflow.com/questions/27028226/python-linspace-in-c

meshgrid: recibe dos arreglos y crea dos matrices. la primera con los valores del primer arreglo puestos en filas
y la segunda con los valores del segundo arreglo puestos en columnas uno al lado del otro.

dot: es el producto punto de dos arrays. 

ones: retorna una arreglo de la forma y tipo indicados, llena de unos.

transpose: Permute the dimensions of an array.(calcula la transpuesta)

trapz: integra a lo largo del eje dado, usando la regla trapezoidal compuesta.
"""     

def vonFoerster(dt, t, tau, nt, tmps, hmrs, pnu, fhnu, pdes, fhrates, pinput):
    """
    VONFOERSTER outputs evolution of an stage population based on vonFoster
    dt:      paso discreto de tiempo
	t:       vector de instantes de tiempos
	tau:     vector de instantes de tiempos
	nt:      n�mero total de espacios de tiempo
	tmps:    vector de temperaturas para todo instante discreto de tiempo
	hmrs:    vector de humedades relativas para todo instante discreto de tiempo
	pnu:     parametros de la funci�n varianza()
	fhnu:    handler de la funci�n varianza()
	pdes:    par�metros de la funci�n de desarrollo
	fhrates: handler de la funci�n de desarrollo
	pinput:  vector de densidades de poblaci�n correspondiente al estado 
    """
    tol     = 0.0001;                      # von Foerster model tolerance
    T,Tau   = np.meshgrid(t,tau);          # T,Tau matrices
    print "meshgrid"
    print "T"
    print T
    print "Tau"
    print Tau
    Pttau   = np.zeros((nt,nt));           # initialize matrix which will hold the convolution kernel
    wts     = np.zeros(nt);                # initialize vector which will hold weight for normalization
    rates   = np.zeros(nt);                # initialize vector of rates for each instant t
    hmrsnul = np.sum(abs(hmrs) > 0);
    
#    evalua con humedad y sin humedad
    if hmrsnul != 0:
        print "IFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
       # calcula las tasas de desarrollo para temperaturas y humedades relativas dadas 
        for i in range(nt):
            rates[i] = feval(fhrates, tmps[i], hmrs[i], pdes);       
    else:
        # 
        print "ELSEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
       # calcula las tasas de desarrollo  para temperaturas dadas
        for i in range(nt):
            rates[i] = feval(fhrates, tmps[i], pdes); 

    print "rates:\n"
    print rates
    # recibe: 1ro un matriz de nt filas y una columna, 2do. 
    print "ones * dt"
    print dt * np.ones((nt,1)) 
    print [cumtrapz(rates)]
    RT    = np.dot(dt * np.ones((nt,1)), [cumtrapz(rates)]);  # create a matrix which is the cumulative development
    print "///////////////////////////////////////n"
    print type(RT)
    print RT
    RTau  = np.transpose(RT); 
    print "transpose"
    print RTau                                   # create transpose of cumulative matrix for use in kernel
    vexp  = np.vectorize(exp);   
    vsqrt = np.vectorize(sqrt);    
    
    if (pnu[0] == 0) and (pnu[1] == 0):
        nu = pnu[2];
        print "nuuuuu"
        print nu
        Pttau  = (T>Tau)*vexp(-(1-(RT-RTau))**2/(4*nu*(abs(T-Tau)+tol)))\
                           /vsqrt(4*pi*nu*(abs(T-Tau)**3+tol));         # extended von foerster kernel
        print "IF PTTAUUUUUUUUUUUUUUUUUUUUUUUUUU"

        print "-------------------------------------------------------------\n\n abs"
        print abs(T-Tau)
        print "\n\n + tol"
        print (abs(T-Tau)+tol)

        print "FINAL ROUND \n\n\n"
        print (4*nu*(abs(T-Tau)+tol))

    else:
        #con estos valores no entra a calcular la varianza
        nus   = feval(fhnu, tmps, pnu);      # calcula las varianzas para temperatures dadas
        NU    = np.dot(np.ones(nt,1), nus);# crea una matriz de varianzas en funci�n de las temperaturas
        Pttau = (T>Tau)*exp(-(1-(RT-RTau))**2/(4*NU*(abs(T-Tau)+tol))) \
                           /vsqrt(4*pi*NU*(abs(T-Tau)**3+tol));         # extended von foerster kernel
        print "ELSE PTTAUUUUUUUUUUUUUUUUUUUUUUUUUU"
    
    ints = dt*np.trapz(Pttau,axis=1);    # integrate in columns to normalize
    wts  = (ints>tol)*ints+(ints<=tol);  # calculate a  weighting factor; 
                                         # make it one if the integral is too  small (< tol)
    pout = dt*np.dot(pinput/np.transpose(wts), Pttau);      # output distribution, normalized by integral of P

# Graficaci�n 3D de la funci�n de probabilidad Pttau en funcion de t y de t*rates
    global idCorrida;
    global codfig;
    if codfig > 7:
        fh = plt.figure(codfig);   
        ax = fh.gca(projection='3d')
        plt.hold (True);
        plt.grid (True, which='both');
        plt.xlabel('Dias');
        plt.ylabel('Tasas de desarrollo acumuladas');
        #plt.zlabel('Probabilidad');
        plt.title('Probabilidad de emergencia segun von Foerster Extendido');
#        [X,Y] = meshgrid(1:nt, 1:nt);
        surf = ax.plot_surface(T, RT, Pttau, rstride=1, cstride=1, \
        cmap=cm.coolwarm, linewidth=0, antialiased=False);
#       ax.set_zlim(-1.01, 1.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fh.colorbar(surf, shrink=0.5, aspect=5)
        plt.show();
#        fig = ['F',np.num2str(codfig), idCorrida, '.eps'];
#        print ('-dpsc2',fig);
#        hold(False); 
#        fh.close;
#        codfig = codfig + 1; 
	
    #Functional return of vonFoerster
    return (T, Tau, rates, RT, Pttau, ints, wts, pout) # Ouput population distribution in time
 
#Generaci�n de datos de prueba de una corrida 
tmin = 30;         # dia inicial de la corrida 30 de enero
tmax = 30 + 100;   # dia final de la corrida 30 de enero mas 100 dias
nd = tmax - tmin;  # numero de dias
td = 4;            # intervalos de tiempo diarios
nt = nd * td;      # numero total de pasos de tiempo
dt = nd/float(nt);        # paso discreto de tiempo
print "dt:"
print dt
t    = np.linspace(tmin,tmax-dt,nt); # vector de tiempos, t
tau  = t;                            # vector of tiempos, tau
T0   = 15;
T1   = 15;
T365 = 15; 
tmps, tmeds = TempSim(t, T0, T1, T365);
print "return TempSim:"
print tmps
print "---------------------------------------"
print tmeds
hmrs  = np.zeros(nt);
# Datos para huevos de Diatrea 
fhrates = "BriereI";
pdes    = np.zeros(3); 
pdes[0] =  0.000131734;
pdes[1] = 10.25740308;  
pdes[2] = 36.65400490;
# Descripci�n de la funci�n de varianza en funci�n de las temperaturas
fhnu    = "Varianza";
pnu     = np.zeros(3); 
pnu[0]  = 0.0;
pnu[1]  = 0.0;  
pnu[2]  = 0.000223;
# Datos de la poblacion de ingreso Pulso 100 individuos 
pinput  =  np.zeros(nt);
for i in range(4):
   pinput[i] = 25;
idCorrida = "vFPy";
codfig = 8;
T, Tau, rates, RT, Pttau, ints, wts, pout = vonFoerster(dt, t, tau, nt, tmps, hmrs, pnu, fhnu, pdes, fhrates, pinput); 
#plt.plot(t, pinput, "b", t, pout, "r");
