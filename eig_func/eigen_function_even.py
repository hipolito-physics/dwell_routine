import time
import configparser
import sys
import matplotlib.pyplot as plt
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy.integrate import quad
from scipy import optimize
from dwell_analytic import *
config = configparser.ConfigParser()
config.read('config.ini')

# +-------------------------------+
# |         Parâmetros            |
# +-------------------------------+

m = float(config['parametros']['m'])
h_bar = float(config['parametros']['h_bar'])
V_B = float(config['potencial']['V_B'])
V_D = float(config['potencial']['V_D'])
n = float(config['delta']['n'])
a = float(config['poco']['a'])
L = float(config['poco']['L'])

# +--------------------------------+
# | FUNÇÃO DE ONDA PAR NORMALIZADA |
# +--------------------------------+

if teste_par != 0:
    def psi_n_p(x):
        if x >= -(L-a)/2 and x <= (L-a)/2:
            return psi_1_p(x)                     #Região I - IV
        if x >= (L-a)/2 and x <= (L+a)/2:
            return psi_2_p(x)                     #Região II
        if x >= (L+a)/2:
            return psi_3_p(x)                     #Região III
        if x >= -(L+a)/2 and x <= -(L-a)/2:
            return psi_2_p(x)                     #Região V
        if x <= -(L+a)/2:
            return psi_3_p(x)                     #Região VI
        
# +-------------------------------+
# | PLOT DAS FUNÇÕES DE ONDA PAR  |
# +-------------------------------+

    print('Função de onda Par para n_quantico = ',n_quantico)

    t=np.linspace(-7,7,1000)
    plt.plot(t, np.array(list(map(psi_n_p, t)))**2,"r-")      #Elevar ao quadrado
    plt.ylabel('Psi²(x)')
    plt.xlabel('x')
    plt.legend(['Função de onda Par'])
    plt.savefig('../graph/graph-01.png')
    plt.show()

else:
    
# +----------------------------------+
# | FUNÇÃO DE ONDA ÍMPAR NORMALIZADA |
# +----------------------------------+

    def psi_n_i(x):
        if x >= -(L-a)/2 and x <= (L-a)/2:
            return psi_1_i(x)                     #Região I - IV
        if x >= (L-a)/2 and x <= (L+a)/2:
            return psi_2_i(x)                     #Região II
        if x >= (L+a)/2:
            return psi_3_i(x)                     #Região III
        if x >= -(L+a)/2 and x <= -(L-a)/2:
            return psi_2_i(x)                     #Região V
        if x <= -(L+a)/2:
            return psi_3_i(x)                     #Região VI

# +--------------------------------+
# | PLOT DAS FUNÇÕES DE ONDA ÍMPAR |
# +--------------------------------+

    print('Função de onda Par para n_quantico = ',n_quantico)

    t=np.linspace(-7,7,1000)
    #plt.figure(figsize=(5,8))
    plt.plot(t, np.array(list(map(psi_n_i, t)))**2,"g-")      #Elevar ao quadrado
    plt.ylabel('Psi²(x)')
    plt.xlabel('x')
    plt.legend(['Função de onda Ímpar'])
    plt.savefig('../graph/graph-01.png')
    plt.show()



    
