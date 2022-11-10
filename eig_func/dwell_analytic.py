#Python 3.8.10 (default, Out 27 2022, 20:18:18)
#Paulo Victor Hipolito Pires
import time
import configparser
import sys
import matplotlib.pyplot as plt
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy.integrate import quad
from scipy import optimize
import warnings
warnings.filterwarnings('ignore')



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

print("......PARÂMETROS DA ROTINA...... \n")
print("m............................... %f\n" % m)
print("h_bar........................... %f\n" % h_bar)
print("V_B............................. %f\n" % V_B)
print("V_D............................. %f\n" % V_D)
print("a............................... %f\n" % a)
print("L............................... %f\n" % L)


# +-------------------------------+
# |        NÚMERO QUÂNTICO        |
# +-------------------------------+

n_quantico = int(input('Digite o número quântico desejado: '))

# +-------------------------------+
# |  CONDICIONAL NÚMERO QUÂNTICO  |
# +-------------------------------+

if n_quantico <= 0 or n_quantico > 10:
    print("Insira valores para o número quântico entre 1 e 10 e tente novamente")
    sys.exit()


# +-------------------------------+
# |       ENERGIAS PRÓPRIAS       |
# +-------------------------------+


#Energias Pares
E_p = np.array([np.nan,-47.270947583542906,np.nan,-40.482219511530566,np.nan,-34.54802237185967,np.nan,-23.6199750008219,np.nan,-8.593518615082614,np.nan])


#Energias Ímpares
E_i = np.array([np.nan,np.nan,-47.20237359739848,np.nan,-39.395511992912326,np.nan,-28.817307986268133,np.nan,-16.753385216577403,np.nan,-1.3049256077337787])


# +-------------------------------+
# |   EQUAÇÕES K_1 , K_2 , K_3    |
# +-------------------------------+

#Equações Pares
q_1_p = np.sqrt((2*m*(V_D-V_B+E_p)/h_bar**2)+0j)
q_2_p = np.sqrt((2*m*(V_D+E_p))/h_bar**2)
q_3_p = np.sqrt((2*m*E_p)/h_bar**2+0j)




#Equações Ímpares
q_1_i = np.sqrt((2*m*(V_D-V_B+E_i)/h_bar**2)+0j)
q_2_i = np.sqrt((2*m*(V_D+E_i))/h_bar**2)
q_3_i = np.sqrt((2*m*E_i)/h_bar**2+0j)


# +-------------------------------+
# |            SIGMAS             |
# +-------------------------------+

#Sigmas Pares
sigma_1_p = q_1_p*(L-a)/2
sigma_2_pos_p = q_2_p*(L+a)/2
sigma_2_neg_p = q_2_p*(L-a)/2
sigma_3_p = np.abs(q_3_p*(L+a)/2)
delta_p = sigma_2_pos_p - np.arctan(sigma_3_p/sigma_2_pos_p) + (2*n+1)*np.pi/2
f_1_p = sigma_2_neg_p*np.tan(sigma_2_neg_p-sigma_2_pos_p+np.arctan(sigma_3_p/sigma_2_pos_p))
f_2_p = sigma_1_p*np.tan(sigma_1_p)

#print(f_1_p-f_2_p)                #Teste verificando os zeros


#Sigmas Ímpares
sigma_1_i = q_1_i*(L-a)/2
sigma_2_pos_i = q_2_i*(L+a)/2
sigma_2_neg_i = q_2_i*(L-a)/2
sigma_3_i = np.abs(q_3_i*(L+a)/2)
delta_i = sigma_2_pos_i - np.arctan(sigma_3_i/sigma_2_pos_i) + (2*n+1)*np.pi/2
f_1_i = sigma_2_neg_i*np.tan(sigma_2_neg_i-sigma_2_pos_i+np.arctan(sigma_3_i/sigma_2_pos_i))
f_2_i = sigma_1_i*1/(np.tan(sigma_1_i))

#print(f_1_i+f_2_i)                #Teste verificando os zeros


# +-------------------------------+
# |   CONSTANTES DE NORMALIZAÇÃO  |
# +-------------------------------+

#Constantes Pares em função de H para a normalização
A_p = (np.sin(sigma_2_neg_p-delta_p)*np.cos(sigma_1_p)-q_2_p/q_1_p*np.cos(sigma_2_neg_p-delta_p)*np.sin(sigma_1_p))
F_p = -((q_2_p)*np.cos(sigma_2_pos_p-delta_p)-np.abs(q_3_p)*np.sin(sigma_2_pos_p-delta_p))/(2*np.abs(q_3_p)*np.exp(-sigma_3_p))

#Constantes Ímpares em função de H para a normalização
B_i = (np.sin(sigma_2_neg_i - delta_i)*np.sin(sigma_1_i)+(q_2_i/q_1_i)*np.cos(sigma_2_neg_i - delta_i)*np.cos(sigma_1_i))
F_i = -((q_2_i)*np.cos(sigma_2_pos_i-delta_i)-np.abs(q_3_i)*np.sin(sigma_2_pos_i-delta_i))/(2*np.abs(q_3_i)*np.exp(-sigma_3_i))


# +-------------------------------+
# |       TESTE DE CHECAGEM       |
# +-------------------------------+

#Checagem se o n_quantico pertence a classe das funções pares ou ímpares

teste_par = n_quantico % 2

if not teste_par == 0:

# +-------------------------------+
# |       NORMALIZAÇÃO - PAR      |
# +-------------------------------+

#Normalização para as autofunções Pares, encontrar valor de H
    
    def regiao_1_p(x):
        return (A_p*np.cos(q_1_p*x))[n_quantico]*np.conjugate((A_p*np.cos(q_1_p*x)))[n_quantico]

    def regiao_2_p(x):
        return (np.sin(q_2_p*x - delta_p))[n_quantico]*np.conjugate(np.sin(q_2_p*x - delta_p))[n_quantico]

    def regiao_3_p(x):
        return (F_p*np.exp(-np.abs(q_3_p)*x))[n_quantico]*np.conjugate(F_p*np.exp(-np.abs(q_3_p)*x))[n_quantico]

    int_1_p = quad(regiao_1_p,0,(L-a)/2)[0]
    #print(int_1_p)
    int_2_p = quad(regiao_2_p,(L-a)/2,(L+a)/2)[0]
    #print(int_2_p)

    int_3_p = quad(regiao_3_p,(L+a)/2,5)[0]
    #print(int_3_p)

    #Função para encontrar H
    def H_p(H):
        return H**2*(int_1_p + int_2_p + int_3_p) - 0.5

    #Solução para H negativo
    root_p = optimize.brentq(H_p,-10,0)
    #print(root_p)
    
# +-------------------------------+
# |  FUNÇÕES DE ONDA POR REGIÃO   |
# +-------------------------------+

#Funções de Onda Pares
    def psi_1_p(x):
        return (root_p * (A_p*np.cos(q_1_p*x))[n_quantico])

    def psi_2_p(x):
        if x > 0:
            return (root_p * np.sin(q_2_p*x - delta_p)[n_quantico])
        if x < 0:
            return (-root_p * np.sin(q_2_p*x + delta_p)[n_quantico])  

    def psi_3_p(x):
        if x > 0:
            return (root_p * (F_p*np.exp(-np.abs(q_3_p)*x)))[n_quantico]
        if x < 0:
            return (root_p * (F_p*np.exp(np.abs(q_3_p)*x)))[n_quantico]

else:
    
# +-------------------------------+
# |      NORMALIZAÇÃO - ÍMPAR     |
# +-------------------------------+

#Normalização para as autofunções Ímpares, encontrar valor de H

    def regiao_1_i(x):
        return (B_i*np.sin(q_1_i*x))[n_quantico]*np.conjugate((B_i*np.sin(q_1_i*x)))[n_quantico]

    def regiao_2_i(x):
        return (np.sin(q_2_i*x - delta_i))[n_quantico]*np.conjugate(np.sin(q_2_i*x - delta_i))[n_quantico]

    def regiao_3_i(x):
        return (F_i*np.exp(-np.abs(q_3_i)*x))[n_quantico]*np.conjugate(F_i*np.exp(-np.abs(q_3_i)*x))[n_quantico]

    int_1_i = quad(regiao_1_i,0,(L-a)/2)[0]
    #print(int_1_i)

    int_2_i = quad(regiao_2_i,(L-a)/2,(L+a)/2)[0]
    #print(int_2_i)

    int_3_i = quad(regiao_3_i,(L+a)/2,5)[0]
    #print(int_3_i)

#Função para encontrar H
    
    def H_i(H):
        return H**2*(int_1_i + int_2_i + int_3_i) - 0.5

#Solução para H negativo
    
    root_i = optimize.brentq(H_i,-10,0)
    #print(root_i)
    
    #Funções de Onda Ímpares
    
    def psi_1_i(x):
        return (root_i * (B_i*np.sin(q_1_i*x))[n_quantico])

    def psi_2_i(x):
        if x > 0:
            return (root_i * np.sin(q_2_i*x - delta_i)[n_quantico])
        if x < 0:
            return (root_i * np.sin(q_2_i*x + delta_i)[n_quantico])  

    def psi_3_i(x):
        if x > 0:
            return (root_i * (F_i*np.exp(-np.abs(q_3_i)*x)))[n_quantico]
        if x < 0:
            return (-root_i * (F_i*np.exp(np.abs(q_3_i)*x)))[n_quantico]



#EOF

