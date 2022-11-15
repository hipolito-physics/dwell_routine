#Python 3.8.10 (default, Out 27 2022, 20:18:18)
import time
import configparser
import sys
import numpy as np
from scipy.integrate import quad
from scipy import optimize
import warnings 

warnings.simplefilter("ignore", np.ComplexWarning) #A função quad dá warnings quando o número é real porem do tipo n + 0j, onde n é real.


#Leitura do arquivo de configuração
config = configparser.ConfigParser()
config.read('config.ini')

# Início da contagem do tempo
inicio_tempo = time.time() #mudar para inicio da rotina


# +-------------------------------+
# |         Parâmetros            |
# +-------------------------------+



mass = float(config['parametros']['m'])
h_bar = float(config['parametros']['h_bar'])
V_B = float(config['potencial']['V_B'])
V_D = float(config['potencial']['V_D'])
n = float(config['delta']['n'])
a = float(config['poco']['a'])
L = float(config['poco']['L'])

if not __name__ == '__main__':
    print("......PARÂMETROS DA ROTINA...... \n")
    print("mass............................... %f\n" % mass)
    print("h_bar........................... %f\n" % h_bar)
    print("V_B............................. %f\n" % V_B)
    print("V_D............................. %f\n" % V_D)
    print("a............................... %f\n" % a)
    print("L............................... %f\n" % L)



# +-------------------------------+
# |       ENERGIAS PRÓPRIAS       |
# +-------------------------------+


#Energias Pares
E_p = np.loadtxt('../dwell/dat/eig_energy_p.txt')
#print('Energias Pares \n',E_p)

#Energias Ímpares
E_i = np.loadtxt('../dwell/dat/eig_energy_i.txt')
#print('\n Energias Ímpares \n',E_i)

#Energias Totais
E_t = np.sort(np.concatenate((E_p, E_i)))
#print(E_t)

# +-------------------------------+
# |        NÚMERO QUÂNTICO        |
# +-------------------------------+

n_quantico = int(config['n_quantico']['n_quantico'])


# +-------------------------------+
# |  CONDICIONAL NÚMERO QUÂNTICO  |
# +-------------------------------+

if n_quantico <= 0 or n_quantico > len(E_t):
    print("Insira valores para n_quantico entre 1 e ",len(E_t)," e tente novamente")
    sys.exit()

# +-------------------------------+
# |   EQUAÇÕES Q_1 , Q_2 , Q_3    |
# +-------------------------------+

q_1 = np.sqrt((2*mass*(V_D-V_B+E_t)/h_bar**2)+0j)[n_quantico - 1]
q_2 = np.sqrt((2*mass*(V_D+E_t))/h_bar**2)[n_quantico - 1]
q_3 = np.sqrt((2*mass*E_t)/h_bar**2+0j)[n_quantico - 1]




# +-------------------------------+
# |            SIGMAS             |
# +-------------------------------+

#Sigmas Pares
sigma_1 = q_1*(L-a)/2
sigma_2_pos = q_2*(L+a)/2
sigma_2_neg = q_2*(L-a)/2
sigma_3 = np.abs(q_3*(L+a)/2)
delta = sigma_2_pos - np.arctan(sigma_3/sigma_2_pos) + (2*n+1)*np.pi/2
f_1 = sigma_2_neg*np.tan(sigma_2_neg-sigma_2_pos+np.arctan(sigma_3/sigma_2_pos))

paridade = n_quantico % 2 #Variável para teste se n_quantico é par ou ímpar.

if not paridade == 0:
    f_2 = sigma_1*np.tan(sigma_1)
    #print(f_1 - f_2)                   #Checagem dos zeros das Energias pares

else:
    f_2 = sigma_1*1/(np.tan(sigma_1))
    #print(f_1 + f_2)                   #Checagem dos zeros das Energias ímpares


# +------------------------------------+
# |   CONSTANTES DE NORMALIZAÇÃO PARES |
# +------------------------------------+

F = -((q_2)*np.cos(sigma_2_pos-delta)-np.abs(q_3)*np.sin(sigma_2_pos-delta))/(2*np.abs(q_3)*np.exp(-sigma_3))

#Condicional para n_quantico Par
if not paridade == 0:
    A = (np.sin(sigma_2_neg-delta)*np.cos(sigma_1)-q_2/q_1*np.cos(sigma_2_neg-delta)*np.sin(sigma_1))

#Condicional para n_quantico Ímpar    
else:
    B = (np.sin(sigma_2_neg - delta)*np.sin(sigma_1)+(q_2/q_1)*np.cos(sigma_2_neg - delta)*np.cos(sigma_1))
    
    
# +-------------------------------+
# |         NORMALIZAÇÃO          |
# +-------------------------------+

def regiao_1(x):
    if not paridade == 0:
        return (A*np.cos(q_1*x))*np.conjugate((A*np.cos(q_1*x))) #Região 1 Par
    else:
        return B*np.sin(q_1*x)*np.conjugate(B*np.sin(q_1*x)) #Região 1 Ímpar

def regiao_2(x):
    if not paridade == 0:
        return np.sin(q_2*x - delta)*np.conjugate(np.sin(q_2*x - delta)) #Região 2 Par
    else:
        return np.sin(q_2*x - delta)*np.conjugate(np.sin(q_2*x - delta)) #Região 2 Ímpar

def regiao_3(x):
    if not paridade == 0:
        return F*np.exp(-np.abs(q_3)*x)*np.conjugate(F*np.exp(-np.abs(q_3)*x)) #Região 3 Par
    else:
        return F*np.exp(-np.abs(q_3)*x)*np.conjugate(F*np.exp(-np.abs(q_3)*x)) #Região 3 Ímpar


int_1 = quad(regiao_1,0,(L-a)/2)[0]
#print(int_1)
int_2 = quad(regiao_2,(L-a)/2,(L+a)/2)[0]
#print(int_2)

int_3 = quad(regiao_3,(L+a)/2,5)[0]
#print(int_3)
    
#Função para encontrar H
def normal_constant(H):
    return H**2*(int_1 + int_2 + int_3) - 0.5

#Solução para H negativo
root = optimize.brentq(normal_constant,-10,0)
#print(root)

# +-------------------------------+
# |  FUNÇÕES DE ONDA POR REGIÃO   |
# +-------------------------------+

#Região 1
def psi_1(x):
    
    if not paridade == 0:
        return root*A*np.cos(q_1*x) #Função de onda 1 Par
    else:
        return root*B*np.sin(q_1*x) #Função de onda 1 Ímpar

#Região 2   
def psi_2(x):
    
    func_2 = root*np.sin(q_2*x - delta) 
    func_2_neg = root*np.sin(q_2*x + delta)
    
    if not paridade == 0:
        if x > 0:
            return func_2 #Função de onda 2 Par positiva
        if x < 0:
            return -func_2_neg #Função de onda 2 Par negativa
    else:
        if x > 0:
            return func_2 #Função de onda 2 Ímpar positiva
        if x < 0:
            return func_2_neg #Função de onda 2 Ímpar negativa

#Região 3
def psi_3(x):
    func_3 = root*F*np.exp(-np.abs(q_3)*x)
    func_3_neg = root*F*np.exp(np.abs(q_3)*x)
    
    if not paridade == 0:
        if x > 0:
            return func_3 #Função de onda 3 Par positiva
        if x < 0:
            return func_3_neg #Função de onda 3 Par negativa
    else:
        if x > 0:
            return func_3 #Função de onda 3 Ímpar positiva
        if x < 0:
            return -func_3_neg #Função de onda 3 Ímpar negativa

# +-------------------------------+
# |    Teste de funcionamento     |
# +-------------------------------+

if __name__ == '__main__':
    print('A rotina está funcionando. Execute o arquivo eigen_function.py.')
    
fim_tempo = time.time()
variacao_tempo = fim_tempo - inicio_tempo

if variacao_tempo <= 60:
    print("Tempo de execução: %f" % (float(variacao_tempo)),"segundos")


#EOF

