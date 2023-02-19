from dwell_analytic import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

# +--------------------------------+
# | FUNÇÃO DE ONDA PAR NORMALIZADA |
# +--------------------------------+

def psi_n(x):
    if x >= -(L-a)/2 and x <= (L-a)/2:
        return psi_1(x)                     #Região I - IV
    if x >= (L-a)/2 and x <= (L+a)/2:
        return psi_2(x)                     #Região II
    if x >= (L+a)/2:
        return psi_3(x)                     #Região III
    if x >= -(L+a)/2 and x <= -(L-a)/2:
        return psi_2(x)                     #Região V
    if x <= -(L+a)/2:
        return psi_3(x)                     #Região VI
        
# +-------------------------------+
# | PLOT DAS FUNÇÕES DE ONDA PAR  |
# +-------------------------------+

#if not paridade == 0:
#    print('Função de onda Par para n = ',n_quantico)
#else:
#    print('Função de onda Ímpar para n = ',n_quantico)

t=np.linspace(-2*(L+a)/2,2*(L+a)/2, 2000)
figure(figsize=(10,10))
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 15
plt.plot(t, np.array(list(map(psi_n, t))),"k-")      #Elevar ao quadrado
plt.ylabel('$\psi(x)$', loc='center', fontsize=30)
plt.xlabel('$x$',  loc='right', fontsize=35)
#"$|\Psi_{n}(x)|^{2}$"
#if not paridade == 0:
#    plt.legend(['Função de onda Par'])
#else:
#    plt.legend(['Função de onda Ímpar'])


#if not paridade == 0:
#    plt.legend(['Função de onda Par'])
#else:
#    plt.legend(['Função de onda Ímpar'])
plt.savefig('../dwell/Graph/graph-01.png')
plt.show()

#EOF



    
