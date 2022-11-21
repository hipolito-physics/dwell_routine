from dwell_analytic import *
import matplotlib.pyplot as plt

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

if not paridade == 0:
    print('Função de onda Par para n = ',n_quantico)
else:
    print('Função de onda Ímpar para n = ',n_quantico)

t=np.linspace(-2*L,2*L,1000)
plt.plot(t, np.array(list(map(psi_n, t)))**2,"b-")      #Elevar ao quadrado
plt.ylabel('Psi²(x)')
plt.xlabel('x')

if not paridade == 0:
    plt.legend(['Função de onda Par'])
else:
    plt.legend(['Função de onda Ímpar'])
plt.savefig('../dwell/Graph/graph-01.png')
plt.show()

#EOF



    
