# Rotina desenvolvida por Jiusandro Kuhn
#com ajustes implementados
from scipy.optimize import brentq
import numpy as np



def zbrak(func,xmin,xmax,nsteps):

    dx = abs(xmax-xmin)/nsteps # Define a largura do passo
    ds = dx/4 # Caso a raiz esteja exatamente sobre o intervalo ele gera um pequeno deslocamnto ds


    a = xmin

    r=[] # gera um vetor nulo para armazenar os intervalos contendo as raízes


    # Teste de inicialização. Testa a existencia de uma raiz exatamente no início do intervalo
    if abs(func(a))<=1E-6:
        a -= ds # Desloca o início do intervalo para -ds
    
    n = 1 # Contador para os intervalos avaliados
    while n < nsteps:
        b = xmin+n*dx
        
        # Verifica se o fim do intervalo é uma raíz
        if abs(func(b))<=1E-6:
            b += ds # Desloca o fim do intervalo para +ds

        # Verifica se há mudança de sinal na função
        if func(a)*func(b) <=0 :
            r = r + [[a,b]] # Armazena o intervalo no vetor
        
        a = b
        n+=1

    return r



def find_roots(func,x1,x2,ndiv):
    # Monta uma lista contendo os pontos do lado direito e esquerdo das raízes. Essa lista serve de entrada para a rotina brentq.
    # x1 e x2 são o início e o fim do intervalo respectivamente. ndiv é o número de subdivisões do intervalo definido por x1 e x2.
    r = zbrak(func,x1,x2,ndiv)
    # print(r)

    nroots = len(r) # Conta o número de raízes existentes


    roots=[] # cria uma lista vazia para armazenar as raízes

    m = 0
    while m <= (nroots-1):
        xa = r[m][0] # Raíz do lado esquerdo
        xb = r[m][1] # Raiz do lado direito
        root = brentq(func,xa,xb) # Encontra a raíz

        #   Essa parte do código faz o teste se a raíz encontrada é de fato uma raíz da equação
        root_test = func(root)
        if abs(root_test)<1E-9:
             #print(root_test) # Teste para verificar se o valor encontrado é uma função, desative o comentário para funcionar
             roots.append([root])
        m+=1

    return(roots)

# +---------------------------------------+
# | TESTE DA ROTINA APENAS NESSE ARQUIVO  |
# +---------------------------------------+

if __name__ == '__main__':
    
    def func1(x):
        return x*np.sin(x)

    def func2(y):
        return y**2 + np.tan(y)

    zeros1 = print(find_roots(func1,0,10,1000))
    zeros2 = print(find_roots(func2,0,5,1000))
    
