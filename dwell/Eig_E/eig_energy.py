#Python 3.8.10 (default, Out 27 2022, 20:18:18)
import sys
sys.path.append("..")
from Eig_E import finder_roots
import configparser
import numpy as np

#

config = configparser.ConfigParser()
config.read('../config.ini')


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

# +-------------------------------+
# |       ENERGIAS PRÓPRIAS       |
# +-------------------------------+

# Energia Par:

def func_even(E):
    k_1 = np.sqrt((2*mass*(V_D - V_B + E)/h_bar**2)+0j)
    k_2 = np.sqrt((2*mass*(E+V_D))/h_bar**2)
    k_3 = np.sqrt((2*mass*E/h_bar**2)+0j)
    sigma_1 = k_1*(L-a)/2
    sigma_2_pos = k_2*(L+a)/2
    sigma_2_neg = k_2*(L-a)/2
    sigma_3 = np.absolute((k_3))*(L+a)/2
    f_1 = sigma_2_neg*np.tan(sigma_2_neg - sigma_2_pos + np.arctan(sigma_3/sigma_2_pos))
    f_2 = np.real((sigma_1*np.tan(sigma_1)))
    return f_1 - f_2

# Salva as Energias respectivas

E_p = [finder_roots.find_roots(func_even,-V_D+0.01,-0.01,1000)]
np.savetxt("../dat/eig_energy_p.txt",E_p[0])
if __name__ == '__main__':
    print('\nEnergias Pares: \n',E_p)

# Energia Ímpar:

def func_odd(E):
    k_1 = np.sqrt((2*mass*(V_D - V_B + E)/h_bar**2)+0j)
    k_2 = np.sqrt((2*mass*(E+V_D))/h_bar**2)
    k_3 = np.sqrt((2*mass*E/h_bar**2)+0j)
    sigma_1 = k_1*(L-a)/2
    sigma_2_pos = k_2*(L+a)/2
    sigma_2_neg = k_2*(L-a)/2
    sigma_3 = np.absolute((k_3))*(L+a)/2
    f_1 = sigma_2_neg*np.tan(sigma_2_neg - sigma_2_pos + np.arctan(sigma_3/sigma_2_pos))
    f_2 = np.real((sigma_1*1/(np.tan(sigma_1))))
    return f_1 + f_2

# Salva as Energias respectivas

E_i = [finder_roots.find_roots(func_odd,-V_D+0.01,-0.01,1000)]
np.savetxt("../dat/eig_energy_i.txt",E_i[0])
if __name__ == '__main__':
    print('\n Energias Ímpares: \n',E_i)
      

#EOF

