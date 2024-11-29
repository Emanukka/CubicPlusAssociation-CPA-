
#%%import

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import numpy as np



from StateVariables.State import *

from Methods.Volume import *





#%%4C

# ------------------------------- 4C scheme-----------------------------------



# Tc = [  647.140  ,   304.12   ]         # (K)
# a0      = [ 2.2519*1e5*1e-6    ,  3.5079 *1e5*1e-6    ]        # (bar.L²/mol²)
# b       = np.array([ 0.01556  ,  0.0272  ])*1e-3          # (L/mol)
# c1      = np.array([ 0.6108  ,   0.76  ])   # (-)
# epsilon = np.array([  141.371   ,    0.0    ])*1e5*1e-3          # (bar.L/mol)
# beta    = np.array([ 0.0608   ,    0.0    ])          # (-)


Tc = [  647.140  ,   304.12   ]         # (K)
a0      = [ 2.25189*1e5*1e-6    ,  3.50789 *1e5*1e-6    ]        # (bar.L²/mol²)
b       = np.array([ 0.015559  ,  0.02719  ])*1e-3          # (L/mol)
c1      = np.array([ 0.61079  ,   0.759  ])   # (-)
epsilon = np.array([  141.371   ,    0.0    ])*1e5*1e-3          # (bar.L/mol)
beta    = np.array([ 0.06079   ,    0.0    ])          # (-)

kij = {'k12':-0.1469}
# kij = -0.1469



# Parametros passados no dicionario devem estar no formato de array numpy , Tupla ou Lista 

sig = 1
eps= 0



ncomp =2 


water_co2= {

    # v notation for vector  ( va0  == vector of a0 parameters)

    # Cubic parameters

    'va0' :a0,
    'vb' :b,
    'vtc' :Tc,
    'kij': kij,
    'sig': sig,
    'eps': eps,
    'ncomp': ncomp,
    'vc1': c1,

    'scheme': ['4C','solvation'],

    'epsilon_cross': 14*1e3,
    'beta_cross': 0.00069,  

    'vepsilon': epsilon,
    'vbeta': beta,
    # criar eps cross e beta cross para misturar com pelo menos um composto que é solvatado ( co2 por exemplo )
     
}

cpa_eos = CPA(water_co2)


print(cpa_eos)
#%%2B
# ------------------------------- 2B scheme-----------------------------------



Tc = [  647.140  ,   304.12   ]         # (K)
a0      = [ 2.9813*1e5*1e-6    ,  3.5079 *1e5*1e-6    ]        # (bar.L²/mol²)
b       = np.array([ 0.0145  ,  0.0272  ])*1e-3          # (L/mol)
c1      = np.array([ 0.0001 ,   0.76  ])   # (-)
epsilon = np.array([  285.982   ,    0.0    ])*1e5*1e-3          # (bar.L/mol)
beta    = np.array([ 0.013   ,    0.0    ])          # (-)



kij = 0.157



# Parametros passados no dicionario devem estar no formato de array numpy , Tupla ou Lista 

sig = 1
eps= 0



ncomp =2 


water_co2_2B= {

    # v notation for vector  ( va0  == vector of a0 parameters)

    # Cubic parameters

    'va0' :a0,
    'vb' :b,
    'vtc' :Tc,
    'kij': kij,
    'sig': sig,
    'eps': eps,
    'ncomp': ncomp,
    'vc1': c1,

    'scheme': ['2B','solvation'],

    'epsilon_cross': 14*1e3,
    'beta_cross': 0.0103,  

    'vepsilon': epsilon,
    'vbeta': beta,
    # criar eps cross e beta cross para misturar com pelo menos um composto que é solvatado ( co2 por exemplo )
     
}

# cpa_eos = CPA(water_co2_2B)

#%%4C State

P = 1*1e5

T = 300

state_liq = State(cpa_eos,vx = np.array([0.5,.5]),Pressure= P , T=T,density_initialization='liquid')


#%%colab

p = 33210526.315789472

w= 0.9
# vx=np.array([w,1-w])
VZ=[0.00375788, 0.99624212]

vx=np.array(VZ)
T = 300
# # rx = 59196.36450609903

# # cpa_eos.Pressure(vx,T,rx)

# # y = Calc_Volume(cpa_eos,vx,T,p,state='liquid')

x= Calc_Volume(cpa_eos,vx,T,p,state='vapor')




# y = VolumeMeu(cpa_eos,vx,T,p,state='vapor')



# o volume do liquido bate com o do colab
# ja o vapor diverge: no colab, dá nan( com 2 iterações )
# aqui, converge pra um valor que efetivamente nao é a raiz, porque
# vi no grafico q n tem 2 raizes
# (se nao vi errado ne )

# ta dando nan (mas nao tava me avisando )
# entao o problema ta no calculo da pressao da cpa ()
# deve ta tend oerro por dividir por 0, ln(0), ou algo do tipo...

# iterações 10
# 4.5076560392776386e-05
# iterações 1
# 6.987220161479616e-06

# %%

# fn = np.arctan

# x0 = 100

# newton(fn,x0,f_obj=0,x_pos=0)

# %%
