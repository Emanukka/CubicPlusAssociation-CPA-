import sys
import os

import numpy as np


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from StateVariables.State import *


sig = 0
eps= 1


#%% Propanol ,, hexano
# ------------------------------- 2B scheme-----------------------------------
#                     1-propanol  n-hexane



Tc_ph  = [  536.78  ,   507.6   ]         # (K)
a0_ph      = [ 11.9102*1e5*1e-6    ,  23.681 *1e5*1e-6    ]        # (bar.L²/mol²)

b_ph       = np.array([ 0.06411  ,  0.10789  ])*1e-3          # (L/mol)
c1_ph      = np.array([ 0.91709  ,   0.8313  ])   # (-)
epsilon_ph = np.array([  210.0   ,    0.0    ])*1e5*1e-3          # (bar.L/mol)
beta_ph    = np.array([ 0.0081   ,    0.0    ])          # (-)

kij_ph = 0.

# Parametros passados no dicionario devem estar no formato de array numpy , Tupla ou Lista 

sig = 1
eps= 0



ncomp =2 


propanol_hexane = {

    # v notation for vector  ( va0  == vector of a0 parameters)

    # Cubic parameters

    'va0' :a0_ph,
    'vb' :b_ph,
    'vtc' :Tc_ph,
    'kij': kij_ph,
    'sig': sig,
    'eps': eps,
    'ncomp': ncomp,
    'vc1': c1_ph,

    # Assoc parameters

    # assoc_scheme:

    # Ao estabelecer que ambos componentes são inertes , os termos da cpa se igualam ao da cubica 
    # (associação zera)

    # -> Isso ocorre porque a matriz de sitios terá compasta apenas por 0's, de tal maneira que 
    # a matriz dos sitios nao associativos será composta por 1's (nenhum componente associa)

    # Smatrix: [[0. 0.]
    #           [0. 0.]]


    'scheme': ['inert','inert'],

    'vepsilon': epsilon_ph,
    'vbeta': beta_ph,
     
}

# Parametros
x = np.array([.5,.5])

# talvez seja desnecessário ter classe para os parametros ( nao acrescenta em nada essa 
# abstração  ) , ja que os parametros serao escolhidos facilmente via palavras chave 

# Object Cpa is initialized only by passing parameters from 2 objects (CEosParamns and AssocEosParamans) 

cpa_eos = CPA(propanol_hexane)

# # #

# Debug 

x = np.array([.5,.5])
T = 296.45
P = 1e5

# print(cpa_eos)
state_liq = State( eos=cpa_eos,vx = x,T = T,Pressure=P,density_initialization="liquid")
state_vap = State( eos=cpa_eos,vx = x,T = T,Pressure=P,density_initialization="vapor")

print(state_liq)
# print(state_vap)


#   ------ Tabela para SRK apenas (sem parte associativa) ------ 

#  \\ Verificado com o igor : 

# ----- Liquido -----
# X_Liq_inteira =[[1. 1.]
#  [1. 1.]]
# Densidade = 9469.575072972342


# X_Liq_prod =1.0
# Phi_Liq =[0.69497321 0.19326048]

#  \\ Resultado daqui :

# +----------------+-------------+----------+----------+-----------+
# | Liquid Phase   | Component   |      Phi |   X_prod |   Density |
# |----------------+-------------+----------+----------+-----------|
# | -              | Comp1       | 0.694972 |        1 |   9469.58 |
# | -              | Comp2       | 0.19326  |        1 |   9469.58 |
# +----------------+-------------+----------+----------+-----------+
