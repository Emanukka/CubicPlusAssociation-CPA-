
# Utilização dos parametros
#%%
import numpy as np
import matplotlib.pyplot as plt
# import jax.numpy as np
import timeit
# import numpy as np
import sys
import os



# sys.path.append(os.path.abspath((__file__),))

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))



from StateVariables.State import *


sig = 0
eps= 1

#------------------------------- 2B scheme-----------------------------------
    #                 1-propanol  n-hexane



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

    'scheme': ['2B','inert'],
    'vepsilon': epsilon_ph,
    'vbeta': beta_ph,
     
}

# Parametros

# talvez seja desnecessário ter classe para os parametros ( nao acrescenta em nada essa 
# abstração  ) , ja que os parametros serao escolhidos facilmente via palavras chave 

# Object Cpa is initialized only by passing parameters from 2 objects (CEosParamns and AssocEosParamans) 

cpa_eos_prop_hex= CPA(propanol_hexane)

# # #

# Debug 

x_0 = np.array([.5,.5])
T_0 = 296.45
P_0 = 1e5

# print(cpa_eos)
# state_liq = State( eos=cpa_eos,vx = x_0,T = T_0,Pressure=P_0,density_initialization="liquid")
# state_vap = State( eos=cpa_eos,vx = x_0,T = T_0,Pressure=P_0,density_initialization="vapor")

# print(state_liq.phi())

# print(state_liq.rho)
# print(state_vap.rho)

# +----------------+-------------+-----------+-----------+-----------+
# | Liquid Phase   | Component   |       Phi |    X_prod |   Density |
# |----------------+-------------+-----------+-----------+-----------|
# | -              | Comp1       | 0.0364874 | 0.0382548 |   9733.04 |
# | -              | Comp2       | 0.320284  | 0.0382548 |   9733.04 |
# +----------------+-------------+-----------+-----------+-----------+
# +---------------+-------------+----------+----------+-----------+
# | Vapor Phase   | Component   |      Phi |   X_prod |   Density |
# |---------------+-------------+----------+----------+-----------|
# | -             | Comp1       | 0.900255 | 0.900551 |   43.4594 |
# | -             | Comp2       | 0.97432  | 0.900551 |   43.4594 |
# +---------------+-------------+----------+----------+-----------+

# do colab

# ----- Liquido -----
# X_Liq_inteira =[[0.19558838 1.        ]
#  [0.19558838 1.        ]]
# V_liq = 0.00010274281994591246
# X_Liq =0.03825481435064395
# Phi_Liq =[0.03648741 0.32028413]
# ----- Vapor -----
# V_vap = 0.023009975171692513
# X_vap = 0.9005506077781109
# Phi_Vap = [0.90025524 0.9743199 ]

#%%plot

T = 23.3 + 273.15 # (K)
P = 1.0 *1e5         # (bar)

#                     x1      X1
expdata = np.array([[0.008, 0.730],
                    [0.016, 0.458],
                    [0.030, 0.341],
                    [0.043, 0.277],
                    [0.053, 0.244],
                    [0.070, 0.183],
                    [0.084, 0.149]])

# from N. von Solms, L. Jensen et al., "Measurement and modeling of hydrogen bonding in 1-alkanol + n-alkane binary mixtures", Fluid Phase Equilibria, 261, 1–2, 272 (2007)

#%%plot

# do igor

NP = 101
eps = 1E-20

x1 = np.linspace(0.0 + eps, 0.1, NP)
x2 = 1 - x1


T_1 = 23.3 + 273.15 # (K)
P_1 = 1.0 *1e5         # (bar)

X2B = np.zeros_like(x1)



for n in range(NP):

    x = np.array([x1[n], x2[n]])

    # rho2B = cpa2B.VARBLdensity(R, T, P, x, 'LIQ')[0]


    # liquid phase monomer fraction (X.prod)

    rho_liq = 1/Calc_Volume(cpa_eos_prop_hex,x,T_1,P_1,state='liquid')['x']


    # monomerfrac2B, m2B, error2B, it2B = cpa2B.VARBLsitemonofrac(R, T, rho2B, x)


    x_prod = cpa_eos_prop_hex.X_product(x,T_1,rho_liq)
    # np.prod()

    X2B[n] = x_prod

    # x_val = X_assoc(rho_liq,vx,S,ncomp:int,nsites:int,delta)



    # monomerfrac3B, m3B, error3B, it3B = cpa3B.VARBLsitemonofrac(R, T, rho3B, x)

    # X2B[n] = monomerfrac2B.prod()
    # X3B[n] = monomerfrac3B.prod()

# X2B[0] = np.nan
# X3B[0] = np.nan

#%%

plt.scatter(expdata[:,0], expdata[:,1], color = 'none', edgecolor = 'black', marker = 's', label = 'von Solms {\it et al.} (2007)', zorder = +3)
plt.plot(x1, X2B, color = 'black', linewidth = 1.0, label = 'CPA-2B', ls='dashed')
plt.title(f'Monomer fraction of 1-propanol + n-hexane at $23.3 ^o$C', fontsize = 8)
plt.ylabel('monomer fraction of propanol [ - ]', fontsize = 8)
plt.xlabel('mole fraction of propanol [ - ]', fontsize = 8)
plt.legend(loc = 'best', fontsize = 8)
plt.ylim((0.00, 1.00))
plt.xlim((0.00, 0.10))
plt.tight_layout()
# 


# plt.savefig(fname ='propanol_hexane')

# plt.savefig(r'D:\CPA\CPA_mod\Exemplos\Plots\propanol_hexane_Xproduct.png')

plt.show()

