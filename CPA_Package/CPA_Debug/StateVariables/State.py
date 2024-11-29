import numpy as np
import sys 
import os
from tabulate import tabulate
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from Methods.Volume import *


class NotState(Exception):
    pass



# Modificar para deixar na estrutura do feos
#Criar estado eqv. a iniciliar com solver de volume
#Esse estado possui metodos associados;
#-> Phi,Gres, Ares, Z

class State():
  

  def __init__(self,eos:CPA,vx:np.ndarray,Pressure:float,T: float, density_initialization: str):

    self.Eos = eos

    self.vx = vx
    self.T = T
    self.Pressure = Pressure
    self.phase = density_initialization

    self.rho = 1/Calc_Volume(eos,vx,T,Pressure,state=density_initialization)['x']
           

    self._R = 8.31446261815324



  def __str__(self) -> str:

    table = self.Table()

    return(table)

  def debug_cpa(self):
     
    X_prod , vphi = self.Eos.debug_cpa(self.vx,self.T,self.rho) # type: ignore

    print(
            f'------ Debug ------\n'
            f'Phase: {self.phase}\n'
            f'Density: { self.rho}\n'
            f'X_prod: {X_prod}\n'
            f'vPhi: {vphi}\n'

              )

    return X_prod , vphi
     
  def phi(self,)-> np.ndarray:
     
     vphi = self.Eos.phi(self.vx,self.T,self.rho)

     return vphi
  

  

  def Table(self):
     
    # X_prod , vphi = self.Eos.debug_cpa(self.vx,self.T,self.rho) #
    # 
    #  type: ignore

    X_prod = self.Eos.X_product(self.vx,self.T,self.rho)
    vphi = self.Eos.phi(self.vx,self.T,self.rho)

    # DeltaMatrixProd = self.Eos.DeltaMatrix(self.vx,self.T,self.rho).prod()


    phase = self.phase
    phi_comp1 = vphi[0]
    phi_comp2 = vphi[1]

    comp1_frac = self.vx[0]
    comp2_frac = self.vx[1]

    table=tabulate([["{} {}".format(phase.capitalize(),"Phase"),"Component","Phi","X_prod","Density"],np.array(["-","Comp1",phi_comp1,X_prod,self.rho]),["-","Comp2",phi_comp2,X_prod,self.rho]],
                  headers="firstrow",tablefmt='psql')
    

    return table


  

