

import numpy as np
# import jax.numpy as jnp

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from EquationOfState.CoreEquations.AscEquations import *

# from EosParamns import 


# dtype = complex


class ascEOS():

  def __init__(self, parameters:dict,):

    Keys = ['vb','vepsilon','vbeta','scheme']

    
    # Passando p\ numpy array;

    for k in Keys:
       if isinstance(parameters[k],np.ndarray) is True:
          continue
          
       elif isinstance(parameters[k],(list,tuple,float)) is not True:
          
          raise ValueError("Arguments must be of type np.ndarray , list , tuple or float ;")
       
       else:
          
          parameters[k] = np.array(parameters[k])
    
    # for scheme_type in parameters['scheme']:

    #    if scheme_type == 'inert' and parameters['beta_cross'] == None :
    
    if 'solvation' in parameters['scheme']:
       
       if 'beta_cross' and 'epsilon_cross' not in parameters:
          
          raise ValueError("solvation entre comp1 e comp2 necessita de parâmetros cruzados 'beta_cross' e 'epsilon_cross.")
      
       else:
          
          self.beta_cross=parameters['beta_cross']
          self.epsilon_cross=parameters['epsilon_cross']

          # print("Solvatação Ok")
           
 
      
           
    self.vb= parameters['vb']
    self.scheme= parameters['scheme']
    self.ncomp=parameters['ncomp']
    self.vepsilon=parameters['vepsilon']
    self.vbeta=parameters['vbeta']

    # self.Sij_matrix =  self.S_matrix()


    # print(self.beta_cross)

    self._R = 8.31446261815324
    # self.Xassoc = ascEOS.X_assoc(self,)

    # Xassoc fica constante e impossibilita calcular derivada
    
    

    return
  
  def __str__(self):
      return(
      f"=== Associative Parameters ===\n"
      f"Scheme Type = {self.scheme}\n"
      f"vBeta = {self.vbeta}\n"
      f"vEpsilon =  {self.vepsilon}\n"
      


      # Constante
      f"Constant R = {self._R}\n"
      )
  
  def f_b_mix(self,vx,):


    vb = self.vb

    bmix = f_b_mix(vx,vb)

    return bmix
  


  def f_g(self,vx,rho:float):


    vb = self.vb     

    gmix = f_g(vx,rho,vb)


    return gmix
  


  def f_dlngdrho(self,vx,rho:float)->float:

    vb = self.f_b_mix(vx)
    dlngdrho = f_dlngdrho(vx,rho,vb)

    return dlngdrho 

  def f_dlngdni(self,vx,rho:float)->np.ndarray:  #$$$

    #method of yakoumis, 1998, eq A3, A4 and A2, with g(...) from kontogeorgis 1999 "s"cpa

    ncomp = self.ncomp
    vb = self.vb

    dlng = f_dlngdni(rho=rho,vx=vx,vb=vb,ncomp=ncomp)

    return dlng 


  def delta_AB(self,vx,T,rho,):

    vepsilon = self.vepsilon
    vb = self.vb
    vbeta = self.vbeta
    ncomp = self.ncomp


    scheme_type = self.scheme

    # association_type = self.p.association_type
    # eps_cross = self.p.eps_cross
    # beta_cross = self.p.beta_cross

    _R = self._R


    #  Se True, calcula delta_cross a partir f_delta(epscross,betacross)
    # if scheme_type == 'solvation':
       
      #  delta_matrix = delta_AB_Solvation(vx,T,rho,vepsilon,vb,vbeta,epsilon_cross,beta_cross,_R,ncomp)
    
    vdelta = f_vdelta(vx,T,rho,veps=vepsilon,vb=vb,vbeta=vbeta,_R= self._R) 


    # if SolvationEffect == True:
    # 
    #   for i in range(ncomp):
    #     
    #   scheme_type
    # 
    # 

    # DeltaCross passará a ser vetor ( matrix no caso )
    # 
    # 

    #for i in range(ncomp):
    #   
    #   if SolvationEffet==True:
    #     
    #     
    # 
    #
    #
    # 

    # for i in range(ncomp):
      
    #   SCHEME_i = scheme_type[i]

    #   #  for SCHEME_i in scheme_type:

    #       # i solvatado por outro componente j

    #   if SCHEME_i == 'S':
         
    #      for j in range(ncomp):


    #      eps_cross_ij =
             

             
    if 'solvation' in scheme_type:

      #  calc delta cross by ECR
      # print(f'Há solvatação')
      eps_cross = self.epsilon_cross
      beta_cross = self.beta_cross

      delta_cross = f_delta_cross(vx, T, rho, eps_cross, vb, beta_cross, _R)
      
      # print(delta_cross,'deltacross')

    else:
       #  calc delta cross by ECR
      # print(f'nao h á solvatação')

      delta_cross = f_delta_cross_ECR(vdelta)


    delta = delta_AB(vdelta,delta_cross,ncomp=ncomp)

    return delta
  
  def S_matrix(self,):

      #A,0 negativo ; B,1 positivo
      # A =0 , B = 1
      #Sji, where j -> site and i for component
      #scheme_tpye must be an vector;




      scheme = self.scheme

      ncomp = self.ncomp


      A = 0
      B = 1

      S = np.zeros((2,ncomp))


      for comp_idx,scheme_type in enumerate(scheme):
          
          # print(comp_idx,scheme_type)

          if scheme_type == '1A':
              S[B][comp_idx] = 1
              
          elif scheme_type == '2B':
              S[A][comp_idx] = S[B][comp_idx] = 1
              
          elif scheme_type == '3B':
          
              S[A][comp_idx] = 2
              S[B][comp_idx] = 1

          elif scheme_type == '4C':

              S[A][comp_idx] = S[B][comp_idx] =  2

          
          # S for Solvation between component i and j

          elif scheme_type == 'solvation':


            #  1 sitio negativo

            S[A][comp_idx] =  1

            #  1 sitio positivo
            
            # S[B][comp_idx] =  1
             


          elif scheme_type == 'inert':
             
             continue
          
          else:
             
             raise ValueError("Insira um esquema adequado : '1A' ; '2B ; '3B'  ; '4C' ou 'Sij' . ")
        

      # print('Smatrix',S)

      # return solvation effect

      return S 
  
  def X_assoc(self,vx,T,rho):

    # S[j,i]: j-> sítio positivo ou negativo ; i-> componente

    S_matrix = self.S_matrix()
    ncomp = self.ncomp
 
    delta = self.delta_AB(vx,T,rho)
  
    X_ji = X_assoc(vx,rho,delta,S_matrix,ncomp)

    return X_ji

  def Pressure(self,vx,T,rho) :
    
    # rho = self.s.rho
    _R = self._R
    vb = self.vb


    fAscMolec= ascEOS.assoc_molecules_frac(self,vx,T,rho)
    
    Passoc = Pressure(vx,T,rho,_R,vb,fAscMolec)

    return Passoc 

  def assoc_molecules_frac(self,vx,T,rho):


    S = self.S_matrix()
    ncomp = self.ncomp

    Xassoc = self.X_assoc(vx,T,rho)

    sum1 = assoc_molecules_frac(vx,T,rho,Xassoc,S,ncomp)

    return sum1

  def phi(self,vx,T:float,rho:float)-> np.ndarray :

    ncomp = self.ncomp


    S= self.S_matrix()


    fAscMolec= ascEOS.assoc_molecules_frac(self,vx,T,rho)

    Xassoc = ascEOS.X_assoc(self,vx,T,rho)


    dlng = ascEOS.f_dlngdni(self,vx,rho)
    
    vphi = f_phi(vx,T,rho,Xassoc,S,ncomp,fAscMolec,dlng)

    return vphi

 