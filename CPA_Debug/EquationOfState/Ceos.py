
import numpy as np
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# from EosParamns import 

from EquationOfState.CoreEquations.CubicEquations import *




# Pasta - > Modulo - > Classe
     
class cEOS():

  def __init__(self, parameters:dict ):

    Keys = ['va0','vb','vc1','vtc']

    # i= 0 

    for k in Keys:
       if isinstance(parameters[k],np.ndarray):
          continue
          
       elif isinstance(parameters[k],(list,tuple,float)) is False:
          
          raise ValueError("Arguments must be of type np.ndarray , list , tuple or float ;")
       
       else:
          
          parameters[k] = np.array(parameters[k])



    # Array like objects
    self.va0=parameters['va0']
    self.vb=parameters['vb']
    self.vc1=parameters['vc1']
    self.vtc=parameters['vtc']


    # Int obj
    self.ncomp=parameters['ncomp']
    self.sig=parameters['sig']
    self.eps=parameters['eps']


    


    if 'kij' in parameters:

      KijDict = parameters['kij']
       
      if isinstance(KijDict,dict):
         

         
         KijMatrix = KijMatrixConstructor(KijDict,self.ncomp)

         self.kij = KijMatrix

      else:
         
         raise TypeError("kij argument must be a dict object.")
      
    else:
       

       self.kij = np.zeros((self.ncomp,self.ncomp))



    self._R = 8.31446261815324

    return

  def __str__(self):
      
      return(
      f"=== Cubic Parameters ===\n"
      f"va0 = {self.va0}\n"
      f"vb = {self.vb}\n"
      f"vc1 = {self.vc1}\n"
      f"vtc = {self.vtc}\n"
      f"kij = {self.kij}\n"
      f"sig = {self.sig}\n"
      f"eps = {self.eps}\n"
      f"Constant R = {self._R}\n"
      )

  def f_dadn(self,vx,T:float,rho:float):

    kij= self.kij
    ncomp = self.ncomp

    

    va= cEOS.f_a_soave(self,T)
    am= cEOS.f_a_mix(self,vx,T,rho)

    dadn=np.zeros(ncomp)


    dadn = calculate_dadn(vx,T,rho,va,am,kij,ncomp)

    # vetor

    return dadn
  
  def f_a_soave(self,T:float):


    va0= self.va0
    vc1 = self.vc1
    vtc = self.vtc

    va = calculate_a_soave(T,va0,vc1,vtc)

    return va

  def f_b_mix(self,vx):


      vb = self.vb

      bmix = calculate_b_mix(vx,vb)

      return bmix

  def f_a_mix(self,vx,T:float,rho:float):

      #(section 3.1, eqs. 3.2 e 3.3; Thermod. M. for Ind. Ap.)
      #kij=>parâmetro de interação binaria: interação i com j, onde k[i,j]=k[j,i]
      # para o resto dos elementos, k=0 (delta de kroniger?)
      
      V=1/rho


      kij= self.kij
      ncomp = self.ncomp
  

      # Calc

      va = cEOS.f_a_soave(self,T)

      am = calculate_a_mix(vx,T,rho,va,kij,ncomp)
      return am


  def Pressure(self,vx,T:float,rho:float)-> float:

    # v = 1/self.s.rho

    # print(v)

    # V=1/self.s.rho
    # rho = self.s.rho


    

    # va0= self.p.va0
    # vc1 = self.p.vc1
    # vtc = self.p.vtc
    # kij= self.p.kij
    # ncomp = self.p.ncomp
    # vx = self.s.vx
    # vb = self.p.vb

    sig = self.sig
    eps= self.eps

    _R= self._R
    

    am= cEOS.f_a_mix(self,vx,T,rho)

    bm= cEOS.f_b_mix(self,vx)

    P = calculate_pressure(vx, T, rho, am, bm, _R, sig, eps)

    return P
  
  # def A_res(self,rho):
     
  #   bm = cEOS.f_b_mix(self,)
  #   am =cEOS.f_a_mix(self,)

  #   R = self.p._R
  #   T = self.s.T 

  #   A=  -np.log(1-bm*rho) - am/((_R*T*)) 
     
     

  def phi(self,vx,T:float,rho:float)->np.ndarray:
        # v = 1/self.s.rho

    # print(v)
    
    V=1/rho

    ncomp = self.ncomp
    vb = self.vb

    sig = self.sig
    eps= self.eps

    _R= self._R
    


    am= cEOS.f_a_mix(self,vx,T,rho)
    bm= cEOS.f_b_mix(self,vx)
    da= cEOS.f_dadn(self,vx,T,rho)


    rho=1/V

    db=vb*1


    P = cEOS.Pressure(self,vx,T,rho)

    phi = calculate_phi(vx, T, rho, am, bm, da, db, P, _R, sig, eps, ncomp)

    return phi
  
def KijMatrixConstructor(Kij:dict,ncomp:int):
  n = ncomp
  keys = Kij.keys()
  # vals = [ a[1:] for  ]
  Aij = np.zeros((n,n))
  for kij_key in keys:
      
      index_ij = kij_key[1:]
      i,j= int(index_ij[0]) , int(index_ij[1])
      Aij[i-1,j-1] = Aij[j-1,i-1] = Kij[kij_key]

  return Aij
 