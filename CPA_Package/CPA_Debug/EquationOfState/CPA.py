import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from AscEos import *

from Ceos import *




class CPA():


    def __init__(self,parameters):


        self.p = parameters

        self.Ceos = cEOS(parameters)
        self.Asc = ascEOS(parameters)

        # cubic

        self.va0=self.Ceos.va0
        self.vb=self.Ceos.vb
        self.vc1=self.Ceos.vc1
        self.vtc=self.Ceos.vtc

        self.kij=self.Ceos.kij

        self.ncomp=self.Ceos.ncomp

        self.sig=self.Ceos.sig
        self.eps=self.Ceos.eps

        # assoc

        self.scheme= self.Asc.scheme
        self.vepsilon=self.Asc.vepsilon
        self.vbeta=self.Asc.vbeta

        # self.Sij_matrix =  self.S_matrix()


        # print(self.beta_cross)

        self._R = 8.31446261815324

    def __str__(self):
        # Configuração para exibir arrays numpy com precisão e sem notação científica
        np.set_printoptions(precision=4, suppress=True)

        return (
            f"{'==='*10} Associative Parameters {'==='*10}\n"
            f"Scheme Type               = {np.array_str(self.Asc.scheme)}\n"
            f"vBeta                     = {np.array_str(self.Asc.vbeta)}\n"
            f"vEpsilon                  = {np.array_str(self.Asc.vepsilon)}\n"
            f"Constant R                = {self.Asc._R:.4f}\n\n"
            
            f"{'==='*10} Cubic Parameters {'==='*10}\n"
            f"kij                       = {np.array_str(self.Ceos.kij,)}\n"
            f"va0                       = {np.array_str(self.Ceos.va0)}\n"
            f"vb                        = {np.array_str(self.Ceos.vb)}\n"
            f"vc1                       = {np.array_str(self.Ceos.vc1)}\n"
            f"vtc                       = {np.array_str(self.Ceos.vtc)}\n"
            f"sig                       = {self.Ceos.sig}\n"
            f"eps                       = {self.Ceos.eps}\n"
            f"Constant R                = {self.Ceos._R:.4f}\n"
        )

    def f_b_mix(self,vx,):


        vb = self.vb

        bmix = f_b_mix(vx,vb)

        return bmix
    
    def Pressure(self,vx,T,rho):



        P_asc = self.Asc.Pressure(vx,T,rho)
        P_cub = self.Ceos.Pressure(vx,T,rho)
        
        # print(f'Pasc:{P_asc}')
        # print(f'Pcub:{P_cub}')

        # if np.isnan(P_asc):
            # P_asc = 0
            # devido a nao convergencia do Xassoc
            #  pra que o newton convirja até uma raiz boa de volume
            # pra que o xassoc funcione
            
        P_total = P_asc + P_cub
                
        return P_total
    
    def phi(self,vx,T,rho)->np.ndarray: 


        Pcpa = self.Pressure(vx,T,rho)
        
        Pcub = self.Ceos.Pressure(vx,T,rho)

        phiassoc = self.Asc.phi(vx,T,rho)
        
        phicub = self.Ceos.phi(vx,T,rho)


        phi_cpa = phicub*phiassoc*(Pcub/Pcpa)

        return phi_cpa
    
    def g_res_partial(self,vx,T,rho):
        
        _R= self._R


        vphi = self.phi(vx,T,rho)

        vlnphi = np.log(vphi)

        vGres = (_R * T) * vlnphi


        return vGres , vphi
    
    def g_res_mix(self,vx,T,rho):
        
        vgres , vphi= self.g_res_partial(vx,T,rho)

        G_res_mix = np.sum(vgres*vx)

        return G_res_mix , vphi
        
    def X_product(self,vx,T,rho):


        X_prod = (self.Asc.X_assoc(vx,T,rho)).prod()

        return X_prod
    
    def DeltaMatrix(self,vx,T,rho):

        delta = self.Asc.delta_AB(vx,T,rho)
        
        return delta
    
     
    

