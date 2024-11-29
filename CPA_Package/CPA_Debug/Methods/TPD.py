
#%%Imports

import os 
import sys

import matplotlib.pyplot as plt

import scipy.optimize as opt

# sys.path.append("E:\\CPA")

# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..','..')))

# print(sys.path)
# print(sys.path)

from Methods.Volume import *

# from Eos import *:
def TPD(eos:CPA|cEOS,T:float,P:float,y:np.ndarray,xguess:np.ndarray ,
        phasey:str,phasex:str, it_max = 200,tol=1e-10):
    
    # tol reduzida pra 1e-7
    # calc h 

    ansy= Calc_Volume(eos , vx = y , T=T , P=P , state = phasey )
    
    # Inserir alerta de que deu nan
    
    voly = ansy['x']

    rhoy = 1/voly

    phiy = eos.phi(vx=y , T = T, rho = rhoy)

    # print(phiy)

    hy = np.log(y) + np.log(phiy) 


    W = xguess*1
    x = xguess*1

    # print(f'Inicialização:\n'
    #       f'ansy = {ansy}\n'
    #       f'phiy = {phiy}\n'
    #       f'y= {y}\n'
    #       f'hy = {hy}\n')
    
    res = 1

    it_max = it_max
    it = 0

    success = True

    while ( res > tol and it<it_max):

        W_0 = W *1

        # calc lnphi_x

        # print('iteração', it)

        # print(phasex)
        # print(x)
        
        ans = Calc_Volume(eos,vx = x , T=T , P=P ,state = phasex)

        volx = ans['x']

        rhox = 1/volx

        # print(f"V filha = {ans}")

        
        
        # print(rhox)
        phix = eos.phi(vx = x , T =T, rho= rhox )

        lnphi_x = np.log(phix)

        # print(f'phix  = {phix}')

        # Atualiza W 

        W = np.exp( hy - lnphi_x )

        x = W / np.sum(W)

        res = np.linalg.norm(W - W_0)

        it += 1

    hx = lnphi_x + np.log(x)

    DG = np.sum( W * ( hx - hy )  )

    if it ==it_max:
       

       success = False

       DG = DG*np.nan
       x = x*np.nan

       print( f'Result TPD = { success} ')

    # sucsses flag for it ==itmax

    # print(f'DG = {DG}')
    # print(f'Iteração= {it}')
    # print(f'x = {x}')

    return DG , x


# primeira parte 

# 
def pure_water_dew_point(eos:CPA,T,P,y_dry,index_water,ywguess=2000*1e-6,
                         xw= 0.99999,method='lm',tol=1e-6)-> float :

  # Aproximação

  N = len(y_dry)

  def fn_res(yw0: list):

    # Fase vapor : mistura com água
    # Fase líquida : apenas água 

    # 1) Inicializa com um chute p\ fração de água;

    yw = 10**( yw0[0] )

    # print(yw)
    
    # 2) Converte p\ base úmida 

    y = y_dry/np.sum(y_dry)
    y=y*(1-yw)
    y[index_water]=yw

    # print('antes',y)

    y=y*(1-yw)
    y[index_water]=yw
    # print('depois',y)
    # print(f'yw inside pure = {yw}')
    # print(f'y = {y}')


  
    # Da água pura:

    x = np.ones(N)*(1-xw)/( N -1)

    x[index_water] = xw

    # print(f'x = {x}')
    

    # Da mistura:

    ansy_vap = Calc_Volume(eos,vx=y,T=T,P=P,state='vapor')
    voly_vap = ansy_vap['x']

    rhoy_vap = 1/voly_vap


    G_vap, phiy_mix_vap = eos.g_res_mix(vx = y ,T=T,rho=rhoy_vap)

    phiy_water = phiy_mix_vap[index_water]

    # Da água pura no liquido:

    # print(phiy_water)

    ansx = Calc_Volume(eos,vx= x,T=T,P=P,state='liquid')

    volx = ansx['x']

    rhox = 1/volx

    phix_pure_water = eos.phi(vx = x , T=T, rho = rhox)[index_water]


    res = np.log(yw) + np.log(phiy_water) - np.log(phix_pure_water)

    # print(res)

    return [res]
        

  # yw0 = Antoine(); implementar antoine?

  # Chamamento da função da função resíduo 


  yw0_log10 = np.log10(ywguess)

  # print('fn',fn_res([yw0]))

  # print('yw0',yw0,ywguess)
  ans = opt.root(fn_res, [ yw0_log10 ], method=method,)

  # print(f'Result Pure Water:{ans.success}')



  ywguess = 10**(ans.x[0])


  return ywguess

# implementar ydry

def Watcon(eos:CPA,T,P, y_dry , index_water,xwguess = 0.99999 , method=None,tol = 1e-6 ,iMAX_TPD=200, tol_TPD = 1e-9,
           phase_y='vapor',phase_x='liquid',)->list:
    

  ncomp = len(y_dry)

  # print(f'index_water = {index_water}')
  # Corpo função resíduo

  def fun(vetor):

  
    # 1) Inicializa com um chute p\ fração de água;
    
    yw=10**vetor[0]

    # 2) Converte p\ base úmida 

    y = y_dry/np.sum(y_dry)

    # print('antes',y)
    y=y*(1-yw)
    y[index_water]=yw
    # print('depois',y)


    # print(f'y vetor inside watcon  = {y}')
    # print(f'ydry vetor inside watcon  = {y_dry}')
    # print(f'yw inside watcon = {yw}')
    
    # y = np.array([ yw,1-yw] )

    # ymod = np.linalg(y)

    # y = y/ymod

    xguess= np.ones(ncomp)*(1-xwguess)/(ncomp-1)
    xguess[index_water] = xwguess

    delta_g , x = TPD(eos,T,P, y = y , xguess= xguess ,phasey=phase_y , phasex= phase_x,it_max=200)


    return [delta_g]
  
  # Chamamento da função resíduo

  ywguess = pure_water_dew_point(eos,T,P,y_dry,index_water=index_water)

  # print(ywguess)
  # print(f'ywguess = {ywguess}')

  # v0 = np.log10(ywguess)
  ywguess_log10 = np.log10(ywguess)

  # print(f'ywguess_log= {ywguess}')

  ans = opt.root( fun, [ywguess_log10], )

  print(f'Result Watcon:{ans.success}')

  yw = 10**ans.x[0]

  suc = ans.success


  ywet = y_dry*(1-yw)

  ywet[index_water] = yw

  ywet_norm = ywet/np.sum(ywet)

  return [ yw , ywet_norm ] 







