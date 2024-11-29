
#%%

import sys
import os




sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

# from Tools.Newton import *
from Tools.diff import *

from EquationOfState.CPA import *





_R = 8.31446261815324

# Modificado pro solver do michelssen
    
def Calc_Volume(Eos:cEOS|CPA,vx:np.ndarray,T:float,P:float,state: str,
                it_max=50,tol=1e-6):
    
    """
    Retorna dicionário:

    result['x'] : float ; volume,
    result['sucsses'] : bool ; se convergiu ou não,
    result['nit'] : int ; n de iterações.

    """
   

    def VolumeSolver(initial_guess):


        bm = Eos.f_b_mix(vx)
        
        rhomax = 1/bm


        def F ( s: float ):


            rho = s*rhomax

            # estado ja ta sendo definido aqui

            P_at_s = Eos.Pressure(vx,T,rho)
            

            F = ( 1 - s )*( P_at_s - P )

            return F




  # esse chute sera chute de s

        def solver(initial_guess):
            
            # limites 
            s_Min = 0
            s_Max = 1

            d = 1e-5

            TOL=tol
            #tolerância do resíduo

            MAX=it_max

            it=0 #contador de iteraçõesm iniciado como zero, incrementado em 1 após cada iteração concluíd

            d = 1e-5
            s1 = initial_guess

            F_0 = F(s1)


            # F_at_s1 = F_0*1.02

            res = 1

            F_at_s0 = F_0*1

            fn_dFds =  grad_1(F,pos=0)

            # print(f's0 = {s1}\n'
            #       f'F0 = {F_0}\n'
            #       f'rhomax = {rhomax}')



            fora=0

            while(np.abs(F_at_s0)>TOL and it<MAX): #Em caso de problema em que o método iterativo crie um loop infinit


                # print(f'ITERAÇÃO{it}')

                it += 1 #incremento no contador de iterações

                s0 = s1*1

                
                
                # 1) -----  Estabelece novos limites de busca ----- 

                F_at_s0 = F(s0)

                if F_at_s0 > 0 :


                    s_Max = s0*1

                elif F_at_s0 <0:
                

                    s_Min = s0*1


                # 2)  -----  Atualiza s ----- 

                dFds = fn_dFds([s0])

                # F_mais = F( s0 + d )
                # F_menos = F( s0 - d   )

                # dFds = 1*( F_mais-F_menos ) / (2*d)


                s1 = s0 - F_at_s0/dFds
                # print(f'Iteração = {it}\n'
                #       f's1 = {s1}')
                
                # 3) ----- Verifica se s_novo está no intervalo

                if s1>=s_Min and s1<=s_Max :

                # print('entre')

                    continue

                # elif s1<s_Min and s1>s_Max:
                else:

                    s1 = (s_Max + s_Min )/2

                # print('foraaaaa')

                # print('s1 apos novo ',s1)

                if it==MAX:

                    s1 = np.nan*1
            
            if np.isnan(s1):

                success = False

            else:
                success = True
            
            ans = {}

            # result = [s1,success,it]

            ans['x'] = s1
            ans['success'] = success
            ans['nit'] = it

            return ans


        
        ans = solver(initial_guess)

        # s_root = ans[0]
        s_root = ans['x']


        #  s = rho/rhomax

        rho = rhomax*s_root

        vol = 1/rho

        ans['x'] = vol

        return ans 

    bmix = Eos.f_b_mix(vx)

    if state=="liquid":

        liquid_density_initial =  (1/(bmix*1.01))
        
        liquid_zeta_initial = 0.99

        # print(liquid_density_initial)


        # State = np.array([vx,T,liquid_zeta_initial])
        zeta_0 =  liquid_zeta_initial

    elif state=="vapor":

        
        # vapor_density_initial =  1/(_R*T/P)

        vapor_zeta_initial = bmix/(bmix + (_R*T)/P )


        # State = np.array([vx,T,vapor_zeta_initial])

        zeta_0 = vapor_zeta_initial
    
    else:
        
        raise TypeError('Insira "liquid" ou "vapor".')
    


    # pressure_fn = Eos.Pressure

    # print(*State)

    # result = newton(pressure_fn,*State,f_obj=P,x_pos=2)

    
    # density_root = result[0]
    # it = result[1]

    # print('iterações',it)
    # Volume = 1/density_root

    result = VolumeSolver(initial_guess= zeta_0)


    return result


# def VolumeMeu(Eos:Ceos|Cpa,vx:np.ndarray,T:float,P:float,state: str,
#                 it_max=50,tol=1e-6):
    
#     bmix = Eos.f_b_mix(vx)

#     if state=="liquid":

#         liquid_density_initial =  (1/(bmix*1.01))
        
#         # liquid_zeta_initial = 0.99

#         # print(liquid_density_initial)


#         State = [vx,T,liquid_density_initial]
#         # zeta_0 =  liquid_zeta_initial

#     elif state=="vapor":

        
#         vapor_density_initial =  1/(_R*T/P)

#         # vapor_zeta_initial = bmix/(bmix + (_R*T)/P )


#         State = [vx,T,vapor_density_initial]

#         # zeta_0 = vapor_zeta_initial
    
#     else:
        
#         raise TypeError('Insira "liquid" ou "vapor".')
    


#     pressure_fn = Eos.Pressure

#     # print(*State)

#     result = newton(pressure_fn,State,f_obj=P,x_pos=2)

    
#     # density_root = result[0]
#     # it = result[1]

#     # print('iterações',it)
#     # Volume = 1/density_root

#     # result = newton()


#     return result




