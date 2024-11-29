
import sys 
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from typing import Callable
import numpy as np
import scipy.optimize as opt

from Bissection import fn_Alpha

from diff import *

#  args_0 = *State_0 

def newton(f: Callable ,args_0:list,f_obj:float,x_pos:int,a=1.):

    """"

    args_0: todos os argumentos necessários da função ; como é um método 
    unidimensional , apenas o argumento na posição dada por " x_pos "
    será modificado; 
    
    f_obj: o objetivo do método (se for 0, será encontrado a raiz no entorno de x0)

    
    
    """


    # Precisa ser cópia

    rtd_args = args_0*1

    # print(rtd_args)

    fn_dfdx = grad_1(f,pos=x_pos)

    f_x1 = f_obj


    x0 = rtd_args[x_pos]


    it_max = 10

    it = 0

    err=1000.

    tol = 1e-6

    success = True

    while err>= tol and it<it_max:
        

        
        # 1) Calc. f x (k-1)
        # print('args',rtd_args)

        f_x0 = f(*rtd_args)

        dfdx = fn_dfdx(rtd_args)

        res0 = f_x1 - f_x0

        h = res0 / dfdx


        result = fn_Alpha(F = f,args=rtd_args,h=h,x_pos=x_pos)

        Alpha , it_bis ,suc_bis = result

        # print(
        #     f'Iteração: { it }\n'
        #     f' x0: {x0}\n'
        #     f' Der: {dfdx}\n'
        #     f'Alpha: {Alpha}\n'
        #     f'Iteração biss: {it_bis}\n'
        #     f'Sucesso biss: {suc_bis}'
        #     f' h: {h}\n'
        #     f'Erro : {err}\n'
        #     )
        

        # 2) Atualiza x

        if np.isnan(Alpha):
            Alpha = 1


        x1 = x0 +  Alpha * h

        # 3) Atualiza a lista de argumentos na posição da variável x
        # do problema 1D f(x)=0

        rtd_args [x_pos] =  x1*1
        
        # 

        x0 = x1*1
        
        err = abs( res0 )

        it += 1




    if  (it==it_max):
    
        print('Não foi encontrada raiz na localidade x ')

        x1 = np.nan*1


        
    if np.isnan(x1):

        success = False

    # modificar pra avisar algum erro (algo dando nan)
    
    ans = [x1,it,success]

    return ans