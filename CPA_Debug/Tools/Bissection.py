from diff import *

# apenas pra função 1D
# pensado pra otimizar o Alpha F( x + alpha *h)
from math import copysign

def fn_Alpha(F:Callable,args:list,h:float,x_pos:int , a=1e-6,b=1):

    """
    a,b : intervalo de busca

    h : direção de busca

    """

    xval = args[x_pos]

    # função cuja raiz possibita o passo ideal na busca 
    # da raiz de uma função F no método de newton

    def f(t)->float:

        
        Alpha = -np.log(t)


        args[x_pos] = xval + Alpha*h

        # é necessário
        f = F(*args)



        return f



    res = 1

    tol=1e-5
    it_max = 100
    it= 0

    # Chute inicial

    t = (a+b)*.5

    
    # print('f(a)',f(a))
    # print('f(b)',f(b))

    succses = True

    while res>tol and it<it_max:

        it +=1



        # print(f't:{t}')

        
        f_val= f(t)

        f_inf = f(a)
        f_sup = f(b)

        # 1) Verifica se f(a) e f(b) possuem diferentes sinais

        print(f'c={t}')

        print(f'- ln(t) = {-np.log(t)}')
        print(f'fval={f_val}')
        print(f'fa={f(a)}')
        print(f'fb={f(b)}')
        print(f'a={a}')
        print(f'b={b}')
        print('iteração',it)

        if f_sup>0 and f_inf<0:
            # Verifica se f_val>0 ; se sim, b <-- c (diminui o intervalo)



            if f_val>0:

                b = t
            else:

                a = t
            

        elif f_sup<0 and f_inf>0:

            # Verifica se f_val<0 ; se sim, b <-- c (diminui o intervalo)

            if f_val<0:

                b = t
            else:

                a = t


        # elif f_val>0:

        #     b = t
        
        # elif f_val<0:

        #     a = t


        else:

            
            

            # print(f'c={t}')
            # print(f'fval={f_val}')
            # print(f'fa={f(a)}')
            # print(f'fb={f(b)}')
            # print(f'a={a}')
            # print(f'b={b}')
            # print('iteração',it)
            # print(f"Não há raízes nesse intervalo ( f(a) e f(b) ) posseum mesmo sinal.")

            t = np.nan*1

            succses = False

            break


        # print(f'c={c}')
        # print(f'fval={f_val}')
        # print(f'fa={f(a)}')
        # print(f'fb={f(b)}')
        # print(f'a={a}')
        # print(f'b={b}')

        
        # passo de newton pra atualizar f


        t = (a+b)*.5

        res = abs(f_val)

    
    
    if it==it_max:

        f_val = f_val*np.nan

        succses = False 

    
    Alpha = -np.log(t)



    return [Alpha,it,succses]

        