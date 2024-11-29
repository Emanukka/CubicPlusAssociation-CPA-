import numpy as np
import numba

@numba.jit
def f_b_mix(vx, vb):

    return np.sum(vx * vb)
@numba.jit
def f_g(vx, rho:float, vb:np.ndarray)-> float:
    bm = f_b_mix(vx, vb)

    # print('rho',rho)
    # print('bm',bm)

    return 1 / (1 - 1.9 * (rho * bm / 4))

@numba.jit
def f_dlngdrho(vx, rho, vb):

    g = f_g(vx, rho, vb)

    # print('rho',rho)
    # print('g',rho)

    bm = f_b_mix(vx, vb)
    
    return (1 / g) * (-1) * g**2 * (-1.9 * bm / 4)

@numba.jit
def f_dlngdni(rho,vx,vb,ncomp: int):

  # na dedução , o termo de 1/n que aparece na derivada  dlngdni
  # cancela quando multiplica dlngdni com h
  # como esperado, ja que é adimensional

  g_mix = f_g(vx=vx,rho=rho,vb=vb)


  v_dlngdn_i = g_mix * 1.9 * (rho*vb/4)

  return v_dlngdn_i

# def f_dlngdrho_i(vx: np.ndarray, rho: float, vb: np.ndarray, ncomp: int) -> np.ndarray:

#     g_mix= f_g(vx, rho, vb)

#     dgdni=np.zeros(ncomp)
#     dgdxi=np.zeros(ncomp)


#     #g=1/(1-1.9*(y/4)); y=
#     dgdxi[:]= (-1)*g_mix**2 * (-1.9*(rho*vb[:]/4))


#     for i in range(ncomp):
#         sxdgdx=0
#         for j in range(ncomp):
#             sxdgdx += dgdxi[j]* vx[j] #incluindo j=i
#             #no caso da regra classica, os x[i]b[i] vao se recombinar em igual a rho dgdrho (usando bm)
#         dgdni[i]=dgdxi[i]#-sxdgdx+rho*dgdrho
# #         print('sxdgdx,rho*dgdrho',sxdgdx,rho*dgdrho, ) #igual OK

#     dlngdrhoi=(1/g_mix)*dgdni[:] *1/(rho)

#     return dlngdrhoi

    
@numba.jit

def f_vdelta(vx, T, rho, veps, vb, vbeta, _R) -> np.ndarray:

    gmix = f_g(vx, rho, vb)
    
    return gmix * (np.exp(veps / (_R * T)) - 1) * vb * vbeta


@numba.jit

def f_delta_cross(vx:np.ndarray, T:float, rho:float, eps_cross:float, vb:np.ndarray, beta_cross:float, _R:float) -> float:

    gmix = f_g(vx, rho, vb)
    
    bmix = f_b_mix(vx,vb)
    
    delta_cross = gmix * (np.exp(eps_cross / (_R * T)) - 1) *  bmix* beta_cross

    return delta_cross

@numba.jit

def f_delta_cross_ECR(vdelta: np.ndarray) -> float :
#  Elliott’s Combining Rule
    delta_cross = np.sqrt( vdelta[0] * vdelta[1] )

    return delta_cross



# @numba.jit

# F (vdelta , delta_cross)
# vai funcionar apenas comp 
def delta_AB(vdelta:np.ndarray,delta_cross:float,ncomp,):

    delta = np.zeros((2, ncomp, 2, ncomp))  # Matriz de 4 dimensões


    # vdelta = np.zeros(ncomp) 

    # vdelta[:]=f_vdelta(vx, T, rho, vepsilon, vb, vbeta, _R)

    # if (association_type == 0 or association_type == None):

    #     print('tipo 0')

    # delta_cross = np.sqrt(vdelta[0] * vdelta[1])

    # if(association_type == 1):

        # print('tipo 1')
        
        # delta_cross = f_delta(vx,T,rho,eps_cross,vb,beta_cross,_R)

    # A: negativo ; B: positivo 

    A = 0 
    B = 1
    for comp1 in range(ncomp):
        for comp2 in range(ncomp):
            # Auto-associação
            if comp1 == comp2:
                delta[A, comp1, B, comp2] = vdelta[comp1]
                delta[B, comp1, A, comp2] = vdelta[comp1]

            # Associação Cruzada
            
            else:
                # print('deltacross',delta_cross)
                delta[A, comp1, B, comp2] = delta_cross
                delta[B, comp1, A, comp2] = delta_cross

    return delta

@numba.jit
def X_assoc(vx,rho,delta,S_matrix,ncomp):

    #método de substituição sucessivas p\ det. fração X_ji
    # S[j,i]: número de sítios j presentes na molécula i
    # x[i]: fração molar do componente/molécula i


    S = S_matrix

    #Matriz X_A
    X_A=np.zeros((2,ncomp))    #linhas=> sítios;
    #Matriz residual (X_A_new - X_A_old)
    dif=np.zeros((2,ncomp))
    #Chute p\ iniciar a iteração do método



    X_A[:,:]= .5

    it_max=1000
    it=0
    tol=1e-7

    res1=1.



    F = 0

    # print(delta)

    while (res1>tol and it<it_max):
        it+=1

        X_A_old=X_A*1 #cria cópia
        for i in range(ncomp):
            for j in range(2):
            # Fazer range(0,1) para reduzir n de iterações?
                #if S[j,i]==0:
                # X_A[j,i]=1  #não há sítios associativos j em i
                    #dif[j,i]= np.abs(X_A[j,i]-X_A_old[j,i])
                    #continue # vai para j+1 (uma vez que j não irá interagir c\ sítios de outras moléculas)
                sum1=0.
                for k in range(ncomp):
                    sum2=0.
                    for l in range(2):
                            #continue
                        
                        # print('Sji',S[l,k])
                        # print('Xold',X_A_old[l,k])
                        # print('delta',delta[l,k,j,i])


                        sum2= sum2 + S[l,k]*X_A_old[l,k]*delta[l,k,j,i]
                        # print(delta[l,k,j,i],l,j)
                        # Se S[l,k]==0 , delta[k,l,i,j]=0
                    # print(sum1)

                    sum1= sum1 + vx[k]*sum2

                    # print('sum1',sum1)
                    # print('sum2_convergido',sum2)

                X_A[j,i]= 1./(1.+rho*sum1)
                dif[j,i]= np.abs((X_A[j,i]-X_A_old[j,i])/X_A[j,i])

        res1=np.max(dif) #todas iterações feitas => matriz X_A pronta
        
        # if ( res1 > res0 and t < 16) :

        #     t = t +1


        #     print(2**-t )
        #     Alpha = 2**-t 

        #     print('alpha value',Alpha,t)

        #     X_A[:,:]= Alpha
        
        # if t>16:

        #     X_A[:,:] = 1

        #     break

        
        # else:
        #     continue

        # # print('F',F)

        # # if (res>res0) and F!=10:

        # #     F+= 1 

        # # elif F==10:
        # #     # ajusta o chute
        # #     # mas ai eu ajusto pra onde? 
        # #     # pra cima ou pra baixo?

        # #     F = 0
        # #     # agua
        # #     X_A=np.zeros((2,ncomp)) 


        # #     X_A[:,0] = .7

        # #     # co2

        # #     X_A[0,1] = .7

        # #     X_A[1,1] = 1

        # #     print(X_A)
        # #     res = 1
        # #     # res0 = 1

        # # else:
        # #     continue

        # print('res',res)
        if it==it_max:

            # não houve convergência do método
            # print('Too many steps in X_assoc')

            # X_A[:,:]=1
            X_A[:,:]=np.nan*1

            # break
    
    # print('X',X_A)
    # print(it,'XASSOC')

    return X_A
@numba.jit
def assoc_molecules_frac(vx: np.ndarray, T: float, rho: float, Xassoc: np.ndarray, S: np.ndarray, ncomp: int) -> float:
    
    sum1 = 0.
    # equivale ao h
    #  
    for i in range(ncomp):
        sum2 = 0.
        for j in range(2):

            sum2 += (1 - Xassoc[j, i]) * S[j, i]
        sum1 += vx[i] * sum2

    # print(Xassoc)
    return sum1

@numba.jit
def Pressure(vx,T,rho,_R,vb,fAscMolec) :
    
    # rho = self.s.rho

    
    dg=f_dlngdrho(vx,rho,vb)




    # print('test')
    # print('fAscMolec:',fAscMolec)
    
    
    Passoc=-_R*T*((1/(2*(1/rho)))*((fAscMolec)*(1+(1/(1/rho))*dg))) 
        #kontogeorgis 2006;

    # except:
    #     Passoc = 0

    # Passoc=-_R*T*((1/(2*(1/rho)))*((fAscMolec)*(1+(1/(1/rho))*dg))) #kontogeorgis 2006;
                                                  #eq. 9.1, section 9.2 (Thermod. M. for Ind. Ap.)

    return Passoc 

@numba.jit
def f_phi(vx: np.ndarray, T: float, rho: float, Xassoc: np.ndarray, S: np.ndarray, ncomp: int, fAscMolec: float, dlngdni: np.ndarray) -> np.ndarray:
    lnphi = np.zeros(ncomp)
    for i in range(ncomp):
        slogX = np.sum(np.log(Xassoc[:, i]) * S[:, i])

        lnphi[i] = slogX - (fAscMolec / 2) * dlngdni[i]

    return np.exp(lnphi)
