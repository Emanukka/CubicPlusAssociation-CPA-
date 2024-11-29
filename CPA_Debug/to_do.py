


# Modificar a parte associativa p ncomp>2;
# solvatação é função do comp. i e j
# valores de beta cross e eps cross também ( no caso que há solvatação)
# talves sera interessante pedir pra scheme ser passado como dicionario,
# tal como o kij
# scheme = [ 4C , S , S]

# se componente é S, entao solvata

# dai , é verificado nos dicionarios

# eps_cross = {eps_cross12: k , eps_cross13 : l}

#%%

import numpy as np


scheme_type = ['4C','S21','S31']

# a ordem dos indexes devem ser os mesmos

eps_cross = {'12':20,'13':10}
beta_cross = {'12':0.2,'13':.1}

# def EpsCrossMatrixConstructor

def DeltaCross(BetaCross:dict,EpsCross:dict,ncomp:int):
  
  n = ncomp


# 1) Transformar os dicionarios na forma de chaves de index apenas


#   keysEps = list(EpsCross.keys())
#   keysBeta = list(BetaCross.keys())
#   sepEps = 'EpsCross'
#   sepBeta = 'BetaCross'
 
#   EpsCrossNew = {}
#   BetaCrossNew = {}

#   for EpsCrossKey in keysEps:
     
#      idx = EpsCrossKey.split(sepEps) [-1]

#      EpsCrossNew[idx] = EpsCross[EpsCrossKey]

#   print(EpsCrossNew)
#   for BetaCrossKey in keysBeta:
     
#      idx = BetaCrossKey.split(sepBeta) [-1]

#      BetaCrossNew[idx] = EpsCross[BetaCrossKey]


#   print(EpsCrossNew,BetaCrossNew)

   


#   {key:EpsCross.values['EpsCross'] for key in EpsCrossKey.split(sep) [-1]}



  # vals = [ a[1:] for  ]

#   index_ij_array = np.zeros_like(keys)
#   print(index_ij_array)


  Aij = np.zeros((n,n))

#   for idx,EpsCrossKey in enumerate(keys):
    
#     #   print(EpsCrossKey.split(sep) [-1])
#       index_ij = EpsCrossKey.split(sep) [-1]

    
#     #   print(index_ij_array)
#     #   print(index_ij)
#       index_ij_array[idx] = index_ij

#       i,j= int(index_ij[0]) , int(index_ij[1])

#     #   print(Aij)
#       Aij[i-1,j-1] = Aij[j-1,i-1] = EpsCross[EpsCrossKey]



  return 0


aij = DeltaCross(beta_cross,eps_cross,3)

def Delta(ncomp,scheme_type,epscross:dict):

    keys = epscross.keys()

    # keys = epscross.keys()

    EpsMatrix ,index_ij_array= EpsCrossMatrixConstructor(eps_cross,3)

    print(EpsMatrix)

    for i in range(ncomp):
      
      SCHEME_i = scheme_type[i]

      #  for SCHEME_i in scheme_type:

          # i solvatado por outro componente j

      if SCHEME_i[0] == 'S':

        ij = SCHEME_i[1:]


        print('solv')

        print(index_ij_array)


        # garante que nao seja varrido 0

        if ij in index_ij_array or ij[::-1] in index_ij_array:
           
           i,j = int(ij[0]), int(ij[1])

           epscross_ij = EpsMatrix[i-1,j-1]

           print(epscross_ij,i,j)
        
        else:
           
           print('nao solvata',ij)
           

           

           


        # Aqui ja extrai o index pelo dicionario

        # fazer match de index 
        

        # print('SOLVATA')

        # # DeltaCrossarray



        # for j in range(ncomp):
           
        #    if i!=j:
              

        #       EpsCross=EpsMatrix[i,j]

        #       print(f'EpsCross= {EpsCross},i={i},j={j}')
              

# Delta(3,scheme_type,eps_cross)
           





         
        
             

             

