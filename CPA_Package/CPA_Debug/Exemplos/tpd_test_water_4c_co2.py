
#%%imports

import sys 

# sys.path.append(r"D:\\CPA\\CPA_mod")


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


# from Methods.Volume import *



from Methods.TPD import *
# import Methods.TPD as a

from Exemplos.Water_4C_CO2 import cpa_eos

#%%Teste TPD



T = 300

P = 1* 1e5

yw=0.4

z=np.array([ yw,1-yw] )

wguess = np.array([.999,.001])


DG , j = TPD(eos=cpa_eos,T=T,P=P,y=z,xguess=wguess,phasey='vapor',phasex='liquid',it_max=5)

print(f'DG = {DG}\n'
      f'xfrac = {j}')



#%%watcon


T= 304
P = 20*1e5
y_dry = np.array([0.,1.])
y = Watcon(cpa_eos,T,P,y_dry,index_water=0)

# ta batendo com o meu do colab

# 0.0025326795832623246

# print(y)







#%%Plot 373 K
data = [ 25.33, 2406, 50.66, 1696, 101.33, 3853, 202.65, 4445, 405.3, 5036, 506.63, 5076, 531.96, 5023, 557.29, 5050]

# # Reshape the data into a 2-column array
a1 = np.reshape(data, (len(data) // 2, 2))

# # Convert the second column (water fraction) from ppm to percentage
a1[:, 1] = a1[:, 1] * 1e-6*100  

l1 = np.array([	11.07	,88020 * 1e-6 * 100])
l2 = np.array([24.26,	46910* 1e-6 * 100])
l3 = np.array([70.88,	21020* 1e-6 * 100])
l4 = np.array([102.35	,17470* 1e-6 * 100])
l5= np.array([133.52,	16750* 1e-6 * 100])
l6 = np.array([170.7	,17120* 1e-6 * 100])

vT = np.array([304.15,348.15,373.15])

v1 = np.array([1.01,25.33,101.33,111.46,126.66,151.99,202.65,405.3,607.95,709.28])

v2 = np.array([318230,11243,8758,8574,9034,10099,11967,13939,14728,14794])*1e-6*100

a2 = np.array([v1,v2]).T

vp = np.linspace(1,700,50)

a0 = np.array([l1,l2,l3,l4,l5,l6])


exp = [a1,a2,a0]

y_dry = np.array([0.,1.])

y_list = []

colors = ['blue','red','green']

for t in vT:
    # print('t',t)
    
    yw = np.zeros(len(vp))

    for index,p in enumerate(vp):

        p = p*1e5


        yw[index] , _= Watcon(cpa_eos,t,p,y_dry,index_water=0)

    y_list.append(yw)


    

for j in range(len(vT)):
    print(j)

    T = vT[j]




    # results_yw=yw_results
    yw = y_list[j]
    
    yw =yw*100

    # plt.figure()
    # plt.scatter(results_p,results_yw)




    plt.plot(vp,yw,label = f'{T}K' , c=colors[j])

    a = exp[j]
    
    pexp = a[:,0]
    yexp = a[:,1]

    # plt.scatter(results_p,results_yw,label=f'T={vT[j]}K')
    plt.scatter(pexp,yexp, c="black")

    # plt.xlim(0,550)

    plt.ylim(0,5)

    plt.xlim(vp[0],vp[-1])


    

    plt.xlabel(" P (bar)")

    # plt.title(f"phi phi agua pura no liquido")

    plt.ylabel(" Water mol fraction ( % mol )")
    plt.legend()


# plt.savefig(r'D:\CPA\CPA_mod\Exemplos\Plots\water_4C_373_2B.png')

plt.show()

#%%Plot 304 K


data = [ 25.33, 2406, 50.66, 1696, 101.33, 3853, 202.65, 4445, 405.3, 5036, 506.63, 5076, 531.96, 5023, 557.29, 5050]

# # Reshape the data into a 2-column array
a = np.reshape(data, (len(data) // 2, 2))

# # Convert the second column (water fraction) from ppm to percentage
a[:, 1] = a[:, 1] * 1e-6*100  


vT = np.array([304.15])
              

vp = np.linspace(1,750,100)

# a = np.array([l1,l2,l3,l4,l5,l6])
pexp = a[:,0]
yexp = a[:,1]


y_list = []


for t in vT:
    # print('t',t)
    
    yw = np.zeros(len(vp))

    for index,p in enumerate(vp):

        p = p*1e5


        yw[index] = Watcon(cpa_eos,t,p)

    y_list.append(yw)


    

for j in range(len(vT)):
    print(j)

    T = vT[j]




    # results_yw=yw_results
    yw = y_list[j]
    
    yw =yw*100

    # plt.figure()
    # plt.scatter(results_p,results_yw)




    plt.plot(vp,yw,)

    # plt.scatter(results_p,results_yw,label=f'T={vT[j]}K')

    # plt.xlim(0,550)

    plt.ylim(0,5)

    plt.xlim(vp[0],vp[-1])


    plt.scatter(pexp,yexp)
    

    plt.xlabel(" P (bar)")

    plt.title(f"Temperatura = {T}K ; H20 2B ; TPD")

    plt.ylabel(" Water mol fraction ( % mol )")
    plt.legend()

# plt.scatter(pexp,yexp)

# plt.savefig(r'D:\CPA\CPA_mod\Exemplos\Plots\water_2B_304_tpd.png')

plt.show()



# plt.show()

# derivada dando 0; coisas estranhas acontecendo 


#%% teste fn_dew_point



# P = 50*1e5

# vp = np.linspace(1,700,100)

# yw = [P for ]
# T=300

# yw=[]
# for p in vp:

    # y = pure_water_dew_point(cpa_eos,T,p)
    # yw.append(y)


# plt.plot(vp,yw)

T= 294

P = 1*1e5

print(pure_water_dew_point(cpa_eos,T,P))



#%%Plot1 ( delta g contra yw á p=1bar t=300)

# 0.00201538

# vy = np.linspace(1e-6,1,100)

# delta_g= []

# T = 300
# P = 1*1e5


# # vT =[304,348,373.15]
# vT =[304]


# for t in vT:

#   delta_g_0 = []

#   for yw in vy:

#     z=np.array([ yw,1-yw] )
#     wguess = np.array([.999,.001])

#     dg,_ = TPD(T,P,y=z,xguess=wguess,phasey='vapor',phasex='liquid',eos=cpa_eos)

#     delta_g_0.append(dg)


#   delta_g.append(delta_g_0)


# for j in range(len(delta_g)):

#   plt.plot(vy,delta_g[j],label = f'T={vT[j]}K')

#   # plt.plot(vy,delta_g[j]
#   # plt.label(f'T={vT[j]}K')






# plt.hlines(0,0,0.8,linestyles='dashed')

# # plt.ylim(-1,1)

# # plt.xlim(0,.4)

# # plt.xlim(.6,.8)

# plt.legend()
# plt.xlabel('water fraction yw')
# plt.ylabel('deltaG/(eps*R*T)')
# plt.title(f'pressao ={P/1e5} bar')

# por que na temperatura nao tem 0, sendo que isso indicaria que nao existe
# yw no vapor tal que haja formação de fase liquida


# tudo bem, aumentar temperatura resulta no aumento da solubilidade da agua
# no vapor; logo, a energia  associada à formação de uma fase incipiente
# se torna cada vez menos espontanea, ao ponto de que em uma faixa de fração de agua
#  nao há mais espontaneidade nesse processo

#%%
data = [ 25.33, 2406, 50.66, 1696, 101.33, 3853, 202.65, 4445, 405.3, 5036, 506.63, 5076, 531.96, 5023, 557.29, 5050]

# # Reshape the data into a 2-column array
a1 = np.reshape(data, (len(data) // 2, 2))

# # Convert the second column (water fraction) from ppm to percentage
a1[:, 1] = a1[:, 1] * 1e-6*100  

l1 = np.array([	11.07	,88020 * 1e-6 * 100])
l2 = np.array([24.26,	46910* 1e-6 * 100])
l3 = np.array([70.88,	21020* 1e-6 * 100])
l4 = np.array([102.35	,17470* 1e-6 * 100])
l5= np.array([133.52,	16750* 1e-6 * 100])
l6 = np.array([170.7	,17120* 1e-6 * 100])

vT = np.array([373.15])
              

vp = np.linspace(1,750,100)

