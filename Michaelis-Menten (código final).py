# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 17:15:08 2022

@author: carva
"""


from sympy.physics.secondquant_stochastic import B, Bd, BKet, BBra, KroneckerDelta, n

from sympy import *

import sympy

import numpy as np

import matplotlib.pyplot as plt

import time

from scipy.linalg import expm

from random import *

from collections import Counter

inicio = time.time() 


Ns = 1
Ne = 1
Ni = Ne

n = 100000

freq = []

def ADC(S,E,I,C1,P,C2,C3):
    
    if C1 >= 0 and C2 >= 0 and P >= 0 and E >= 0 and S >= 0 and C3 >= 0 and I >= 0:
        freq.append((S,E,I,C1,P,C2,C3))

for i in range(Ns+1):
    
    for j in range(n):
        
        E = randint(0,Ne)
        I = randint(0, Ni)
        S = randint(0, Ns)
        P = i
        
        C1 = Ne - E - Ni + I
        C3 = Ns - S - P - C1
        C2 = Ni - I - C3

        ADC(S,E,I,C1,P,C2,C3)

lista = list(Counter(freq).keys())
lista.sort(reverse = True)

stateK = []
stateB = []

m = len(lista)

for i in range(0,m):
    stateK.append(BKet(lista[i]))

for i in range(0,m):
    stateB.append(BBra(lista[i]))

print(stateK)
#Número de estados:
state_number = len(stateK)

#t = Symbol('t')


x = 2
y = 1
z = 1
k = 10
w = 5
u = 9
o = 4
p = 10
r = 5
q = 0




def media(stateK, H, j):
    x, y = 0, 0
    for i in range(len(stateK)):
        if stateK[i][j] != 0:
            x += stateK[i][j]*H[i][0]
            y += (stateK[i][j]**2)*H[i][0]
    return x, y-x**2





def Hamiltoniana(stateKet,stateBra):
    
    H1 = x*stateKet[0]*stateKet[1]*(KroneckerDelta(stateBra[0],stateKet[0]-1)*KroneckerDelta(stateBra[1],stateKet[1]-1)*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3]+1)*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H2 = y*stateKet[3]*(KroneckerDelta(stateBra[0],stateKet[0]+1)*KroneckerDelta(stateBra[1],stateKet[1]+1)*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3]-1)*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H3 = z*stateKet[3]*(KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1]+1)*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3]-1)*KroneckerDelta(stateBra[4],stateKet[4]+1)*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H4 = k*stateKet[1]*stateKet[0]*stateKet[2]*(KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1]-1)*KroneckerDelta(stateBra[2],stateKet[2]-1)*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5]+1)*KroneckerDelta(stateBra[6],stateKet[6]) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H5 = w*stateKet[5]*stateKet[0]*(KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1]+1)*KroneckerDelta(stateBra[2],stateKet[2]+1)*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5]-1)*KroneckerDelta(stateBra[6],stateKet[6]) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H6 = u*stateKet[3]*stateKet[2]*(KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2]-1)*KroneckerDelta(stateBra[3],stateKet[3]-1)*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]+1) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H7 = o*stateKet[6]*(KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2]+1)*KroneckerDelta(stateBra[3],stateKet[3]+1)*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]-1) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H8 = p*stateKet[5]*stateKet[0]*(KroneckerDelta(stateBra[0],stateKet[0]-1)*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5]-1)*KroneckerDelta(stateBra[6],stateKet[6]+1) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H9 = r*stateKet[6]*(KroneckerDelta(stateBra[0],stateKet[0]+1)*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5]+1)*KroneckerDelta(stateBra[6],stateKet[6]-1) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))
    
    H10 = q*stateKet[6]*(KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1]+1)*KroneckerDelta(stateBra[2],stateKet[2]+1)*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4]+1)*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]-1) - KroneckerDelta(stateBra[0],stateKet[0])*KroneckerDelta(stateBra[1],stateKet[1])*KroneckerDelta(stateBra[2],stateKet[2])*KroneckerDelta(stateBra[3],stateKet[3])*KroneckerDelta(stateBra[4],stateKet[4])*KroneckerDelta(stateBra[5],stateKet[5])*KroneckerDelta(stateBra[6],stateKet[6]))    
    
    return H1 + H2 + H3 + H4 + H5 + H6 + H7 + H8 + H9 + H10


#Criando a matriz H:
m = []

for j in range(state_number):
     linha = []    
     stateBra = stateB[j]     
     for i in range(state_number):
         stateKet  = stateK[i]
         H = -Hamiltoniana(stateKet,stateBra)
         linha.append(float(H))
     m.append(linha)

m = np.array(m)




S = []
E = []
I = []
C1 = []
p = []
C2 = []
C3 = []
temp = []
soma = []


S_var = []
E_var = []
I_var = []
C1_var = []
p_var = []
C2_var = []
C3_var = []


#O cálculo H=QJQ^-1:
t = 0

while t <= 12:
    
    P = expm(-m*t)
    
    P = P.tolist()

    S.append(re(media(stateK, P, 0)[0]))
    E.append(re(media(stateK, P, 1)[0]))
    I.append(re(media(stateK, P, 2)[0]))
    C1.append(re(media(stateK, P, 3)[0]))
    p.append(re(media(stateK, P, 4)[0]))
    C2.append(re(media(stateK, P, 5)[0]))
    C3.append(re(media(stateK, P, 6)[0]))
    
    S_var.append(re(media(stateK, P, 0)[1]))
    E_var.append(re(media(stateK, P, 1)[1]))
    I_var.append(re(media(stateK, P, 2)[1]))
    C1_var.append(re(media(stateK, P, 3)[1]))
    p_var.append(re(media(stateK, P, 4)[1]))
    C2_var.append(re(media(stateK, P, 5)[1]))
    C3_var.append(re(media(stateK, P, 6)[1]))
    
    temp.append(t)
    soma.append(np.sum(P, axis=0))
    
    
    #t += 0.5
    
    if t <= 0.1:
        t += 0.001
    
    else:
        t += 0.01
    
    



## GRÁFICOS

l = int(stateK[0][0])

## GRÁFICOS DAS MÉDIAS

fig, ax = plt.subplots(figsize=(10, 7), facecolor=(1, 1, 1))

plt.yticks(np.arange(0, 11, 1))
plt.ylim(0, l)
plt.xticks(np.arange(0,17,1))
plt.xlim(-0.1, t)

ax.plot(temp,E,'r-',linewidth=3,label= r'$\langle E \rangle$')
ax.plot(temp,S,'m--',linewidth=3,label= r'$\langle S \rangle$')
ax.plot(temp,C1,'b:',linewidth=3,label= r'$\langle C_1 \rangle$')
ax.plot(temp,p,'g--',linewidth=3,label= r'$\langle P \rangle$')
if k != 0:
    ax.plot(temp,C2,'c-',linewidth=3,label= r'$\langle C_2 \rangle$')
if u != 0 or k != 0:
    ax.plot(temp,I,'b-',linewidth=3,label= r'$\langle I \rangle$')
if u != 0:
    ax.plot(temp,C3,'y',linewidth=3,label= r'$\langle C_3 \rangle$')
#ax.plot(temp,soma,'k--',linewidth=3,label='SOMA')


ax.legend(loc=1, bbox_to_anchor=(.9,1),fontsize=14,frameon = False, 
          ncol=2, handletextpad=0.2, columnspacing = 1)
#ax.legend(loc='best',fontsize=14,frameon = False)
plt.tick_params(labelsize=15)
plt.rcParams['legend.fontsize'] = 18
plt.xlabel('t',fontsize=20)



#ax1 = plt.axes([0.2, 0.6, 0.25, 0.25])
ax1 = plt.axes([0.2, 0.6, 0.25, 0.25])
        
plt.yticks(np.arange(0, float(stateK[0][0]) + 1, 2))
plt.xlim(-0.002, 0.1)
plt.xticks(np.arange(0,0.12,0.02))

ax1.plot(temp,E,'r-',linewidth=3)
ax1.plot(temp,S,'m--',linewidth=3)
ax1.plot(temp,C1,'b:',linewidth=3)
ax1.plot(temp,p,'g--',linewidth=3)
if k != 0:
    ax1.plot(temp,C2,'c-',linewidth=3)
if u != 0 or k != 0:
    ax1.plot(temp,I,'b-',linewidth=3)
if u != 0:
    ax1.plot(temp,C3,'y',linewidth=3)
#ax.plot(temp,soma,'k--',linewidth=3,label='SOMA')


plt.tick_params(labelsize=15)
plt.rcParams['legend.fontsize'] = 18
#plt.xlabel('t',fontsize=20)


plt.show()
fig.savefig('Médias (varias) Mixed - I3 maior - (k1+=k3+=k5+=10, k1-=k3-=k5-=5, k4+=9, k4-=4, k2=1).pdf', bbox_inches='tight',format='pdf', dpi=400)








## GRÁFICOS DAS VARIÂNCIAS

fig, ax = plt.subplots(figsize=(10, 7), facecolor=(1, 1, 1))

plt.yticks(np.arange(0,11,1))
plt.ylim(0, 4)
plt.xticks(np.arange(0,17,1))
plt.xlim(-0.1, t)

ax.plot(temp,E_var,'r-',linewidth=3,label= r'$ Var(E) $')
ax.plot(temp,S_var,'m--',linewidth=3,label= r'$ Var(S) $')
ax.plot(temp,C1_var,'b:',linewidth=3,label= r'$ Var(C_1) $')
ax.plot(temp,p_var,'g--',linewidth=3,label= r'$ Var(P) $')
if k != 0:
    ax.plot(temp,C2_var,'c-',linewidth=3,label= r'$ Var(C_2) $')
if u != 0 or k != 0:
    ax.plot(temp,I_var,'b-',linewidth=3,label= r'$ Var(I) $')
if u != 0:
    ax.plot(temp,C3_var,'y',linewidth=3,label= r'$ Var(C_3) $')
#ax.plot(temp,soma,'k--',linewidth=3,label='SOMA')

ax.legend(loc=1, bbox_to_anchor=(.2,1),fontsize=14,frameon = False, 
          ncol=1, handletextpad=0.2, columnspacing = 1)
#ax.legend(loc='best',fontsize=14,frameon = False)
plt.tick_params(labelsize=15)
plt.rcParams['legend.fontsize'] = 18
plt.xlabel('t',fontsize=20)




ax1 = plt.axes([0.57, 0.55, 0.3, 0.3])     

plt.yticks(np.arange(0,11,0.3))
plt.ylim(0,1.2)
plt.xticks(np.arange(0,0.2,0.01))
plt.xlim(-0.002, 0.04)


ax1.plot(temp,E_var,'r-',linewidth=3)
ax1.plot(temp,S_var,'m--',linewidth=3)
ax1.plot(temp,C1_var,'b:',linewidth=3)
ax1.plot(temp,p_var,'g--',linewidth=3)
if k != 0:
    ax1.plot(temp,C2_var,'c-',linewidth=3)
if u != 0 or k != 0:
    ax1.plot(temp,I_var,'b-',linewidth=3)
if u != 0:
    ax1.plot(temp,C3_var,'y',linewidth=3)
#ax.plot(temp,soma,'k--',linewidth=3,label='SOMA')


plt.tick_params(labelsize=15)
plt.rcParams['legend.fontsize'] = 18
plt.xlabel('t',fontsize=20)



plt.show()
fig.savefig('Variancia (varias) Mixed - I3 maior - (k1+=k3+=k5+=10, k1-=k3-=k5-=5, k4+=9, k4-=4, k2=1).pdf', bbox_inches='tight',format='pdf', dpi=400)




fim = time.time()
#print("\n\nTempo de execução: ",str(datetime.timedelta(seconds=(fim-inicio))))

print('\n\nTempo de execução:', fim - inicio)
