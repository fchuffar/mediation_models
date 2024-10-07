# loading useful libraries
#matplotlib widget

import numpy as np
#import matplotlib.pyplot as plt
import copy
import time

import sys



# compute contact probability for a given set of gaps and loops
def Pcij(posi, posj, csglfm, glm, llm):
    sij=posj-posi
    deb=csglfm[0]-glm[0]-llm[0]
    fin=csglfm[-1]-glm[-1]-llm[-1]    
    
    if posi<deb+glm[0]: #i is in a gap
        if posj<fin+glm[-1]: #j is in a gap
            if deb==fin:#i & j in the same gap
                return sij**(-1.5)
            else:
                return (deb+glm[0]-posi+posj-fin+1+np.sum(glm[1:-1]+1)-np.sum(llm[0:-1]==0))**(-1.5)
        else: #j is in a loop
            if (llm[-1]-1)==0:
                r2j=0
            else:
                r2j=(posj-(fin+glm[-1]))*(csglfm[-1]-1-posj)/(llm[-1]-1)
            return ((deb+glm[0]-posi+np.sum(glm[1:]+1)-np.sum(llm[0:-1]==0))+r2j)**(-1.5)
    else: #i is in a loop
        if posj<fin+glm[-1]: #j is in a gap
            if (llm[0]-1)==0:
                r2i=0
            else:
                r2i=(csglfm[0]-1-posi)*(posi-(deb+glm[0]))/(llm[0]-1)
            return (r2i+np.sum(glm[1:-1]+1)+posj-fin+1-np.sum(llm[1:-1]==0))**(-1.5)
        else:#j is in a loop
            if deb==fin:#i & j in the same loop
                return (max(1, sij*(posi-posj+llm[0]-1)/(llm[0]-1)))**(-1.5)
            else:
                if (llm[0]-1)==0:
                    r2i=0
                else:
                    r2i=(csglfm[0]-1-posi)*(posi-(deb+glm[0]))/(llm[0]-1)
                if (llm[-1]-1)==0:
                    r2j=0
                else:
                    r2j=(posj-(fin+glm[-1]))*(csglfm[-1]-1-posj)/(llm[-1]-1)
                return (r2i+r2j+np.sum(glm[1:]+1)-np.sum(llm[1:-1]==0))**(-1.5)

# correct the set of loops depending the occupancy and position of a & b            
def correct(ca, cb, posa, posb, idmin, idmax, deb, fin, csglfo, glo, llo):
    gl=np.copy(glo)
    ll=np.copy(llo)
    sab=posb-posa
    if ca==1 and  posa>=deb+glo[idmin]: #one ctcf bound in a and a in a loop
        if cb==1 and  posb>=fin+glo[idmax]:#one ctcf bound in b and b in a loop
            if deb==fin: #same loop
                if (llo[idmin]-1)>0:
                    e=np.random.choice(3, p=[posa-(deb+glo[idmin]), sab, csglfo[idmin]-1-posb]/(llo[idmin]-1))
                    if e==0:
                        ll[idmin]=posa+1-deb-glo[idmin]
                        gl[idmin+1]=glo[idmin+1]+csglfo[idmin]-1-posa
                    elif e==1:
                        ll[idmin]=sab+1
                        gl[idmin]=posa-deb
                        gl[idmin+1]=glo[idmin+1]+csglfo[idmin]-1-posb
                    else:
                        gl[idmin]=posb-deb
                        ll[idmin]=csglfo[idmin]-posb
            else:
                if (llo[idmin]-1)>0:
                    e=np.random.choice(2, p=[posa-(deb+glo[idmin]), csglfo[idmin]-1-posa]/(llo[idmin]-1))
                    if e==0:
                        ll[idmin]=posa+1-deb-glo[idmin]
                        gl[idmin+1]=glo[idmin+1]+csglfo[idmin]-1-posa
                    else:
                        gl[idmin]=posa-deb
                        ll[idmin]=csglfo[idmin]-posa
                if (llo[idmax]-1)>0:
                    e=np.random.choice(2, p=[posb-(fin+glo[idmax]), csglfo[idmax]-1-posb]/(llo[idmax]-1))
                    if e==0:
                        ll[idmax]=posb+1-fin-glo[idmax]
                        gl[idmax+1]=glo[idmax+1]+csglfo[idmax]-1-posb
                    else:
                        gl[idmax]=gl[idmax]-glo[idmax]+posb-fin
                        ll[idmax]=csglfo[idmax]-posb
        else:
            if (llo[idmin]-1)>0:
                e=np.random.choice(2, p=[posa-(deb+glo[idmin]), csglfo[idmin]-1-posa]/(llo[idmin]-1))
                if e==0:
                    ll[idmin]=posa+1-deb-glo[idmin]
                    gl[idmin+1]=glo[idmin+1]+csglfo[idmin]-1-posa
                else:
                    gl[idmin]=posa-deb
                    ll[idmin]=csglfo[idmin]-posa
    elif cb==1 and  posb>=fin+glo[idmax]:#one ctcf bound in b and b in a loop
        if (llo[idmax]-1)>0:
            e=np.random.choice(2, p=[posb-(fin+glo[idmax]), csglfo[idmax]-1-posb]/(llo[idmax]-1))
            if e==0:
                ll[idmax]=posb+1-fin-glo[idmax]
                gl[idmax+1]=glo[idmax+1]+csglfo[idmax]-1-posb
            else:
                gl[idmax]=posb-fin
                ll[idmax]=csglfo[idmax]-posb
    return gl, ll

#compute contact probability for the 4 different possibilities of ctcf binding at a & b for one random iteration
def contact(sab=200, i=10, sij=180, sloop=100, gap=100):
    #sab: genomic distance between a & b, i: position of one locus relative to a, sij:genomic distance between i & j, a<b & i<j
    
    posa=10*(sloop+gap)+np.random.geometric(1/gap)+np.random.geometric(1/sloop) #position of a from the beginning of the larger domain
    posb=posa+sab #position of b from the beginning of the larger domain
    sdom=posb+10*(sloop+gap) #approximative size of the larger domain including a and b
    posi=posa+i #position of i from the begining of the larger domain
    posj=posi+sij #position of j from the begining of the larger domain

    Pij=np.zeros((2, 2)) #contact proba for the 4 possible states (0, 0) ; (1, 0) ; (0, 1) ; (1, 1)
        
    #generate random loops and gaps, with gap first
    s=0
    n=0
    glo=[]
    llo=[]
    while s<=sdom:
        g=np.random.geometric(1/gap)
        l=np.random.geometric(1/sloop)
        s=s+g+l
        n=n+1
        glo.append(g)
        llo.append(l)
    glo=np.array(glo)
    llo=np.array(llo)
    slo=glo+llo 
    csglfo=np.cumsum(slo)
    csglbo=np.cumsum(slo[::-1])[::-1] 
    sdom=csglbo[0]
    idxfo=(csglfo)>(posa)
    idxbo=(csglbo)>=(sdom-posb)
    idxn=np.arange(len(glo))[idxfo & idxbo]
    idmin=np.min(idxn)
    idmax=np.max(idxn)
    deb=csglfo[idmin]-glo[idmin]-llo[idmin]
    fin=csglfo[idmax]-glo[idmax]-llo[idmax]
    for ca in range(2):
        for cb in range(2):
            #correct for ctcf binding
            gl, ll=correct(ca, cb, posa, posb, idmin, idmax, deb, fin, csglfo, glo, llo)
            #keep only the gap/loop of interest (between i & j, included)            
            sl=gl+ll
            csglf=np.cumsum(sl)
            csglb=np.cumsum(sl[::-1])[::-1]
            idxf=(csglf)>(posi)
            idxb=(csglb)>=(sdom-posj)        
            csglfm=csglf[idxf & idxb]
            glm=gl[idxf & idxb]
            llm=ll[idxf & idxb]
            
            Pij[ca, cb]=Pcij(posi, posj, csglfm, glm, llm)
    return Pij



#compute contact probability for many trajectories
def meancontact(niter=10000, sab=200, i=10, sij=180, sloop=100, gap=100):
    Pij=np.zeros((2, 2))
    for it in range(niter):
        Pij=Pij+contact(sab=sab, i=i, sij=sij, sloop=sloop, gap=gap)
        
    return (Pij/niter)



#simulate one Gillespie trajectory for the full system
def simumedfull(alpha0a=1, alpha1a=0, alpha0b=1, alpha1b=0, beta=1, b0a=1, b1a=0, b0b=1, b1b=0, delta=1, gamma0=5, gamma1=5, gamma2=0, Th=1, i=1, lamb=1, f=50, d=0.1, Pij=np.array([[0, 0], [0, 0]]), timemax=100, seed=0): 
    #random initial state with seed
    np.random.seed(seed=seed) #initialize the seed
    curstate=np.random.randint(0, 2, size=6) #ma, mb, ca, cb, nmrna, p
    curstate[4]=np.random.poisson(gamma0/d)

    transfer=[[-1, 0, 0, 0, 0, 0], # demethylation of a
              [ 1, 0, 0, 0, 0, 0], # methylation of a
              [ 0, -1, 0, 0, 0, 0], # demethylation of b
              [ 0, 1, 0, 0, 0, 0], # methylation of b
              [ 0, 0, -1, 0, 0, 0], # unbinding of ctcf at a
              [ 0, 0, 1, 0, 0, 0], # binding of ctcf at a
              [ 0, 0, 0, -1, 0, 0], # unbinding of ctcf at b
              [ 0, 0, 0, 1, 0, 0], # binding of ctcf at b
              [ 0, 0, 0, 0, -1, 0], # degradation of one mRNA
              [ 0, 0, 0, 0, 1, 0], # production of one mRNA
              [ 0, 0, 0, 0, 0, -1], # contact turnover
              [ 0, 0, 0, 0, 0, 1]] # contact formation 
    
    t=0
    logs=np.zeros((1, 7))


    while (t<timemax):
        cso=np.copy(curstate)
        [ma, mb, ca, cb, nrna, pcontact]=np.copy(curstate)

        log = np.append(t, cso)
#         print(log)
        logs=np.vstack([logs, log])

        # print(ma, mb, ca, cb)
        event=[ma*beta, #demethylation of a
               (1-ma)*(alpha0a*(1-ca)+alpha1a*ca), #methylation of a
               mb*beta, #demeth of b
               (1-mb)*(alpha0b*(1-cb)+alpha1b*cb), #meth of b
               ca*delta, #unbinding of a
               (1-ca)*(b0a*(1-ma)+b1a*ma), #binding of a
               cb*delta, #unbind of b
               (1-cb)*(b0b*(1-mb)+b1b*mb), #bind of b
               d*nrna, #mRNA degradation
               (gamma0*(1-ca*np.heaviside(Th-np.abs(i), 1))+gamma1*ca*np.heaviside(Th-np.abs(i), 1))+gamma2/f*pcontact, #mRNA production
               lamb*pcontact, #contact turnover
               lamb*(1-pcontact)*f*Pij[ca, cb]/(1-f*Pij[ca, cb])] #contact formation
        
        #print(event)
        rtot=np.sum(event)
        t=t-np.log(np.random.rand())/rtot #time of the next reaction (exponential law)
        e=np.random.choice(12, p=event/rtot)
        
        curstate=curstate+transfer[e]
        #print(t, curstate)

    # if seed<10:
    #     file_path = "sim_res_seed"+str(seed)+"_"+str(i)+"_"+str(sij)+"_"+str(gamma1)+"_"+str(gamma2)+"_"+str(lamb)+"_"+str(f)+"_"+str(d)+".txt"
    # #     print(file_path)
    # #     print(np.shape(logs))
    #     # Stockage de la matrice dans le fichier texte
    #     with open(file_path, 'w') as file:
    #         for log in logs:
    #             log_str = ' '.join(map(str, log))  # Conversion des éléments de la ligne en chaînes de caractères et les joindre avec un espace
    #             file.write(log_str + '\n')  # Écriture de la ligne dans le fichier avec un saut de ligne après chaque ligne

    return cso



#simulate many Gillespie trajectories
def simutot(ncell=100, alpha0a=1, alpha1a=0, alpha0b=1, alpha1b=0, beta=1, b0a=1, b1a=0, b0b=1, b1b=0, delta=1, gamma0=0, gamma1=0, gamma2=5, Th=1, lamb=1, f=50, d=0.1, niter=10000, sab=200, i=10, sij=180, sloop=100, gap=100, timemax=100, seed=0):
    np.random.seed(seed=seed) #initialize the seed
    statecell=np.zeros((ncell, 6))
    tm0=time.time()
    Pij=meancontact(niter=niter, sab=sab, i=i, sij=sij, sloop=sloop, gap=gap)
    print(Pij)
    for it in range(ncell):
        statecell[it, :]=simumedfull(alpha0a=alpha0a, alpha1a=alpha1a, alpha0b=alpha0b, alpha1b=alpha1b, beta=beta, b0a=1, b1a=b1a, b0b=b0b, b1b=b1b, delta=delta, gamma0=gamma0, gamma1=gamma1, gamma2=gamma2, Th=Th, i=i, lamb=lamb, f=f, d=d, Pij=Pij, timemax=timemax, seed=it)
        if it % 100 ==0:
            print('it:', it, '/', ncell)
    print('duration:', time.time()-tm0)  
    return statecell




# GO!
#ncell=1000
#i=1
#sij=100
#gamma1 = 0
#gamma2 = 50
#d=1
#lamb=1.
#f=50

print(sys.argv[1])

exec(open(sys.argv[1]).read())

out = simutot(ncell=ncell, gamma0=0, gamma1=gamma1, gamma2=gamma2, lamb=lamb, f=f, seed=seed, i=i, sij=sij, d=d)
# Chemin du fichier texte où la matrice sera stockée
file_path = "sim_res_"+str(ncell)+"_"+str(i)+"_"+str(sij)+"_"+str(gamma1)+"_"+str(gamma2)+"_"+str(lamb)+"_"+str(f)+"_"+str(d)+"_rep"+str(seed)+".txt"
print(file_path)
# Stockage de la matrice dans le fichier texte
with open(file_path, 'w') as file:
    for row in out:
        row_str = ' '.join(map(str, row))  # Conversion des éléments de la ligne en chaînes de caractères et les joindre avec un espace
        file.write(row_str + '\n')  # Écriture de la ligne dans le fichier avec un saut de ligne après chaque ligne

