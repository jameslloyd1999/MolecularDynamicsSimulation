import numpy as np
import matplotlib.pyplot as plt

def rp(ri,rj,N,L,p):
    if p == False:
        return rj - ri
    rij=np.zeros(N)
    for i in range(0,N):
        xij = rj[i]%L - ri[i]%L
        xij = xij - L*np.round(xij/L)
        rij[i] = xij
    return rij

def acceleration(x,j,n,N,L,p):
    a=np.zeros(N)
    xj=x[j]
    for i in range(0,n):
        if i != j:
            if p:
                xij = rp(x[i],xj,N,L,p)
            else:
                xij = xj - x[i]
            r= np.sqrt( np.dot(xij,xij) )
            s=(1/r)**14-(1/r)**8
            a=a+s*xij
    return 12*a

def boundary(x,v,N,L):
    for i in range(0,N):
        if x[i] <= 0 or x[i] >= L:
            v[i]=-v[i]
    return v

def timeStep(r,v,a,N,n,L,dt,p):
    rNext=np.zeros((n,N))
    vNext=np.zeros((n,N))
    aNext=np.zeros((n,N))
    for i in range(0,n):
        rNext[i]=r[i]+v[i]*dt+0.5*a[i]*dt**2
    for i in range(0,n):
        aNext[i]=acceleration(rNext.copy(),i,n,N,L,p)
    for i in range(0,n):
        vNext[i]=v[i]+0.5*(aNext[i]+a[i])*dt
        if not p:
            vNext[i]=boundary(rNext[i].copy(),vNext[i].copy(),N,L)
    return rNext,vNext,aNext

def simulate(N,n,T,L,dt,temp,p,R):
    r=np.zeros((T,n,N))
    v=np.zeros((T,n,N))
    a=np.zeros((T,n,N))
    t=np.arange(0,T*dt,dt)
    
    '''starting positions'''
    r0=np.zeros((n,N))
    for i in range(0,n):
        collision=True
        while(collision==True):
            #randCord=np.random.random(N)*(L-2*R)+R
            randCord=np.random.random(N)*L
            collision=False
            for j in range(0,i):
                r1=randCord
                r2=r0[j]
                distVect = rp(r1,r2,N,L,p)
                distSquared=np.dot(distVect,distVect)
                dist=np.sqrt(distSquared)
                if dist<=2*R:
                    collision=True
        r0[i]=randCord
    '''dim=int(L)
    for i in range(0,dim):
        for j in range(0,dim):
            r0[dim*i+j]=[0.5+i,0.5+j]'''
    r[0]=r0

    '''starting velocities'''
    vMax=np.sqrt(N*n*temp)
    v0=np.random.random(N)*2-1
    v0=vMax * v0 / np.sqrt( np.dot(v0,v0) )
    v[0][0]=v0
    
    '''starting accelerations'''
    for j in range(0,n):
        a[0][j]=acceleration(r[0],j,n,N,L,p)
    
    for i in range(0,T-1):
        r[i+1],v[i+1],a[i+1] = timeStep(r[i],v[i],a[i],N,n,L,dt,p)
        if 100.0*i/T==int(100.0*i/T):
            print(str(100.0*i/T)+"%")
    return r,v,a,t

def path(r,L,n):
    x=[]
    y=[]
    for j in range(0,n):
        for i in range(0,len(r)):
            x.append(r[i][j][0])
            y.append(r[i][j][1])
        plt.plot(x,y)
        x=[]
        y=[]
    plt.show()

def distance(r,t,N,L,p):
    x=np.zeros(len(t))
    for i in range(0,len(t)):
        x01 = rp(r[i][0],r[i][1],N,L,p)
        d = np.sqrt( np.dot(x01,x01) )
        x[i]=d
    plt.plot(t,x)
    plt.title('Particles distance apart')
    plt.xlabel('$\\tau$')
    plt.ylabel('R')
    plt.show()

def energy(r,v,t,n,L,p):
    kinetic = np.zeros(len(t))
    potential = np.zeros(len(t))
    
    for i in range(0,len(t)):
        kSum = 0
        for j in range(0,n):
            kSum = kSum + 0.5 * np.dot(v[i][j],v[i][j])
        kinetic[i] = kSum
        
    for i in range(0,len(t)):
        pSum = 0
        for j in range(0,n):
            for k in range(j+1,n):
                rjk = rp(r[i][j],r[i][k],N,L,p)
                rSquared = np.dot(rjk,rjk)
                pSum = pSum + (1/rSquared)**6 - 2 * (1/rSquared)**3
        potential[i] = pSum
    
    #plt.plot(t,kinetic,label='kinetic')
    #plt.plot(t,potential,label='potential')
    plt.plot(t,potential+kinetic,label='total')
    plt.title('Total energy of the system')
    plt.xlabel('$\\tau$')
    plt.ylabel('$E/\epsilon$')
    #plt.legend()
    plt.show()

def pairCorrelation(x,N,n,L,T,p,temp):
    fs=10
    dr=L/200
    factor=2/(T*n*(n-1)*np.pi*dr)
    V=4/3 *np.pi * L**3/8
    A=np.pi*L**2/4
    if N==2:
        factor = factor*A/2
    if N==3:
        factor = factor*V/4
    r=np.arange(dr/2,L/2+dr/2,dr)
    g=np.zeros(r.size)
    for i in range(0,T):
        for j in range(0,n):
            for k in range(0,n):
                if j != k:
                    xj = x[i][j]
                    xk = x[i][k]
                    dist = rp(xj,xk,N,L,p)
                    distMod = np.sqrt( np.dot(dist,dist) )
                    if distMod < L/2:
                        ri = int(distMod/dr)
                        #g[ri] = g[ri] + 1/T *2*V/(n*(n-1)*4*np.pi*r[ri]**2*dr)
                        #g[ri] = g[ri] + 1/T *2*A/(n*(n-1)*2*np.pi*r[ri]*dr)
                        g[ri] = g[ri] + 1/r[ri]**(N-1)
        if 100.0*i/T==int(100.0*i/T):
            print(str(100.0*i/T)+"%")
    g=g*factor
    #return g,r
    plt.bar(r,g,dr)
    plt.xlabel('R',fontsize=fs)
    plt.ylabel('g(R)',fontsize=fs)
    plt.title('Pair correlation function for K='+str(temp)+' and '+str(n)+' particles',fontsize=fs)
    plt.show()
    
def bigPairCorrelation(xs,ns,Ks):
    fs=30
    fig, ax = plt.subplots(2,2)
    gs=[]
    rs=[]
    dr=L/200
    for i in range(0,4):
        g,r= pairCorrelation(xs[i],N,ns[i],L,T,p,Ks[i])
        gs.append(g)
        rs.append(r)
    for j in range(0,2):
        for k in range(0,2):
            i=2*j+k
            ax[j,k].bar(rs[i],gs[i],dr)
            ax[j,k].set_title('n='+str(ns[i])+' K='+str(Ks[i]),fontsize=fs)
            ax[j,k].set_xlabel('R',fontsize=fs,x=1)
            ax[j,k].set_ylabel('g(R)',fontsize=fs)
    fig.suptitle('Pair correlations',fontsize=fs)
    plt.show()
            

def velocityDistribution(v,n,k):
    vMod=[]
    for i in range(0,len(v)):
        for j in range(0,len(v[0])):
            vSquared=np.dot(v[i][j],v[i][j])
            vMod.append(np.sqrt(vSquared))
        print(str(100.0*i/len(v))+'%')
    '''plt.hist(vMod,100,density=True)
    plt.xlabel('V')
    plt.ylabel('$\\rho(V)$')
    plt.xlim(0,5)
    plt.ylim(0,0.8)
    plt.title('Velocity distribution for '+str(n)+' particles in a box of length '+str(int(L))+'$r_0$ at K='+str(k))
    '''
    x = np.linspace(0,max(vMod),1000)
    y = np.sqrt(2/(np.pi*k**3)) * x**2 * np.exp(-x**2/(2*k))
    #plt.plot(x,y)
    
    #plt.show()
    return vMod,x,y

def bigVelocityDistributionPlot(vs,ns,Ks):
    fs=30
    fig, ax = plt.subplots(2,2)
    vDist=[]
    x=[]
    y=[]
    for i in range(0,4):
        vMod,xi,yi = velocityDistribution(vs[i], ns[i], Ks[i])
        vDist.append(vMod)
        x.append(xi)
        y.append(yi)
    for j in range(0,2):
        for k in range(0,2):
            i=2*j+k
            ax[j,k].set_title(str(ns[i])+' particles, K='+str(Ks[i]),fontsize=fs)
            ax[j,k].set_xlim(0,5)
            ax[j,k].set_ylim(0,0.8)
            ax[j,k].set_xlabel('V',fontsize=fs,x=1)
            ax[j,k].set_ylabel('$\\rho(V)$',fontsize=fs)
            ax[j,k].hist(vDist[i],100,density=True)
            ax[j,k].plot(x[i],y[i])
    plt.show()
    
def plotHeatmap(r,L,N):
    x=[]
    y=[]
    for i in range(0,len(r)):
        for j in range(0,len(r[0])):
            xi=r[i][j][0] %L
            yi=r[i][j][1] %L
            x.append(xi)
            y.append(yi)
    
    return x,y
    


def bigDensityMap(r,L,N,p):
    fs=30
    fig, ax = plt.subplots(1,2)
    x,y = plotHeatmap(r,L,N)
    if p==True:
        s=' and periodic boundary conditions'
    else:
        s=' and hard boundary conditions'
    fig.suptitle(str(N)+'D particle distribution with '+str(n)+' particles'+s,fontsize=fs)
    
    ax[0].hist2d(x,y,bins=100,range=[[0,L],[0,L]])
    ax[0].set_xlabel('X',fontsize=fs)
    ax[0].set_ylabel('Y',fontsize=fs)
    ax[0].set_title('2D slice',fontsize=fs)
    #fig.colorbar(,ax=ax[0])
    
    ax[1].hist(x,100,density=True,range=[0,L])
    ax[1].set_xlabel('X',fontsize=fs)
    ax[1].set_ylabel('$\\rho(X)$',fontsize=fs)
    ax[1].set_title('1D slice',fontsize=fs)
    
    plt.show()

def diffusion(r,t,n,tTherm):
    dRMS=np.zeros(len(t))
    for i in range(len(t)):
        dSTotal=0
        for j in range(n):
            displacement=r[i][j]-r[0][j]
            dS = np.dot(displacement,displacement)
            dSTotal = dSTotal + dS
        dRMS[i]=np.sqrt(dSTotal/n)

    #t[0]=t[1]
    #dRMS[0]=dRMS[1]
    for i in range(0,tTherm):
        t[i]=t[tTherm]
        dRMS[i]=dRMS[tTherm]

    return t,dRMS

def diffPlot(N,n,T,L,R,tTherm,sims):
    fs=30
    times=[]
    dRMSs=[]
    
    for i in range(0,sims):
        r,v,a,t = simulate(N,n,T,L,dt,temp,p,R)
        ti,drms = diffusion(r,t,n,tTherm)
        times.append(ti)
        dRMSs.append(drms)
    
    fig, ax = plt.subplots(1,2)
    
    halfx=[np.log(times[0][0]), np.log(times[0][-1])]
    halfy=[np.log(dRMSs[0][0]), np.log(dRMSs[0][0])+0.5*(halfx[1] - halfx[0])]
    ax[1].plot(halfx,halfy,linestyle=':',color='grey')
    
    onex=[halfx[0], halfx[0] + 0.5*(halfx[1]-halfx[0])]
    oney=[halfy[0], halfy[1]]
    ax[1].plot(onex,oney,linestyle=':',color='grey')
    
    for i in range(0,sims):
        ax[0].plot(times[i],dRMSs[i])
        ax[1].plot(np.log(times[i]),np.log(dRMSs[i]))
    
    
    
    fig.suptitle('$D_{RMS}$ vs $\\tau$',fontsize=fs)
    ax[0].set_xlabel('$\\tau$',fontsize=fs)
    ax[0].set_ylabel('$D_{RMS}$',fontsize=fs)
    ax[1].set_xlabel('$ln(\\tau)$',fontsize=fs)
    ax[1].set_ylabel('$ln(D_{RMS})$',fontsize=fs)
    
    plt.show()
            
N=2 #number 0f dimensions
n=10 #number of particles
T=1000 #number of times
L=5.0 #Length of volume
temp = 1.0 #temperature
p = False #periodic
R=0.5

savestring = ' N'+str(N)+' n'+str(n)+' L'+str(L)+' K'+str(temp)+' p'+str(p)+' T'+str(T)


dt=0.01

'''rs=[]
ns=[50,50,50,10]
Ks=[0.2,2.0,20.0,1.0]

for i in range(0,4):
    r,v,a,t=simulate(N,ns[i],T,L,dt,Ks[i],p,R)
    rs.append(r)'''

r,v,a,t = simulate(N,n,T,L,dt,temp,p,R)

'''save data'''
#np.savetxt('r'+savestring+'.txt', r.reshape(T*n,N))
#np.savetxt('v'+savestring+'.txt', v.reshape(T*n,N))
#np.savetxt('t'+savestring+'.txt',t)

'''read data'''
#r = np.loadtxt('r'+savestring+'.txt').reshape(T,n,N)
#v = np.loadtxt('v'+savestring+'.txt').reshape(T,n,N)
#t = np.loadtxt('t'+savestring+'.txt')

path(r,L,n)
#distance(r,t,p,L,p)
#energy(r,v,t,n,L,p)
#pairCorrelation(r,N,n,L,T,p,temp)
#velocityDistribution(v,n,L,temp)
#bigVelocityDistributionPlot(vs,ns,Ks)
#bigDensityMap(r,L,N,p)
#bigPairCorrelation(rs,ns,Ks)
#diffPlot(N,n,T,L,R,500,5)