import numpy as np
import matplotlib.pyplot as plt

def vPair(xk,xl,vk,vl,R):
    dx=xk-xl
    e=dx/np.sqrt(np.dot(dx,dx))
    dv=vk-vl
    vk2=vk-e*np.dot(dv,e)
    vl2=vl+e*np.dot(dv,e)
    return vk2,vl2

def wallCollisions(r0,v0,t0,n,L,R,N):
    tWall=np.inf
    v2=np.zeros(N)
    particle=-1
    
    for k in range(0,n):
        for l in range(0,N):
            x=r0[k][l]
            vx=v0[k][l]
            if vx != 0:
                t=(L-R-x)/vx
                if t>1e-10 and t<tWall:
                    tWall=t
                    v2=v0[k].copy()
                    v2[l]=-v2[l]
                    particle=k
                    
                t=(R-x)/vx
                if t>1e-10 and t<tWall:
                    tWall=t
                    v2=v0[k].copy()
                    v2[l]=-v2[l]
                    particle=k
        
    return t0+tWall,v2,particle

def pairCollision(r0,v0,t0,n,L,R):
    tPair=np.inf
    vNext=[[0,0],[0,0]]
    pair=[-1,-1]
    for k in range(0,n):
        for l in range(k+1,n):
            xk=r0[k]
            xl=r0[l]
            vk=v0[k]
            vl=v0[l]
            dx=xk-xl
            dv=vk-vl
            dxdv=np.dot(dx,dv)
            dxdx=np.dot(dx,dx)
            dvdv=np.dot(dv,dv)
            gamma=dxdv**2-dvdv*(dxdx-4*R**2)
            if gamma>0 and dxdv<0:
                t=-(dxdv+np.sqrt(gamma))/(dvdv)
                if t>1e-10 and t<tPair:
                    tPair=t
                    vNext=vPair(xk+vk*t,xl+vl*t,vk,vl,R)
                    pair=[k,l]
    return t0+tPair,vNext,pair

def plotParticles(r,v,t,n,T,L,R):
    plotSize=int(np.sqrt(T))
    fig, ax = plt.subplots(plotSize,plotSize)
    
    for i in range(0,T):
        row=int(i/plotSize)
        colum=i%plotSize
        for j in range(0,n):
            c=plt.Circle((r[i][j][0],r[i][j][1]),R,color='y')
            ax[row,colum].add_artist(c)
            if np.dot(v[i][j],v[i][j]) != 0:
                ax[row,colum].arrow(r[i][j][0],r[i][j][1],v[i][j][0],v[i][j][1],width=0.05,color='black')
            ax[row,colum].set_title('$\\tau$='+str(round(t[i],4)))
            ax[row,colum].set_xlim(0,L)
            ax[row,colum].set_ylim(0,L)
            ax[row,colum].get_xaxis().set_visible(False)
            ax[row,colum].get_yaxis().set_visible(False)
    fig.set_size_inches(7,7)
    plt.show()
    
def plotHeatmap(r,L,N):
    x=[]
    y=[]
    for i in range(0,len(r)):
        for j in range(0,len(r[0])):
            xi=r[i][j][0]
            yi=r[i][j][1]
            x.append(xi)
            y.append(yi)
    
    '''plt.hist2d(x,y,bins=100,density=True,range=[[0,L],[0,L]])
    plt.title('Heatmap ' + str(N) + 'D')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar()
    plt.show()'''
    return x,y
    
def distributionMap(r,t,L,n):
    x=[]
    for i in range(0,len(t)):
        for j in range(0,n):
            xi=r[i][j][0]
            #yi=r[i][j][1]
            #if yi>0.45*L and yi<0.55:
            x.append(xi)
    '''plt.hist(x,100,density=True,range=[0,L])
    plt.title('Heatmap slice ' + str(N) + 'D')
    plt.xlabel('X')
    plt.ylabel('position density')
    plt.show()'''
    return x

def bigDensityMap(r,L,N):
    fs=30
    fig, ax = plt.subplots(1,2)
    x,y = plotHeatmap(r,L,N)
    fig.suptitle(str(N)+'D particle distribution with '+str(n)+' particles',fontsize=fs)
    
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
    
def velocityDistribution(v,n,L):
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
    plt.title('Time averaged velocity distribution for '+str(n)+' particles in a box of length '+str(int(L))+' $r_0$')
    '''
    x = np.linspace(0,max(vMod),1000)
    y = np.sqrt(2/np.pi) * x**2 * np.exp(-x**2/2)
    ''''plt.plot(x,y)'''
    return vMod,x,y
    #plt.show()
    
def bigVelocityDistributionPlot(vs,ns,Ls):
    fs=30
    fig, ax = plt.subplots(2,2)
    vDist=[]
    x=[]
    y=[]
    for i in range(0,4):
        vMod,xi,yi = velocityDistribution(vs[i], ns[i], Ls[i])
        vDist.append(vMod)
        x.append(xi)
        y.append(yi)
    for j in range(0,2):
        for k in range(0,2):
            i=2*j+k
            ax[j,k].set_title(str(ns[i])+' particles',fontsize=fs)
            ax[j,k].set_xlim(0,5)
            ax[j,k].set_ylim(0,0.8)
            ax[j,k].set_xlabel('V',fontsize=fs,x=1)
            ax[j,k].set_ylabel('$\\rho(V)$',fontsize=fs)
            ax[j,k].hist(vDist[i],100,density=True)
            ax[j,k].plot(x[i],y[i])
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
        r,v,t = simulate(N,n,T,L,R)
        #r,v,t = thermalise(r,v,t,tTherm)
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
    

def path(r,t,L):
    x=[]
    y=[]
    for i in range(0,t):
        x.append(r[i][0][0])
        y.append(r[i][0][1])
    plt.plot(x,y)
    plt.xlim(0,L)
    plt.ylim(0,L)
    plt.show()

def pairCorrelation(x,N,n,L,T):
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
    for i in range(0,len(x)):
        for j in range(0,n):
            for k in range(0,n):
                if j != k:
                    xj = x[i][j]
                    xk = x[i][k]
                    dist = xj-xk
                    distMod = np.sqrt( np.dot(dist,dist) )
                    if distMod < L/2:
                        ri = int(distMod/dr)
                        #g[ri] = g[ri] + 1/T *2*V/(n*(n-1)*4*np.pi*r[ri]**2*dr)
                        #g[ri] = g[ri] + 1/T *2*A/(n*(n-1)*2*np.pi*r[ri]*dr)
                        g[ri] = g[ri] + 1/r[ri]**(N-1)
        #if 100.0*i/T==int(100.0*i/T):
            print(str(100.0*i/T)+"%")
    g=g*factor
    #return g,r
    plt.bar(r,g,dr)
    plt.xlabel('R',fontsize=fs)
    plt.ylabel('g(R)',fontsize=fs)
    plt.title('Pair correlation function for '+str(n)+' particles',fontsize=fs)
    plt.show()

def energy(v,t,n,T):
    k=np.zeros((n,T))
    kTot=np.zeros(T)
    for i in range(0,n):
        for j in range(0,T):
            k[i][j] = 0.5*np.dot(v[j][i],v[j][i])
            kTot[j]=kTot[j]+k[i][j]
    
    plt.plot(t,kTot,label='Total')
    for i in range(0,n):
        plt.plot(t,k[i])
    plt.title('Energy of hard disk model for 10 particles')
    plt.xlabel('$\\tau$')
    plt.ylabel('$E/k_B T$')
    plt.legend()
    plt.show()
    
def interpolate(r,v,t,T):
    '''tDiff=np.zeros(T-1)
    for i in range(0,T-1):
        tDiff[i] = t[i+1] - t[i]
    tStep=min(tDiff)'''
    tStep = 0.1*t[-1]/T
    
    time=np.arange(0,t[-1]+tStep,tStep)
    pos=[]
    vel=[]
    
    for i in range(0,T-1):
        ti=t[i]
        while(ti<t[i+1]):
            p=r[i]+v[i]*(ti-t[i])
            pos.append(p)
            vel.append(v[i])
            ti=ti+tStep
        if 100.0*i/T==int(100.0*i/T):
            print(str(100.0*i/T)+'%')
    
    return time,pos,vel

def thermalise(r,v,t,tTherm):
    time=np.zeros(len(t)-tTherm)
    pos=[]
    vel=[]
    
    for i in range(0,len(t)-tTherm):
        time[i]=t[i+tTherm]-t[tTherm]
        pos.append(r[i+tTherm])
        vel.append(v[i+tTherm])
    return pos,vel,time
    
def simulate(N,n,T,L,R):
    r=np.zeros((T,n,N))
    v=np.zeros((T,n,N))
    t=np.zeros(T)
    '''starting positions'''
    r0=np.zeros((n,N))
    for i in range(0,n):
        collision=True
        while(collision==True):
            randCord=np.random.random(N)*(L-2*R)+R
            collision=False
            for j in range(0,i):
                r1=randCord
                r2=r0[j]
                distSquared=np.dot(r1-r2,r1-r2)
                dist=np.sqrt(distSquared)
                if dist<=2*R:
                    collision=True
        r0[i]=randCord
    #r0=[[2,3.3],[6,6]]
    r[0]=r0
    print('initial positions set')
    
    '''starting velocities'''
    vMax=np.sqrt(N*n)
    v0=np.random.random(N)*2-1
    v0=vMax * v0 / np.sqrt( np.dot(v0,v0) )
    #v0=[2,2]
    v[0][0]=v0
    
    c,w=0,0
    for i in range(0,T-1):
        r0=r[i]
        v0=v[i]
        t0=t[i]
        tWall,vWall,kWall=wallCollisions(r0.copy(),v0.copy(),t0,n,L,R,N)
        tPair,vkl,kl=pairCollision(r0.copy(),v0.copy(),t0,n,L,R)
        v[i+1]=v0
        if tWall<tPair:
            tNext=tWall
            v[i+1][kWall]=vWall
            w=w+1
        else:
            tNext=tPair
            k=kl[0]
            l=kl[1]
            vk=vkl[0]
            vl=vkl[1]
            v[i+1][k]=vk
            v[i+1][l]=vl
            c=c+1
        
        r1=r0+v0*(tNext-t0)
        r[i+1]=r1
        t[i+1]=tNext
        if 100.0*i/T==int(100.0*i/T):
            print(str(100.0*i/T)+"%")
    
    print("Pair collisions: " + str(c) + ", Wall collisions: " + str(w))
    rho=[0,0]
    rho[0] = n * np.pi * R**2 / L**2
    rho[1] = n * 4.0/3.0 *np.pi * R**3 / L**3
    print(str(rho[N-2]*100) + "% space occupied")

    return r,v,t



N=2 #number 0f dimensions
n=10 #number of particles
T=16 #number of collisions
L=20.0 #Length of volume
R=1.0 #Radius of spheres

#vs=[]
#ns=[5,10,20,30]
#Ls=[10,10,10,10]

#for i in range(0,4):
#    r,v,t=simulate(N,ns[i],T,Ls[i],R)
#    vs.append(v)

r,v,t = simulate(N,n,T,L,R)

#time,pos,vel=interpolate(r,v,t,T)

plotParticles(r,v,t,n,T,L,R)
#plotHeatmap(pos,L,N)
#distributionMap(pos,time,L,n)
#velocityDistribution(v,n,L)
#animate(pos,vel,time,tStep,n,L,R)
#diffusion(r,t,n,200)
#diffPlot(N,n,T,L,R,100,1)
#path(r,T,L)
#bigVelocityDistributionPlot(vs, ns, Ls)
#bigDensityMap(pos,L,N)
#pairCorrelation(pos,N,n,L,len(pos))
#energy(v,t,n,T)