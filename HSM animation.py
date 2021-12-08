import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

'''functions'''

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
            c=plt.Circle((r[i][j][0],r[i][j][1]),R,color='grey')
            ax[row,colum].add_artist(c)
            if np.dot(v[i][j],v[i][j]) != 0:
                ax[row,colum].arrow(r[i][j][0],r[i][j][1],v[i][j][0],v[i][j][1],width=0.01,color='black')
            ax[row,colum].set_title('t='+str(round(t[i],4)))
            ax[row,colum].set_xlim(0,L)
            ax[row,colum].set_ylim(0,L)
            ax[row,colum].get_xaxis().set_visible(False)
            ax[row,colum].get_yaxis().set_visible(False)
    plt.show()

def interpolate(r,v,t,T,tStep):
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
    
    return time,pos,vel

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
    r[0]=r0
    print('initial positions set')
    
    '''starting velocities'''
    vMax=np.sqrt(N*n)
    v0=np.random.random(N)*2-1
    v0=vMax * v0 / np.sqrt( np.dot(v0,v0) )
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
        print(str(100.0*i/T)+'%')
    
    print("Pair collisions: " + str(c) + ", Wall collisions: " + str(w))
    rho=[0,0]
    rho[0] = n * np.pi * R**2 / L**2
    rho[1] = n * 4.0/3.0 *np.pi * R**3 / L**3
    print(str(rho[N-2]*100) + "% space occupied")

    return r,v,t







'''initial code'''

N=2 #number 0f dimensions
n=20 #number of particles
T=1000 #number of collisions
L=20.0 #Length of volume
R=1.0 #Radius of spheres

r,v,t = simulate(N,n,T,L,R)

tStep=0.01
time,pos,vel=interpolate(r,v,t,T,tStep)






'''animation'''

fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(6,6)

ax = plt.axes(xlim=(0,L), ylim=(0,L))
particles=[]
for i in range(0,n):
    particles.append(plt.Circle(r, R, fc='y'))

def init():
    for i in range(0,n):
        particles[i].center = r[0][i]
        ax.add_patch(particles[i])
    return particles

def animate(i):
    
    
    for j in range(0,n):
        x=pos[i][j][0]
        y=pos[i][j][1]
        particles[j].center = (x,y)
    
    return particles

anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=len(time), 
                               interval=0.1,
                               blit=True)

plt.show()
