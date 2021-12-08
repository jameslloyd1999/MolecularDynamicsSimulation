import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def rp(ri,rj,N,L):
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
                xij = rp(x[i],xj,N,L)
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






'''initial code'''

N=2 #number 0f dimensions
n=20 #number of particles
T=100000 #number of times
L=10.0 #Length of volume
temp = 1.0 #temperature
periodic=False

R = 0.5

dt=0.01

r=np.zeros((T,n,N))
v=np.zeros((T,n,N))
a=np.zeros((T,n,N))

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

'''starting velocities'''
vMax=np.sqrt(N*n*temp)
v0=np.random.random(N)*2-1
v0=vMax * v0 / np.sqrt( np.dot(v0,v0) )
v[0][0]=v0

'''starting accelerations'''
for j in range(0,n):
    a[0][j]=acceleration(r[0],j,n,N,L,periodic)






'''animation'''

fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(6,6)

if periodic:
    aR=0
else:
    aR=R

ax = plt.axes(xlim=(-aR,L+aR), ylim=(-aR,L+aR))
particles=[]
for i in range(0,n):
    particles.append(plt.Circle(r, R, fc='lime'))

def init():
    for i in range(0,n):
        particles[i].center = r[0][i]
        ax.add_patch(particles[i])
    return particles

def animate(i):
    global r
    global v
    global a
    r[i+1],v[i+1],a[i+1]=timeStep(r[i],v[i],a[i],N,n,L,dt,periodic)
    
    for j in range(0,n):
        x=r[i+1][j][0]
        y=r[i+1][j][1]
        if periodic:
            x=x%L
            y=y%L
        particles[j].center = (x,y)
    
    return particles

anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=T, 
                               interval=0.1,
                               blit=True)

plt.show()
