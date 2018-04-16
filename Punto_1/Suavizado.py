import numpy as np
import imageio
import matplotlib.pylab as plt
def Fourier2D(image):
    Shape=np.shape(image)
    Transform=np.zeros(Shape,dtype=complex)
    u=range(Shape[0])
    v=range(Shape[1])
    Norm=1/(Shape[0]*Shape[1])
    for i in u:
        for j in v:
            T=0
            for l in u:
                for m in v:
                    C=image[l,m]
                    T+=C*np.exp(-1j*2*np.pi*((u[i]*u[l]/Shape[0])+(v[j]*v[m]/Shape[1])))
            Transform[i,j]=T*Norm
    return(Transform)
def InvFourier2D(image):
    Shape=np.shape(image)
    Transform=np.zeros(Shape,dtype=complex)
    u=range(Shape[0])
    v=range(Shape[1])
    for i in u:
        for j in v:
            T=0
            for l in u:
                for m in v:
                    C=image[l,m]
                    T+=C*np.exp(1j*2*np.pi*((u[i]*u[l]/Shape[0])+(v[j]*v[m]/Shape[1])))
            Transform[i,j]=T
    return(Transform)
def Gaussian2D(Sigma,X,Y):
    M=len(X)
    N=len(Y)
    Gauss=np.zeros([M,N])
    A=(1)/(2*np.pi*(Sigma**2))
    for i in range(M):
        for j in range(N):
            Gauss[i,j]=A*np.exp(-(((X[i]-(M/2))**2 + (Y[j]-(N/2))**2)/(2*Sigma**2)))
    return(Gauss)
def Conc(R,G,B):
    Shape=np.shape(R)
    Im=np.zeros([Shape[0],Shape[1],3])
    for i in range(Shape[0]):
        for j in range(Shape[1]):
            Im[i,j]=np.array([R[i,j],G[i,j],B[i,j]])
    return (Im)
Imagen=imageio.imread('vm.png')
CanalR=Imagen[:,:,0]
CanalG=Imagen[:,:,1]
CanalB=Imagen[:,:,2]
SH=Imagen.shape
X=np.array(range(SH[0]))
Y=np.array(range(SH[1]))
FTR=Fourier2D(CanalR)
blurR=FTR*Gaussian2D(np.sqrt(SH[1]*SH[0]),X,Y)
InvR=np.real(invFourier2D(blurR))
FTG=Fourier2D(CanalG)
blurG=FTG*Gaussian2D(np.sqrt(SH[1]*SH[0]),X,Y)
InvG=np.real(invFourier2D(blurG))
FTB=Fourier2D(CanalB)
blurB=FTB*Gaussian2D(np.sqrt(SH[1]*SH[0]),X,Y)
InvB=np.real(invFourier2D(blurB))
imageio.imwrite('Suave.png',Conc(InvR,InvG,InvB))

