
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
   
 
def espectro(Hs,fs,smax,Omax,Omin,fmax,precisf,precisO):
    
    fp = fs / 1.05
    
    x = np.arange(0.01,fmax,precisf)
    
    y = np.arange(-np.pi,np.pi,precisO)
    
    XX, YY = np.meshgrid(x,y)
    
    Z = z(Hs,x,y,fp,fs,smax,Omax,Omin)
        
    plt.contourf(XX,YY,Z)
    
    plt.xlabel("$Frecuencia$")
    
    plt.ylabel("$Angulos\ de\ incidencia$")
    
    plt.title("$Espectro\ de\ oleaje$")
    
    plt.colorbar()
    
    plt.contour(XX,YY,Z)
    
    plt.grid()
    
    return Z
    

def z(Hs,x,y,fp,fs,smax,Omax,Omin):
    
    Z = np.zeros((np.size(y),np.size(x)))
    
    for j in range(0,np.size(y)-1):
        
        for i in range(0,np.size(x)-1):
            
                s = ((x[i] / fp)**5) * smax * (x[i] <= fp) + ((x[i] / fp)**(-2.5)) * smax * (x[i] > fp )                
                              
                g = lambda aux,s: (np.cos(aux/2))**(2*s)
                                               
                G0 = quad(g,Omax,Omin,args=s)

                F = (0.257 * (Hs**2) * (fs**4) * x[i]**(-5) * np.exp(-1.03*( (fs**(-1))*x[i])**(-4)))

                G = (G0[-1]**(-1))*np.cos(y[j]/2)**(2*s)
                
                Z[j,i] = F * G
        
    return Z