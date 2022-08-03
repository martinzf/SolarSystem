import numpy as np
import matplotlib.pyplot as plt
a=20
e=0.9
fullrot = np.linspace(0, 2 * np.pi, 100)
ellipse = np.vstack((a * (np.cos(fullrot) - e),
                     a * np.sqrt(1 - e ** 2) * np.sin(fullrot),
                     np.zeros(100)))
w=np.pi/3
I=w
W=w
matrix=lambda x:np.array([[np.cos(x),-np.sin(x),0],[np.sin(x),np.cos(x),0]])
roll=np.vstack((matrix(w),np.array([[0,0,1]])))
yaw=np.vstack((matrix(W),np.array([[0,0,1]])))
matrix=lambda x:np.array([[0,np.cos(x),-np.sin(x)],[0,np.sin(x),np.cos(x)]])
pitch=np.vstack((np.array([[1,0,0]]),matrix(I)))
#ellipse=roll@ellipse
print(np.deg2rad(1e-6))