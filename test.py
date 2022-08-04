import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots(subplot_kw={'projection':'3d'})
xdata1, ydata1, zdata1 = [], [], []
xdata2, ydata2, zdata2 = [], [], []
ln1, = ax.plot([], [], [], 'o')
ln2, = ax.plot([], [], [], 'k', lw=1)

def plotorbit(E):
    a = 2; e = .8; w = np.pi/6; I = np.pi/3; W = np.pi/4
    x = np.array([[a * (np.cos(E) - e)],
                  [a * np.sqrt(1 - e ** 2) * np.sin(E)]])
    fullrot = np.linspace(0, 2 * np.pi, 100)
    ellipse = np.vstack((a * (np.cos(fullrot) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(fullrot)))
    cw = np.cos(w); sw = np.sin(w)
    cI = np.cos(I); sI = np.sin(I)
    cW = np.cos(W); sW = np.sin(W)
    rot = np.array([[cw * cW - sw * sW * cI, -sw * cW - cw * sW * cI],
                    [cw * sW + sw * cW * cI, -sw * sW + cw * cW * cI],
                    [sw * sI, cw * sI]])
    xecl = rot @ x
    ellipsecl = rot @ ellipse
    return xecl, ellipsecl

def init():
    ax.set_xlim3d(-3, 3)
    ax.set_ylim3d(-3, 3)
    ax.set_zlim3d(-3, 3)
    return ln1, ln2

def animate(frame):
    x, el = plotorbit(frame)
    ln1.set_data(x[0], x[1])
    ln1.set_3d_properties(x[2])
    ln2.set_data(el[0], el[1])
    ln2.set_3d_properties(el[2])
    return ln1, ln2

ani = FuncAnimation(fig, animate, frames=np.linspace(0, 2*np.pi, 128), init_func=init, blit=True)
plt.show()