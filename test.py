import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation
import os

fig, ax = plt.subplots(figsize=(20.5,10), subplot_kw={'projection':'3d'})
lns1 = []
lns2 = []

planet1 = pd.Series([2, .8, np.pi/6, np.pi/3, np.pi/4])
planet2 = pd.Series([1, .6, np.pi/3, np.pi/5, np.pi/6])
planets = [planet1, planet2]
for planet in planets:
    lobj, = ax.plot([],[],[],'o',label=planet.name)
    lns1.append(lobj)
    lobj, = ax.plot([], [], [], 'k', lw=1)
    lns2.append(lobj)

def calcorbit(planet, E):
    a, e, w, I, W = planet[0], planet[1], planet[2], planet[3], planet[4]
    x = np.array([[a * (np.cos(E) - e)],
                  [a * np.sqrt(1 - e ** 2) * np.sin(E)]])
    fullrot = np.linspace(0, 2 * np.pi, 100)
    ellipse = np.vstack((a * (np.cos(fullrot) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(fullrot)))
    cw, sw = np.cos(w), np.sin(w)
    cI, sI = np.cos(I), np.sin(I)
    cW, sW = np.cos(W), np.sin(W)
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
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    for ln1 in lns1:
        ln1.set_data([],[])
    for ln2 in lns2:
        ln2.set_data([],[])
    return *lns1, *lns2

def animate(frame):
    x = np.array([[],[],[]])
    el = np.array([[],[],[]])
    for planet in planets:
        xecl, ellipsecl = calcorbit(planet, frame)
        x = np.hstack((x,xecl))
        el = np.hstack((el,ellipsecl))
    for idx, ln1 in enumerate(lns1):
        ln1.set_data([x[0,idx]], [x[1,idx]])
        ln1.set_3d_properties([x[2,idx]])
    for idx, ln2 in enumerate(lns2):
        ln2.set_data(el[0,100*idx:100*(idx+1)], el[1,100*idx:100*(idx+1)])
        ln2.set_3d_properties(el[2,100*idx:100*(idx+1)])
    return *lns1, *lns2

ani = FuncAnimation(fig, animate, frames=np.linspace(0, 2*np.pi, 128), init_func=init, blit=True)
ani.save('test.gif', writer='pillow', fps=30, dpi=100)
os.system('"test.gif"')
