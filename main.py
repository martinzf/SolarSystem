import numpy as np
import scipy as scp
import pandas as pd
import datetime as dat
import os
import numba
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use(['science', 'notebook', 'grid'])

# E.M. Standish (1992)'s data for simplified orbit propagation
keplerel = pd.read_csv('keplerel.csv', index_col=0)             # Planets' Kepler elements
jovians = pd.read_csv('jovians.csv', index_col=0)               # Jovian planets' correction terms
ellipsepoints = 100
anidpi = 100
anifps = 30
fig, ax = plt.subplots(figsize=(20.5,10), subplot_kw={'projection':'3d'})
lns1 = []
lns2 = []
for planet, _ in keplerel.iterrows():
    lobj, = ax.plot([],[],[],'o',label=planet)
    lns1.append(lobj)
    lobj, = ax.plot([], [], [], 'k', lw=1)
    lns2.append(lobj)

def getdate():
    J2000 = dat.date(2000,1,1)                                  # Reference epoch (1st July 2000, 12:00 GMT)
    print('Input a date between 3000 BC and 3000 AD')
    try:
        TT = dat.date.fromisoformat(input('Date yyyy-mm-dd: ')) # Terrestrial time, approx. JPL ephemeris time
    except ValueError:
        print('Invalid date')
        exit()
    if TT > dat.date(3000,12,31):
        print('Date outside specified range')
        exit()
    start = (dat.date.today() - J2000).days                     # Days between today and J2000
    end = (TT - J2000).days                                     # Days between input final date and J2000
    return [start, end]

def keplerelements(planet, dt):
    # Propagation of Kepler orbital elements using Standish's linear fit, dt in centuries
    keplerel.loc[planet.name, 'a'] = planet['a'] + dt * planet['da/dt']
    e = planet['e'] + dt * planet['de/dt']
    keplerel.loc[planet.name, 'e'] = e
    keplerel.loc[planet.name, 'w'] = planet['w'] + dt * planet['dw/dt']
    keplerel.loc[planet.name, 'I'] = planet['I'] + dt * planet['dI/dt']
    keplerel.loc[planet.name, 'W'] = planet['W'] + dt * planet['dW/dt']
    M = planet['M'] + dt * planet['dM/dt']
    if planet.name in jovians.index:                            # Jovian planets correction terms
        jovian = jovians.loc[planet.name]
        b = jovian['b']; c = jovian['c']; s = jovian['s']; f = jovian['f']
        M = M + b * dt ** 2 + c * np.cos(f * dt) + s * np.sin(f * dt)
    M = M % (2 * np.pi)                                         # Modulo M between -180º and 180º
    if M > np.pi:
        M = M - 2 * np.pi
    keplerel.loc[planet.name, 'M'] = M
    K = lambda E: E - e * np.sin(E) - M                         # Kepler's equation = 0
    dKdE = lambda E: 1 - e * np.cos(E)                          # Derivative
    E = scp.optimize.newton(K, planet['E'], fprime=dKdE, tol=1.75e-8)
    keplerel.loc[planet.name, 'E'] = E

def calcorbit(planet):
    a, e, w, I, W, E = planet['a'], planet['e'], planet['w'], planet['I'], planet['W'], planet['E']
    # Coords. in heliocentric orbital plane
    x = np.array([[a * (np.cos(E) - e)],
                  [a * np.sqrt(1 - e ** 2) * np.sin(E)]])
    fullrot = np.linspace(0, 2 * np.pi, ellipsepoints)
    ellipse = np.vstack((a * (np.cos(fullrot) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(fullrot)))
    # Coords. in J2000 ecliptic plane, rotation matrix (anticlock z (w), anticlock x (I), anticlock z (W))
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
    ax.set_xlim3d(-25, 25)
    ax.set_ylim3d(-25, 25)
    ax.set_zlim3d(-25, 25)
    ax.set_xlabel(r'$\hat{X}$: Vernal Equinox (AU)', labelpad=10)
    ax.set_ylabel(r'$\hat{Y}=\hat{Z}\times\hat{X}$ (AU)', labelpad=10)
    ax.set_zlabel(r'$\hat{Z}$: North Ecliptic Pole (AU)', labelpad=10)
    for ln1 in lns1:
        ln1.set_data([],[])
    for ln2 in lns2:
        ln2.set_data([],[])
    return *lns1, *lns2

def animate(dt):
    x = np.array([[],[],[]])
    el = np.array([[],[],[]])
    for _, planet in keplerel.iterrows():
        keplerelements(planet, dt)
        xecl, ellipsecl = calcorbit(planet)
        x = np.hstack((x, xecl))
        el = np.hstack((el,ellipsecl))
    for idx, ln1 in enumerate(lns1):
        ln1.set_data([x[0, idx]], [x[1, idx]])
        ln1.set_3d_properties([x[2, idx]])
    for n, ln2 in enumerate(lns2):
        idx = np.arange(ellipsepoints * n, ellipsepoints * (n + 1))
        ln2.set_data(el[0, idx], el[1, idx])
        ln2.set_3d_properties(el[2, idx])
    return *lns1, *lns2

def main():
    '''start, end = getdate()
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    dt = start
    datetime = dat.datetime.combine(dat.date.today(), dat.time(12))'''
    ani = FuncAnimation(fig, animate, frames=np.ones(10), init_func=init, blit=True)
    ani.save('test.gif', writer='pillow', fps=anifps, dpi=anidpi)
    os.system('"test.gif"')

if __name__=="__main__":
    main()