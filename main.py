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
J2000 = dat.date(2000,1,1)                                      # Reference epoch (1st July 2000, 12:00 GMT)
ellipsepoints = 50                                              # Orbit line resolution
anidpi = 100
anifps = 50
aniframes = 150
fig, ax = plt.subplots(figsize=(20.5,10), subplot_kw={'projection':'3d'})
lns1 = []
lns2 = []
for planet in keplerel.index:
    lobj, = ax.plot([], [], [], 'o', label=planet)
    lns1.append(lobj)
    lobj, = ax.plot([], [], [], 'k', lw=1)
    lns2.append(lobj)

def getdate():
    print('Input a date between 3000 BC and 3000 AD')
    try:
        TT = dat.date.fromisoformat(input('Date yyyy-mm-dd: ')) # Terrestrial time, approx. JPL ephemeris time
    except ValueError:
        print('Invalid date')
        exit()
    if TT > dat.date(3000,12,31):
        print('Date outside specified range')
        exit()
    start = (dat.date.today() - J2000).days / 36525             # Centuries between today and J2000
    end = (TT - J2000).days / 36525                             # Centuries between input final date and J2000
    return start, end

def calcorbit(planet, t):
    # Propagation of Kepler orbital elements using Standish's linear fit, t in centuries
    a = planet['a'] + t * planet['da/dt']                       # Astronomical units
    e = planet['e'] + t * planet['de/dt']                       # Radians
    w = planet['w'] + t * planet['dw/dt']
    I = planet['I'] + t * planet['dI/dt']
    W = planet['W'] + t * planet['dW/dt']
    M = planet['M'] + t * planet['dM/dt']
    if planet.name in jovians.index:                            # Jovian planets correction terms
        jovian = jovians.loc[planet.name]
        b, c, s, f = jovian['b'], jovian['c'], jovian['s'], jovian['f'] # Radians
        M += b * t ** 2 + c * np.cos(f * t) + s * np.sin(f * t)
    M = M % (2 * np.pi)                                         # Modulo M between -180º and 180º
    if M > np.pi:
        M -= 2 * np.pi
    elif M < -np.pi:
        M += 2 * np.pi
    K = lambda E: E - e * np.sin(E) - M                         # Kepler's equation = 0
    dKdE = lambda E: 1 - e * np.cos(E)                          # Derivative
    E = scp.optimize.newton(K, M, fprime=dKdE, tol=1.75e-8)
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
        ln1.set_data([], [])
    for ln2 in lns2:
        ln2.set_data([], [])
    return *lns1, *lns2

def animate(t):
    x = np.array([[], [], []])
    el = np.array([[], [], []])
    for _, planet in keplerel.iterrows():
        xecl, ellipsecl = calcorbit(planet, t)
        x = np.hstack((x, xecl))
        el = np.hstack((el,ellipsecl))
    for idx, ln1 in enumerate(lns1):
        ln1.set_data([x[0, idx]], [x[1, idx]])
        ln1.set_3d_properties([x[2, idx]])
        ln1.set_label(keplerel.index[idx])
    for n, ln2 in enumerate(lns2):
        idx = np.arange(ellipsepoints * n, ellipsepoints * (n + 1))
        ln2.set_data(el[0, idx], el[1, idx])
        ln2.set_3d_properties(el[2, idx])
    ax.legend(bbox_to_anchor=(1.5, .7))
    ax.set_title(f'Solar System, ECLIPJ2000 reference frame, {J2000 + dat.timedelta(days=t*36525)} (TT)')
    return *lns1, *lns2

def main():
    start, end = getdate()
    ani = FuncAnimation(fig, animate, frames=np.linspace(start, end, aniframes), init_func=init, interval=50, blit=True)
    ani.save('solar_system.gif', writer='pillow', fps=anifps, dpi=anidpi)
    os.system('"solar_system.gif"')

if __name__=="__main__":
    main()