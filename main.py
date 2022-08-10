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
J2000 = dat.date(2000, 1, 1)                                    # Reference epoch (1st July 2000, 12:00 GMT)
ellipsepoints = 50                                              # Orbit line resolution
anidpi = 100
anifps = 50
aniframes = 150
fig, ax = plt.subplots(figsize=(20.5, 10), subplot_kw={'projection':'3d'})
axlims = keplerel['a'].max()
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
        print('Invalid date format')
        exit()
    if TT > dat.date(3000, 12, 31):
        print('Date outside specified range')
        exit()
    era = input('AD/BC: ')
    start = (dat.date.today() - J2000).days / 36525             # Centuries between today and J2000
    if era == 'AD':
        end = (TT - J2000).days / 36525                         # Centuries between input final date and J2000
    elif era == 'BC':
        end = - ((TT - dat.date(1, 1, 1)) + (J2000 -dat.date(1, 1, 1))).days / 36525
    else:
        print('Must input "AD" or "BC"')
        exit()
    return start, end

def calcorbit(planet, t):
    # Propagation of Kepler orbital elements using Standish's linear fit
    # Units: centuries, astronomical units, radians
    a, e, w, I, W, M = np.array(planet[1::2]) + t * np.array(planet[2::2])
    if planet[0] in jovians.index:                              # Jovian planets correction terms
        jovian = jovians.loc[planet[0]]
        b, c, s, f = jovian['b'], jovian['c'], jovian['s'], jovian['f'] # Radians
        M += b * t ** 2 + c * np.cos(f * t) + s * np.sin(f * t)
    M = M % (2 * np.pi)                                         # Modulus M between -π and π rad
    if M > np.pi:
        M -= 2 * np.pi
    elif M < -np.pi:
        M += 2 * np.pi
    K = lambda E: E - e * np.sin(E) - M                         # Kepler's equation = 0
    dKdE = lambda E: 1 - e * np.cos(E)                          # Derivative
    E = scp.optimize.newton(K, M, fprime=dKdE, tol=1.75e-8)
    # Coords. in heliocentric orbital plane (2D)
    x = np.array([[a * (np.cos(E) - e)],
                  [a * np.sqrt(1 - e ** 2) * np.sin(E)]])
    fullrot = np.linspace(0, 2 * np.pi, ellipsepoints)
    ellipse = np.vstack((a * (np.cos(fullrot) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(fullrot)))
    # Coords. in J2000 ecliptic plane (3D), rotation matrix (Rz(W) * Rx(I) * Rz(w) * v)
    cw, sw = np.cos(w), np.sin(w)
    cI, sI = np.cos(I), np.sin(I)
    cW, sW = np.cos(W), np.sin(W)
    rot = np.array([[cw * cW - sw * sW * cI, -sw * cW - cw * sW * cI],
                    [cw * sW + sw * cW * cI, -sw * sW + cw * cW * cI],
                    [sw * sI, cw * sI]])
    xecl = rot @ x
    ellipsecl = rot @ ellipse
    return xecl, ellipsecl

def init():                                                     # Initialise figure
    ax.set_xlim3d(-axlims, axlims)
    ax.set_ylim3d(-axlims, axlims)
    ax.set_zlim3d(-axlims, axlims)
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
    for planet in keplerel.itertuples():                        # Get coordinates for planets and orbits
        xecl, ellipsecl = calcorbit(planet, t)
        x = np.hstack((x, xecl))
        el = np.hstack((el,ellipsecl))
    for idx, ln1 in enumerate(lns1):                            # Plot planets
        ln1.set_data([x[0, idx]], [x[1, idx]])
        ln1.set_3d_properties([x[2, idx]])
        ln1.set_label(keplerel.index[idx])
    for n, ln2 in enumerate(lns2):                              # Plot orbits
        idx = np.arange(ellipsepoints * n, ellipsepoints * (n + 1))
        ln2.set_data(el[0, idx], el[1, idx])
        ln2.set_3d_properties(el[2, idx])
    ax.legend(bbox_to_anchor=(1.5, .8))
    try:                                                        # Display date
        date = J2000 + dat.timedelta(days=t * 36525)
        era = 'AD'
    except OverflowError:
        date = dat.date(1, 1, 1) - dat.timedelta(days=t * 36525) - (J2000 - dat.date(1, 1, 1))
        era = 'BC'
    ax.set_title(f'Solar System, ECLIPJ2000 reference frame, {date} {era} (TT)', fontsize=20)
    return *lns1, *lns2

def main():
    start, end = getdate()
    frames = np.concatenate((start * np.ones(int(anifps / 2)),  # Linger on starting frame
                             np.linspace(start, end, aniframes - anifps),
                             end * np.ones(int(anifps / 2))))   # Linger on final frame
    ani = FuncAnimation(fig, animate, frames=frames, init_func=init, interval=50, blit=True)
    ani.save('solar_system.gif', writer='pillow', fps=anifps, dpi=anidpi)
    os.system('"solar_system.gif"')

if __name__=="__main__":
    main()