import numpy as np
import scipy as scp
import pandas as pd
import datetime as dt
import os
import matplotlib.pyplot as plt
plt.style.use(['science', 'notebook', 'grid'])

keplerel = pd.read_csv('keplerel.csv', index_col=0)             # Planets' Kepler elements
jovians = pd.read_csv('jovians.csv', index_col=0)               # Jovian planets' correction terms

def plotorbit(planet, t, ax):
    # Propagation of Kepler orbital elements using Standish's linear fit
    a = planet['a'] + t * planet['da/dt']
    e = planet['e'] + t * planet['de/dt']
    I = np.deg2rad(planet['I'] + t * planet['dI/dt'])
    L = np.deg2rad(planet['L'] + t * planet['dL/dt'])
    long_peri = np.deg2rad(planet['long.peri.'] + t * planet['dlong.peri./dt'])
    W = np.deg2rad(planet['W'] + t * planet['dW/dt'])
    w = long_peri - W
    M = L - long_peri
    if planet.name in jovians.index:                            # Jovian planets correction terms
        jovian = np.deg2rad(jovians.loc[planet.name])
        b = jovian['b']; c = jovian['c']; s = jovian['s']; f = jovian['f']
        M = M + b * t ** 2 + c * np.cos(f * t) + s * np.sin(f * t)
    K = lambda E: E - e * np.sin(E) - M                         # Kepler's equation = 0
    dKdE = lambda E: 1 - e * np.cos(E)                          # Derivative
    E = scp.optimize.newton(K, 0, fprime=dKdE)
    # Coords. in heliocentric orbital plane
    x = np.array([a * (np.cos(E) - e),
                  a * np.sqrt(1 - e ** 2) * np.sin(E)])
    fullrot = np.linspace(0, 2 * np.pi, 100)
    ellipse = np.vstack((a * (np.cos(fullrot) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(fullrot)))
    # Coords. in J2000 ecliptic plane, rotation matrix (anticlock z (w), anticlock x (I), anticlock z (W))
    rot = np.array([[np.cos(w) * np.cos(W) - np.sin(w) * np.sin(W) * np.cos(I), -np.sin(w) * np.cos(W) - np.cos(w) * np.sin(W) * np.cos(I)],
                    [np.cos(w) * np.sin(W) + np.sin(w) * np.cos(W) * np.cos(I), -np.sin(w) * np.sin(W) + np.cos(w) * np.cos(W) * np.cos(I)],
                    [np.sin(w) * np.sin(I)                                    , np.cos(w) * np.sin(I)                                     ]])
    xecl = rot @ x
    ellipsecl = rot @ ellipse
    ax.plot(xecl[0], xecl[1], xecl[2], 'o', label=planet.name)
    ax.plot(ellipsecl[0], ellipsecl[1], ellipsecl[2], 'k', linewidth=1)

def main():
    J2000 = dt.date(2000, 1, 1)                                # Reference epoch (1st July 2000, ~12:00 GMT)
    print('Input a date between 3000 BC and 3000 AD')
    try:
        TT = dt.date.fromisoformat(input('Date yyyy-mm-dd: ')) # Terrestrial time, approx JPL ephemeris time
    except ValueError:
        print('Invalid date format')
        return
    t = (TT-J2000).days / 36525                                # Centuries between J2000.0 and input final date
    # E.M. Standish (1992)'s data for simplified orbit propagation
    ax = plt.axes(projection='3d')
    for _, planet in keplerel.iterrows():
        plotorbit(planet, t, ax)
    ax.set_xlim3d(-25,25)
    ax.set_ylim3d(-25,25)
    ax.set_zlim3d(-25,25)
    ax.set_xlabel(r'$\hat{X}$: Vernal Equinox (AU)', labelpad=10)
    ax.set_ylabel(r'$\hat{Y}=\hat{Z}\times\hat{X}$ (AU)', labelpad=10)
    ax.set_zlabel(r'$\hat{Z}$: North Ecliptic Pole (AU)', labelpad=10)
    ax.set_title(f'Solar System in J2000 Ecliptic Plane, {TT}(TT)')
    plt.legend(bbox_to_anchor = [1.8, .9])
    plt.savefig('solar_syst.png', dpi=200, bbox_inches='tight')
    os.system('pycharm64 solar_syst.png')

if __name__=="__main__":
    main()