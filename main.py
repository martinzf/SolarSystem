import numpy as np
import scipy as scp
import pandas as pd
import datetime as dt
import os
import matplotlib.pyplot as plt
from PIL import Image
plt.style.use(['science', 'notebook', 'grid'])

# E.M. Standish (1992)'s data for simplified orbit propagation
keplerel = pd.read_csv('keplerel.csv', index_col=0)             # Planets' Kepler elements
jovians = pd.read_csv('jovians.csv', index_col=0)               # Jovian planets' correction terms

def getdate():
    J2000 = dt.date(2000, 1, 1)                                 # Reference epoch (1st July 2000, ~12:00 GMT)
    print('Input a date between 3000 BC and 3000 AD')
    try:
        TT = dt.date.fromisoformat(input('Date yyyy-mm-dd: '))  # Terrestrial time, approx JPL ephemeris time
    except ValueError:
        print('Invalid date format')
        exit()
    start = (dt.date.today() - J2000).days                      # Days between today and J2000
    end = (TT - J2000).days                                     # Days between input final date and J2000
    return [start, end]

def keplerelements(planet, t):
    # Propagation of Kepler orbital elements using Standish's linear fit
    a = planet['a'] + t * planet['da/dt']
    e = planet['e'] + t * planet['de/dt']
    w = np.deg2rad(planet['w'] + t * planet['dw/dt'])
    I = np.deg2rad(planet['I'] + t * planet['dI/dt'])
    W = np.deg2rad(planet['W'] + t * planet['dW/dt'])
    if planet.name in jovians.index:                            # Jovian planets correction terms
        jovian = np.deg2rad(jovians.loc[planet.name])
        b = jovian['b']; c = jovian['c']; s = jovian['s']; f = jovian['f']
        M = M + b * t ** 2 + c * np.cos(f * t) + s * np.sin(f * t)
    K = lambda E: E - e * np.sin(E) - M                         # Kepler's equation = 0
    dKdE = lambda E: 1 - e * np.cos(E)                          # Derivative
    E = scp.optimize.newton(K, planet['E'], fprime=dKdE, tol=1.75e-4, maxiter=100)
    return [a, e, w, I, W, M, E]

def plotorbit(planet, ax):
    a = planet['a']; e = planet['e']; w = planet['w']; I = planet['I']; W = planet['W']; E = planet['E']
    # Coords. in heliocentric orbital plane
    x = np.array([[a * (np.cos(E) - e)],
                  [a * np.sqrt(1 - e ** 2) * np.sin(E)]])
    fullrot = np.linspace(0, 2 * np.pi, 100)
    ellipse = np.vstack((a * (np.cos(fullrot) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(fullrot)))
    # Coords. in J2000 ecliptic plane, rotation matrix (anticlock z (w), anticlock x (I), anticlock z (W))
    cw = np.cos(w); sw = np.sin(w)
    cI = np.cos(I); sI = np.sin(I)
    cW = np.cos(W); sW = np.sin(W)
    rot = np.array([[cw * cW - sw * sW * cI, -sw * cW - cw * sW * cI],
                    [cw * sW + sw * cW * cI, -sw * sW + cw * cW * cI],
                    [sw * sI, cw * sI]])
    xecl = rot @ x
    ellipsecl = rot @ ellipse
    ax.plot(xecl[0], xecl[1], xecl[2], 'o', label=planet.name)
    ax.plot(ellipsecl[0], ellipsecl[1], ellipsecl[2], 'k', lw=1)

def main():
    start, end = getdate()
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    frames = abs(end - start)
    print(keplerel.loc['Mercury', 'E'])
    for n in range(frames):
        plt.cla()
        ax.set_xlim3d(-25, 25)
        ax.set_ylim3d(-25, 25)
        ax.set_zlim3d(-25, 25)
        ax.set_xlabel(r'$\hat{X}$: Vernal Equinox (AU)', labelpad=10)
        ax.set_ylabel(r'$\hat{Y}=\hat{Z}\times\hat{X}$ (AU)', labelpad=10)
        ax.set_zlabel(r'$\hat{Z}$: North Ecliptic Pole (AU)', labelpad=10)
        ax.set_title(f'Solar System in J2000 Ecliptic Plane, noon {dt.date.today()}(TT)')
        t = start + n * np.sign(end - start)
        for idx, planet in keplerel.iterrows():
            keplerel.loc[idx, ['a', 'e', 'w', 'I', 'W','E']] = keplerelements(planet, t)
            plotorbit(planet, ax)
        ax.axes.legend(bbox_to_anchor=[1.8, .9])
        plt.savefig(f"frames/{n}.png")
    images = [Image.open(f"frames/{n}.png") for n in range(frames)]
    images[0].save('solar_system.gif', save_all=True, append_images=images[1:], duration=100, loop=0)
    os.system('pycharm64 solar_system.gif')

if __name__=="__main__":
    main()