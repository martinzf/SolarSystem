import numpy as np
import scipy as scp
import pandas as pd
import datetime as dat
import os
import matplotlib.pyplot as plt
from PIL import Image
plt.style.use(['science', 'notebook', 'grid'])

# E.M. Standish (1992)'s data for simplified orbit propagation
keplerel = pd.read_csv('keplerel.csv', index_col=0)             # Planets' Kepler elements
jovians = pd.read_csv('jovians.csv', index_col=0)               # Jovian planets' correction terms

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
    M = M % (2 * np.pi)                                         # Keep M between -180º and 180º
    if M > np.pi:
        M = M - 2 * np.pi
    keplerel.loc[planet.name, 'M'] = M
    K = lambda E: E - e * np.sin(E) - M                         # Kepler's equation = 0
    dKdE = lambda E: 1 - e * np.cos(E)                          # Derivative
    E = scp.optimize.newton(K, planet['E'], fprime=dKdE, tol=1.75e-8)
    keplerel.loc[planet.name, 'E'] = E

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

def plotssyst(dt, datetime, frame, ax):
    plt.cla()
    ax.set_xlim3d(-25, 25)
    ax.set_ylim3d(-25, 25)
    ax.set_zlim3d(-25, 25)
    ax.set_xlabel(r'$\hat{X}$: Vernal Equinox (AU)', labelpad=10)
    ax.set_ylabel(r'$\hat{Y}=\hat{Z}\times\hat{X}$ (AU)', labelpad=10)
    ax.set_zlabel(r'$\hat{Z}$: North Ecliptic Pole (AU)', labelpad=10)
    for _, planet in keplerel.iterrows():
        keplerelements(planet, dt / 36525)
        plotorbit(planet, ax)
    ax.set_title(f'Solar System in J2000 Ecliptic Plane, {datetime} (TT)')
    ax.axes.legend(bbox_to_anchor=[1.8, .9])
    plt.savefig(f"frames/{frame}.png", bbox_inches='tight')

def main():
    start, end = getdate()
    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    dt = start
    datetime = dat.datetime.combine(dat.date.today(), dat.time(12))
    plotssyst(dt, datetime, 0, ax)
    frames = 31
    dt = (end - start) / (frames - 1)
    for frame in np.arange(1, frames):
        datetime += dat.timedelta(days=dt)
        plotssyst(dt, datetime, frame, ax)
    images = [Image.open(f"frames/{n}.png") for n in range(frames)]
    images[0].save(fp='solar_system.gif', format='GIF', save_all=True, append_images=images[1:], duration=17, loop=0)

if __name__=="__main__":
    main()