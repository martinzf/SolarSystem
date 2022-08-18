import numpy as np
import pandas as pd
import datetime as dat
import webbrowser
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use(['dark_background'])

# E.M. Standish (1992)'s data for simplified orbit propagation
keplerel = pd.read_csv('keplerel.csv', index_col=0)             # Planets' Kepler elements
jovians = pd.read_csv('jovians.csv', index_col=0)               # Jovian planets' correction terms

J2000 = dat.date(2000, 1, 1)                                    # Reference epoch (1st July 2000, 12:00 GMT)
tol = 1.75e-8                                                   # Tolerance for Kepler equation solution (radians)
ellipsepoints = 50                                              # Orbit line resolution (point count)
anidpi = 100
anifps = 50
aniframes = 150
fig, ax = plt.subplots(figsize=(13, 10), subplot_kw={'projection':'3d'})
axlims = keplerel['a'].max()
colours = ['darkgrey', 'darkorange', 'g', 'r', 'darksalmon', 'gold', 'lightblue', 'b']
lns1 = [ax.plot([], [], [], 'o', label=planet, color=colours[idx])[0] for idx, planet in enumerate(keplerel.index)]
lns2 = [ax.plot([], [], [], 'k', lw=1)[0] for planet in keplerel.index]
time_text = fig.text(.612, .889, '', fontsize=20)

def getdate():
    print('Input a date between 3000 BC and 3000 AD')
    # Get terrestrial time (TT), approx. JPL ephemeris time
    while True:
        try:
            TT = dat.date.fromisoformat(input('Date yyyy-mm-dd: '))
            if TT > dat.date(3000, 12, 31):
                print('Date outside specified range')
            else:
                break
        except ValueError:
            print('Invalid date format')
    # Centuries between today and J2000
    start = (dat.date.today() - J2000).days / 36525
    # Centuries between input final date and J2000
    while True:
        era = input('AD/BC: ').lower()
        if era == 'ad' or era == 'bc':
            break
        else:
            print('Must input "AD" or "BC"')
    if era == 'ad':
        end = (TT - J2000).days / 36525
    else:
        end = - ((TT - dat.date(1, 1, 1)) + (J2000 -dat.date(1, 1, 1))).days / 36525
    return start, end

def calcorbit(elements, correction, t):
    # Propagation of Kepler orbital elements using Standish's linear fit
    # Units: centuries, astronomical units, radians
    a, e, w, i, O, M = elements[::2] + t * elements[1::2]
    b, c, s, f = correction
    M += b * t ** 2 + c * np.cos(f * t) + s * np.sin(f * t)
    # Modulus M between -π and π rad
    M = M % (2 * np.pi)
    if M > np.pi:
        M -= 2 * np.pi
    elif M < -np.pi:
        M += 2 * np.pi
    # Newton's method to find roots of Kepler's equation
    E = M
    f = M - (E - e * np.sin(E))
    while np.abs(f) > tol:
        E += f / (1 - e * np.cos(E))
        f = M - (E - e * np.sin(E))
    # Coords. in heliocentric orbital plane (2D)
    x = np.array([[a * (np.cos(E) - e)],
                  [a * np.sqrt(1 - e ** 2) * np.sin(E)]])
    fullrot = np.linspace(0, 2 * np.pi, ellipsepoints)
    ellipse = np.vstack((a * (np.cos(fullrot) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(fullrot)))
    # Coords. in J2000 ecliptic plane (3D), rotation matrix (Rz(O) * Rx(i) * Rz(w) * v)
    cw, sw = np.cos(w), np.sin(w)
    ci, si = np.cos(i), np.sin(i)
    cO, sO = np.cos(O), np.sin(O)
    rot = np.array([[cw * cO - sw * sO * ci, -sw * cO - cw * sO * ci],
                    [cw * sO + sw * cO * ci, -sw * sO + cw * cO * ci],
                    [sw * si, cw * si]])
    xecl = rot @ x
    ellipsecl = rot @ ellipse
    return xecl, ellipsecl

def init():
    # Initialise figure
    ax.set_xlim3d(-axlims, axlims)
    ax.set_ylim3d(-axlims, axlims)
    ax.set_zlim3d(-axlims, axlims)
    ax.set_xlabel(r'$\hat{X}$: Vernal Equinox (AU)', labelpad=10, fontsize=15)
    ax.set_ylabel(r'$\hat{Y}=\hat{Z}\times\hat{X}$ (AU)', labelpad=10, fontsize=15)
    ax.set_zlabel(r'$\hat{Z}$: North Ecliptic Pole (AU)', labelpad=10, fontsize=15)
    ax.set_title('Solar System, ECLIPJ2000 reference frame,', fontsize=20)
    plt.subplots_adjust(left=-0.15)
    ax.legend(bbox_to_anchor=(1.5, .8), fontsize=15)
    return *lns1, *lns2

def animate(t):
    for idx, elements in enumerate(keplerel.to_numpy()):
        # Get coordinates for planets and orbits
        planet = keplerel.index[idx]
        if planet in jovians.index:
            correction = jovians.loc[planet].to_numpy()
        else:
            correction = np.zeros(4)
        xecl, ellipsecl = calcorbit(elements, correction, t)
        # Plot planets
        lns1[idx].set_data(xecl[0], xecl[1])
        lns1[idx].set_3d_properties(xecl[2])
        # Plot orbits
        lns2[idx].set_data(ellipsecl[0], ellipsecl[1])
        lns2[idx].set_3d_properties(ellipsecl[2])
    # Display date
    try:
        date = J2000 + dat.timedelta(days=t * 36525)
        era = 'AD'
    except OverflowError:
        date = dat.date(1, 1, 1) - dat.timedelta(days=t * 36525) - (J2000 - dat.date(1, 1, 1))
        era = 'BC'
    time_text.set_text(f'{date} {era} (TT)')
    return *lns1, *lns2

def main():
    start, end = getdate()
    frames = np.concatenate((start * np.ones(int(anifps / 2)),  # Linger on starting frame
                             np.linspace(start, end, aniframes - anifps),
                             end * np.ones(int(anifps / 2))))   # Linger on final frame
    ani = FuncAnimation(fig, animate, frames=frames, init_func=init, interval=50, blit=True)
    ani.save('solar_system.gif', writer='pillow', fps=anifps, dpi=anidpi)
    webbrowser.open('solar_system.gif')

if __name__=="__main__":
    main()