import numpy as np
import pandas as pd
import datetime as dat
import webbrowser
import matplotlib.lines
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use(['dark_background'])

# E.M. Standish (1992)'s data for simplified orbit propagation
orbitalElements = pd.read_csv('orbitel.csv', index_col=0) 
joviansCorrection = pd.read_csv('jovians.csv', index_col=0) 

# Change degrees to radians
orbitalElements.loc[:, 'I':'d/dt(long.peri.)'] = np.pi / 180 * orbitalElements.loc[:, 'I':'d/dt(long.peri.)']
joviansCorrection = np.pi / 180 * joviansCorrection

# Classical Kepler orbital elements
# Argument of perihelion
orbitalElements['w'] = orbitalElements['long.peri.'] - orbitalElements['long.node.'] 
orbitalElements['dw/dt'] = orbitalElements['d/dt(long.peri.)'] - orbitalElements['d/dt(long.node.)']
# Mean anomaly
orbitalElements['M'] = orbitalElements['L'] - orbitalElements['long.peri.'] 
orbitalElements['dM/dt'] = orbitalElements['dL/dt'] - orbitalElements['d/dt(long.peri.)']

# Separating out the angles actually used for computing positions (Kepler elements), their derivatives, and correction terms
KEPLER_ELEMENTS = orbitalElements[['a', 'e', 'w', 'I', 'long.node.', 'M']].to_numpy()
VARIATIONS = orbitalElements[['da/dt', 'de/dt', 'dw/dt', 'dI/dt', 'd/dt(long.node.)', 'dM/dt']].to_numpy()
CORRECTIONS = np.zeros((len(orbitalElements.index), 4))
for idx, planet_name in enumerate(orbitalElements.index):
        # Get corrections for planetary orbits
        if planet_name in joviansCorrection.index:
            CORRECTIONS[idx, :] = joviansCorrection.loc[planet_name].to_numpy()
        else:
            CORRECTIONS[idx, :] = np.zeros(4)

J2000 = dat.date(2000, 1, 1) # Reference epoch (1st July 2000, 12:00 GMT)
TOL = 1.75e-8 # Tolerance for Kepler equation solution (radians)
ELLIPSEPOINTS = 50 # Orbit line resolution (point count)
ANIDPI = 100
ANIFPS = 50
ANIFRAMES = 150
fig, ax = plt.subplots(figsize=(13, 10), subplot_kw={'projection':'3d'})
AXLIMS = orbitalElements['a'].max()
COLOURS = ['darkgrey', 'darkorange', 'g', 'r', 'darksalmon', 'gold', 'lightblue', 'b']
planet_dots = [ax.plot([], [], [], 'o', label=planet, color=COLOURS[idx])[0] 
               for idx, planet in enumerate(orbitalElements.index)]
orbit_lines = [ax.plot([], [], [], 'k', lw=1)[0] for _ in orbitalElements.index]
time_text = fig.text(.612, .889, '', fontsize=20)

def start_n_end() -> tuple[dat.date]:
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

def calcorbit(angles: np.array, correction: np.array, t: float) -> tuple[np.array]:
    # Propagation of Kepler orbital elements using Standish's linear fit
    # Units: centuries, astronomical units, radians
    a, e, w, I, O, M = angles
    b, c, s, f = correction
    M += b * t ** 2 + c * np.cos(f * t) + s * np.sin(f * t)
    # Modulus M between -π and π rad
    M = M % (2 * np.pi)
    if M > np.pi:
        M -= 2 * np.pi
    if M < -np.pi:
        M += 2 * np.pi
    # Newton's method to find roots of Kepler's equation
    E = M
    f = M - (E - e * np.sin(E))
    while np.abs(f) > TOL:
        E += f / (1 - e * np.cos(E))
        f = M - (E - e * np.sin(E))
    # Coords. in heliocentric orbital plane (2D)
    x = np.array([[a * (np.cos(E) - e)],
                  [a * np.sqrt(1 - e ** 2) * np.sin(E)]])
    two_pi_radians = np.linspace(0, 2 * np.pi, ELLIPSEPOINTS)
    ellipse = np.vstack((a * (np.cos(two_pi_radians) - e),
                         a * np.sqrt(1 - e ** 2) * np.sin(two_pi_radians)))
    # Coords. in J2000 ecliptic plane (3D), rotation matrix (Rz(O) * Rx(i) * Rz(w))
    cw, sw = np.cos(w), np.sin(w)
    ci, si = np.cos(I), np.sin(I)
    co, so = np.cos(O), np.sin(O)
    rot = np.array([[cw * co - sw * so * ci, -sw * co - cw * so * ci],
                    [cw * so + sw * co * ci, -sw * so + cw * co * ci],
                    [sw * si, cw * si]])
    xecl = rot @ x
    ellipsecl = rot @ ellipse
    return xecl, ellipsecl

def get_coords(
    t: float, 
    elements: np.array, 
    variations: np.array, 
    corrections: np.array
) -> tuple[np.array]:
    angles = elements + t * variations
    xecl, ellipsecl = calcorbit(angles[0, :], corrections[0, :], t)
    for i, angles in enumerate(angles[1:, :]):
        # Get coordinates for planets and orbits
        a, b = calcorbit(angles, corrections[i, :], t)
        xecl = np.hstack((xecl, a))
        ellipsecl = np.dstack((ellipsecl, b))
    return xecl, ellipsecl

def init_func() -> tuple[matplotlib.lines.Line2D]:
    # Initialise figure
    ax.set_xlim3d(-AXLIMS, AXLIMS)
    ax.set_ylim3d(-AXLIMS, AXLIMS)
    ax.set_zlim3d(-AXLIMS, AXLIMS)
    ax.set_xlabel(r'$\hat{X}$: Vernal Equinox (AU)', labelpad=10, fontsize=15)
    ax.set_ylabel(r'$\hat{Y}=\hat{Z}\times\hat{X}$ (AU)', labelpad=10, fontsize=15)
    ax.set_zlabel(r'$\hat{Z}$: North Ecliptic Pole (AU)', labelpad=10, fontsize=15)
    ax.set_title('Solar System, ECLIPJ2000 reference frame,', fontsize=20)
    plt.subplots_adjust(left=-0.15)
    ax.legend(bbox_to_anchor=(1.5, .8), fontsize=15)
    return *planet_dots, *orbit_lines

def animate(t: float) -> tuple[matplotlib.lines.Line2D]:
    xecl, ellipsecl = get_coords(t, KEPLER_ELEMENTS, VARIATIONS, CORRECTIONS)
    for idx in range(orbitalElements.shape[0]):
        # Plot planets
        planet_dots[idx].set_data_3d(xecl[0, idx], xecl[1, idx], xecl[2, idx])
        ## Plot orbits
        orbit_lines[idx].set_data_3d(ellipsecl[0, :, idx], ellipsecl[1, :, idx], ellipsecl[2, :, idx])
    # Display date
    try:
        date = J2000 + dat.timedelta(days=t*36525)
        era = 'AD'
    except OverflowError:
        date = dat.date(1, 1, 1) - dat.timedelta(days=t*36525) - (J2000 - dat.date(1, 1, 1))
        era = 'BC'
    time_text.set_text(f'{date} {era} (TT)')
    return *planet_dots, *orbit_lines

def main():
    start, end = start_n_end()
    frames = np.concatenate((start * np.ones(int(ANIFPS / 2)),  # Linger on starting frame
                             np.linspace(start, end, ANIFRAMES - ANIFPS),
                             end * np.ones(int(ANIFPS / 2))))   # Linger on final frame
    ani = FuncAnimation(fig, animate, frames=frames, init_func=init_func, interval=50, blit=True)
    ani.save('solar_system.gif', writer='pillow', fps=ANIFPS, dpi=ANIDPI)
    webbrowser.open('solar_system.gif')

if __name__=="__main__":
    main()