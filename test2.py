import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d as ax3d
from matplotlib import animation
x=[1,2,3,4,5]
y=[1,2,3,4,5]
z=[1,2,3,4,5]

xj = [1,2,3,4,5]
yj = [2,4,6,8,10]
zj = [2,4,6,8,10]

N=len(x)
fig, ax = plt.subplots(subplot_kw={'projection':'3d'})

ax.set_xlim(0,10)
ax.set_ylim(0,10)
ax.set_zlim(0,10)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

lines = []
for index in range(2):
    lobj = ax.scatter3D([],[],[],'o')[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([],[])
        line.set_3d_properties([])
    return lines

def animate(i):
    xlist = [x[:i], xj[:i]]
    ylist = [y[:i], yj[:i]]
    zlist = [z[:i], zj[:i]]

    for lnum,line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum])
        line.set_3d_properties(zlist[lnum])
    time.sleep(1)
    return lines

anim = animation.FuncAnimation(fig, animate,init_func=init, frames=N , interval=5, blit=True)

plt.show()