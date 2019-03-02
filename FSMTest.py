"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time

fig, ax = plt.subplots()
line, = ax.plot(np.random.randn(100))
plt.show(block=False)
fig.canvas.draw()

tstart = time.time()
num_plots = 0
while time.time()-tstart < 5:
    line.set_ydata(np.random.randn(100))
    ax.draw_artist(ax.patch)
    ax.draw_artist(line)
    fig.canvas.draw_idle()
    fig.canvas.flush_events()
    num_plots += 1
print(num_plots/5)
"""
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 6*np.pi, 100)
y = np.sin(x)

# You probably won't need this if you're embedding things in a tkinter plot...
plt.ion()

fig, axes = plt.subplots(nrows=3, ncols=2)
axes = axes.flatten()
graphs = []
for ax in axes:
    graphs.append(ax.plot(x, y, 'r-')[0])

fig.canvas.draw()

a = [0]
b = [1]
for phase in np.linspace(0, 10*np.pi, 500):
    i = 0
    a.append(0)
    b.append(float(b[-1] + 0.03))
    for graph in graphs:
        graph.set_data(b, a)
        axes[i].draw_artist(axes[i].patch)
        axes[i].draw_artist(graph)
        axes[i].set_xlim(0, b[-1] + 5) 
        axes[i].set_ylim(min(a)-5, max(a)+5)
        i += 1
    fig.canvas.blit(ax.bbox)
    fig.canvas.flush_events()