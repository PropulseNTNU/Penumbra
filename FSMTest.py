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

fig, axes = plt.subplots(nrows=5)
graphs = []
for ax in axes:
    graphs.append(ax.plot(x, y, 'r-')[0])

fig.canvas.draw()

a = []
b = []
for phase in np.linspace(0, 10*np.pi, 500)
    a
    for graph in graphs:
        print(np.sin(x + phase))
        graph.set_ydata(np.sin(x + phase))
        axes[i].draw_artist(axes[i].patch)
        axes[i].draw_artist(graph)
        i += 1
    fig.canvas.draw_idle()
    fig.canvas.flush_events()