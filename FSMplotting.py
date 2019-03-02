
import numpy as np
import matplotlib.pyplot as plt

# You probably won't need this if you're embedding things in a tkinter plot...
plt.ion()
axes = None
fig = None
graphs = None

def init(teensyData, cols):
    global axes
    global fig
    global graphs
    rows = int(np.ceil(len(teensyData)/cols))
    fig, axes = plt.subplots(nrows=rows, ncols=cols)
    axes = axes.flatten()
    graphs = []
    for i, (ax, data) in enumerate(zip(axes, teensyData.items())):
        graphs.append(ax.plot([], [], 'black')[0])
        ax.set_title(data[1][0])
    fig.canvas.draw()
    plt.subplots_adjust(hspace=1)
    plt.subplots_adjust(wspace=0.5)
    
def plotData(teensyData, timeData):

    for i, (graph, val) in enumerate(zip(graphs, teensyData.items())):
        data = val[1][1]
        graph.set_data(timeData, data)
        axes[i].draw_artist(axes[i].patch)
        axes[i].draw_artist(graph)
        axes[i].set_xlim(0, timeData[-1] + 2) 
        axes[i].set_ylim(min(data)-50, max(data)+50)
        if val[0] =="c_s":
            axes[i].set_ylim(min(data)-0.001, max(data) + 0.001)
    fig.canvas.blit(axes[0].bbox)
    fig.canvas.flush_events()