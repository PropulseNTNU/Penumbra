import numpy as np
import matplotlib.pyplot as plt

fig = plt.gcf()
fig.set_size_inches((12, 10))
fig.show()
plt.subplots_adjust(hspace=1)
plt.subplots_adjust(wspace=0.5)

def plotData(teensyData, timeData):
    i = 1
    rows = np.ceil(len(teensyData) / 2)
    for key, val in teensyData.items():
        plt.subplot(rows, 2, i)
        plt.plot(timeData[-2:], val[1][-2:], 'black')
        plt.title(val[0])
        plt.xlabel("Time(s)")
        i += 1
    fig.canvas.draw()

