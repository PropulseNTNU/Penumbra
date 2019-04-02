import control
import kalman_filter
import interpolation
import math
from numpy import genfromtxt
import matplotlib.pyplot as plt


dt=0.03
ki=0.01
kp=1
rieman_sum=0
flaps_widt=0.02
flaps_max_extention=0.04

def airbrakes_main(height, acc, dt):

    estimated_h, estimated_v=kalman_filter.kalman(height, acc, dt)
    v_ref=interpolation.get_reference_velocity(estimated_h)
    error=v_ref-estimated_v
    u=control.controller(error, ki, kp, rieman_sum, dt)
    #rieman_sum må oppdateres
    if u<0:
        u=0
    elif u>90*(math.pi)/180:
        u=90*(math.pi)/180
    else:
        u*=((math.pi)/180)
    return control.calculate_area(u), estimated_h, estimated_v, error

def preasureToHeight(preasure):
    T0 = 273.15+6#Temperature at sea level in K
    a = -0.0065 #Temperature decrease coefficient
    P0 = 1001.4 #Preasure at sea level
    R = 287.053 #Ideal gass constant
    g = 9.80665 #Acceleration of gravity
    h0 = 0 #Height at launch

    Analog2 = preasure*5/255
    sourceVoltage = 5.0;
    Analog2 = (1/0.009)*(Analog2/sourceVoltage+0.095)*1000; 
    P = Analog2/100;
    return (T0/a) * ((P/P0)**(-a * R/g) - 1) + h0



def test_kalman():
    data = genfromtxt('data/Launch_file_MAINTM.csv', delimiter=',')
    estimates_h = []
    estimates_v = []
    sensor_h = []
    sensor_acc = []
    prev_time = 0 #just to have a dt at the start
    time = []

    riemansum_v = [0]
    riemansum_h = [0]
    for row in data:
        dt = (row[0] - prev_time)
        time.append(row[0])
        acc_y = ((row[2]*5/255)- 2.55)/0.040
        riemansum_v.append(riemansum_v[-1] + (acc_y * dt))
        riemansum_h.append(riemansum_h[-1] + (riemansum_v[-1] * dt))
        h = preasureToHeight(row[3])
        sensor_acc.append(acc_y)
        sensor_h.append(h)
        est_h, est_v = kalman_filter.kalman(h, acc_y, dt)
        estimates_h.append(est_h.item(0))
        estimates_v.append(est_v.item(0))
        prev_time = row[0]


    plt.figure()
    plt.plot(time, estimates_h)

    plt.plot(time, estimates_v)

    plt.figure()
    plt.plot(time, sensor_acc)

    plt.plot(time, sensor_h)

    plt.figure()
    plt.plot(time, riemansum_v[1:])

    plt.figure()
    plt.plot(time, riemansum_h[1:])
    plt.show()

test_kalman()
