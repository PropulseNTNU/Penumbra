import control
import kalman_filter
import interpolation
import math
from numpy import genfromtxt


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

def test_kalman():
    data = genfromtxt('data/flight1.csv', delimiter=',')
    estimates_h = []
    estimates_v = []
    prev_time = data[0][0] - 30 #just to have a dt at the start
    for row in data:
        print(row[0])
        dt = (row[0] - prev_time) / 1000
        est_h, est_v = kalman_filter.kalman(row[4], row[21], dt)
        estimates_h.append(est_h.item(0))
        estimates_v.append(est_v.item(0))
        prev_time = row[0]

    print("Height: ")
    print(estimates_h)
    print("Velocity: ")
    print(estimates_v)

test_kalman()