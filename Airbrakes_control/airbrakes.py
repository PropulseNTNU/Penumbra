import control
import kalman_filter
import interpolation
import math

ki=0.008
kp=10
rieman_sum=0
prev_u = 0
default_rotation = 30
def airbrakes_main(height, velocity, dt):
    global rieman_sum
    global prev_u
    #estimated_h, estimated_v=kalman_filter.kalman(height, acc, dt)
    v_ref=interpolation.get_reference_velocity(height)
    error=float(velocity-v_ref)
    u, rieman_sum = control.controller(error, kp, ki, rieman_sum, dt)
    u=float(u)
    u += 30
    #prev_u = control.test_modifications(u, prev_u, dt)
    return u, error
