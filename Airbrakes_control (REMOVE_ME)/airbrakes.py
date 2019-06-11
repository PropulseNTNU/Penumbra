import control
import kalman_filter
import interpolation
import math

ki=0.01
kp=1
rieman_sum=0
prev_u = 0
def airbrakes_main(height, acc, dt):
    global rieman_sum
    global prev_u
    estimated_h, estimated_v=kalman_filter.kalman(height, acc, dt)
    v_ref=interpolation.get_reference_velocity(estimated_h)
    error=estimated_v-v_ref
    u, rieman_sum = control.controller(error, kp, ki, rieman_sum, dt)
    prev_u = control.test_modifications(u, prev_u, dt)
    u = prev_u
    if u<0:
        u=0
    elif u>90:
        u=90
 
    return control.calculate_area(u), estimated_h, estimated_v, error
