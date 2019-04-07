import control
import kalman_filter
import interpolation
import math
#dt=0.03
ki=0.01
kp=1
rieman_sum=0
flaps_widt=0.2
flaps_max_extention=0.04
temp_rieman_sum=0;
def airbrakes_main(height, acc, dt):
    global rieman_sum
    estimated_h, estimated_v=kalman_filter.kalman(height, acc, dt)
    v_ref=interpolation.get_reference_velocity(estimated_h)
    error=estimated_v-v_ref
    u, temp_rieman_sum = control.controller(error, ki, kp, rieman_sum, dt)
    rieman_sum=temp_rieman_sum
    #rieman_sum må oppdateres
    if u<0:
        u=0
    elif u>90*(math.pi)/180:
        u=90*(math.pi)/180
    else:
        u*=((math.pi)/180)
    return control.calculate_area(u), estimated_h, estimated_v, error
