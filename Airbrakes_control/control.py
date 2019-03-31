
#include "controll.h"
#include <math.h>
import math

#varibles for testing
flap_width=0.1106 # the with of the air brake flap in meters. This is only for testing, and will not be used during flight
max_extention=0.02 #the max length of the air brake flaps in meters. This is only for testing, and will not be used during flight


def integrate(prev_sum, value, step):
  return prev_sum + (value * step)

def controller(error, ki, kp, riemann_sum, dt):#PI-controller
  #add a lower and upper bound to prevernt overflow
  riemann_sum = integrate(riemann_sum, error, dt)#integrates error
  return kp*error + ki*riemann_sum


#functions for testing
def calculate_area(u):
  return 3*flap_width*max_extention*math.sin(u)

def test_modifications(ref_u, prev_u, dt):#Calculates the actual actuation based on servospeed
  if ref_u==prev_u:
    return ref_u

  if ref_u>prev_u:
    prev_u+=(60/0.13)*dt#Calculates the servo position based on rotationspeed.
    if prev_u>ref_u:#In this case the servo has reaced its reference, and ref_u can be returned
      return ref_u
    return prev_u

  prev_u -= (60/0.13)*dt #Calculates the servo position based on rotationspeed.
  if prev_u<ref_u: #In this case the servo has reaced its reference, and ref_u can be returned
    return ref_u
  return prev_u
