import math

#varibles for testing
flap_width=0.1106 # the with of the air brake flap in meters. This is only for testing, and will not be used during flight
max_extention=0.02 #the max length of the air brake flaps in meters. This is only for testing, and will not be used during flight


def integrate(prev_sum, value, step):
  return prev_sum + (value * step)

def pros_to_deg(prosent):
  return math.asin(prosent/100.0) * 180.0/math.pi


def controller(error, kp, ki, riemann_sum, dt):#PI-controller
  #add a lower and upper bound to prevernt overflow
  temp = integrate(riemann_sum, error, dt)
  if temp < 40 and temp > -40:
    riemann_sum = integrate(riemann_sum, error, dt)#integrates error

  prosent = kp*error + ki*riemann_sum
  if prosent > 86.6:
    return 86.6, riemann_sum

  elif prosent < 0.0:
    return 0, riemann_sum

  else:
    return prosent, riemann_sum


#functions for testing
def calculate_area(u):
  return 3*flap_width*max_extention*math.sin(math.radians(u))

def test_modifications(ref_u, prev_u, dt):#Calculates the actual actuation based on servospeed
  if ref_u==prev_u:
    return ref_u

  if ref_u>prev_u:
    prev_u+=(66.66/0.13)*dt#Calculates the servo position based on rotationspeed.
    if prev_u>ref_u:#In this case the servo has reaced its reference, and ref_u can be returned
      return ref_u
    return prev_u

  prev_u -= (66.66/0.13)*dt #Calculates the servo position based on rotationspeed.
  if prev_u<ref_u: #In this case the servo has reaced its reference, and ref_u can be returned
    return ref_u
  return prev_u
