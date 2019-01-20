import numpy as np
import vpython as vp

rad2deg = 180/np.pi
deg2rad = np.pi/180

def np2vp(np_array):
    return vp.vector(np_array[0], np_array[1], np_array[2])

def Ryzx(pitch, yaw, roll):
    cp = np.cos(pitch)
    sp = np.sin(pitch)
    cy = np.cos(yaw)
    sy = np.sin(yaw)
    cr = np.cos(roll)
    sr = np.sin(roll)
    R = np.array([[cp*cy, sp*sr-cp*sy*cr, sp*cr+cp*sy*sr],
    [sy, cy*cr, -cy*sr], [-sp*cy, cp*sr+sp*sy*cr, cp*cr-sp*sy*sr]])
    return R

def rel_euler(pitch_1, yaw_1, roll_1, pitch_2, yaw_2, roll_2):
    R = Ryzx(pitch_1, yaw_1, roll_1).T @ Ryzx(pitch_2, yaw_2, roll_2)
    pitch = np.arctan2(-R[2,0], R[0,0])
    yaw = np.arcsin(R[1,0])
    roll = np.arctan2(-R[1,2], R[1,1])
    return (pitch, yaw, roll)

# Rocket in world frame
def launch(sample_rate, position, orientation):
    for step in range(len(position)):
        position[step][2] = -position[step][2]

    l = 2
    rad = 0.09
    steps = len(position)

    # Initialization of enviorment
    render = vp.canvas(height=600, width=1200, background=vp.vector(0.8, 0.8, 0.8), forward=vp.vector(1, 0, 0), up=vp.vector(0, 0, 1))
    render.select()
    render.caption = "Loading..."
    render.visible = False

    # Initialization of body and cone
    body = vp.cylinder(pos=vp.vector(-l, 0, 0), axis=vp.vector(l*0.8, 0, 0), radius=rad)
    cone = vp.cone(pos=vp.vector(-l*0.2, 0, 0), axis=vp.vector(l*0.2, 0, 0), radius=rad)

    # Initialization of fins
    a_fins = vp.box(pos=vp.vector(-l + (l*0.05), 0, 0), size=vp.vector(l*0.1, rad*4, rad*0.25))
    b_fins = vp.box(pos=vp.vector(-l + (l*0.05), 0, 0), size=vp.vector(l*0.1, rad*0.25, rad*4))

    inv_box = vp.box(pos=vp.vector(l/2, 0, 0), size=vp.vector(l, 10e-10, 10e-10), visible=False)

    # Initialization of rocket
    rocket = vp.compound([body, cone, a_fins, b_fins, inv_box], pos=vp.vector(0, 0, 1), axis=vp.vector(1, 0, 0), up=vp.vector(0, 0, 1), color=vp.vector(1,1,1))

    a = 4
    c = 0.3
    b = a - c
    sqr_out = [[-a, a],[a, a],[a, -a],[-a, -a], [-a, a]]
    sqr_in1 = [[-b, b - 3*c],[b, b - 3*c],[b, -b],[-b, -b], [-b, b - 3*c]]
    sqr_in2 = [[-b, b],[b, b], [b, b - 2*c],[-b, b - 2*c], [-b, b]]
    ref_box = vp.extrusion(path=[vp.vector(-c/2, 0, 0), vp.vector(c/2, 0, 0)], shape=[sqr_out, sqr_in1, sqr_in2], color=vp.vector(0.5, 0.5, 0.5))

    render.camera.follow(rocket)
    render.autoscale = False
    launch_pad = vp.box(pos=vp.vector(position[0][0], position[0][1], -1), size=vp.vector(16, 16, 2), color=vp.vector(0.2, 0.2, 0.2))
    launch_pad = vp.box(pos=vp.vector(position[-1][0], position[-1][1], -1), size=vp.vector(16, 16, 2), color=vp.vector(0.2, 0.2, 0.2))

    rocket.normal = vp.cross(rocket.up, rocket.axis)

    #vp.attach_arrow(rocket, 'up', color=vp.color.green)
    #vp.attach_arrow(rocket, 'axis', color=vp.color.blue)
    #vp.attach_arrow(rocket, 'normal', color=vp.color.red)

    roll = 0
    vp.attach_trail(rocket, radius=rad/3, color=vp.color.red)

    sqr_out = [[-a, a],[a, a],[a, -a],[-a, -a], [-a, a]]
    sqr_in1 = [[-b, b - 3*c],[b, b - 3*c],[b, -b],[-b, -b], [-b, b - 3*c]]
    sqr_in2 = [[-b, b],[b, b], [b, b - 2*c],[-b, b - 2*c], [-b, b]]

    for step in range(1, steps):
        if not step % 5:
            a = 4
            c = 0.3
            b = a - c
            ref = vp.extrusion(path=[vp.vector(0, 0, 0), vp.vector(c, 0, 0)], shape=[sqr_out, sqr_in1, sqr_in2], axis=vp.vector(1, 0, 0), up=vp.vector(0, 1, 0), pos=vp.vector(position[step][0], position[step][1], position[step][2]))
            ref.up = vp.vector(0, 0, 1)
            ref.rotate(angle=orientation[step][0], axis=vp.cross(ref.axis, ref.up))
            ref.rotate(angle=orientation[step][1], axis=ref.up)
            ref.rotate(angle=orientation[step][2], axis=ref.axis)

    render.caption = "Running"
    render.visible = True
    for step in range(steps):
        x = position[step][0]
        y = position[step][1]
        z = position[step][2]

        rocket.normal = vp.cross(rocket.axis, rocket.up)

        pitch, yaw, roll = orientation[step][0], orientation[step][1], orientation[step][2],

        rocket.rotate(angle=pitch, axis=rocket.normal)
        rocket.rotate(angle=yaw, axis=rocket.up)
        rocket.rotate(angle=roll, axis=rocket.axis)

        rocket.pos = vp.vector(x, y, z)
        vp.sleep(1/sample_rate)
        if step + 1 != steps:
            rocket.rotate(angle=-roll, axis=rocket.axis)
            rocket.rotate(angle=-yaw, axis=rocket.up)
            rocket.rotate(angle=-pitch, axis=rocket.normal)

    render.caption = "Done"

    # Vloop
