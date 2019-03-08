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
def launch(sample_rate, position, orientation, COM, COP, thrust, gravity, lift, drag):
    print("Initializing...")
    for step in range(len(position)):
        position[step][2] = -position[step][2]

    force_scale = 0.01
    l = 2
    rad = 0.09

    n = 1;
    while sample_rate // n > 30:
        n += 1
    sample_rate = sample_rate // n

    print("Adjusting resolution")

    position = position[::n]
    orientation = orientation[::n]
    COM = COM[::n]
    COP = COP[::n]
    thrust = thrust[::n]
    gravity = gravity[::n]
    lift = lift[::n]
    drag = drag[::n]

    steps = len(position)

    # Initialization of enviorment
    print("Initializing render enviorment")
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

    #COM_sphere = vp.sphere(radius = 0.04, color=vp.vector(0, 0, 0))
    #COP_sphere = vp.sphere(radius = 0.04, color=vp.vector(0, 0, 0))
    #thrust_pointer = vp.arrow(shaftwidth = 0.1, color=vp.vector(1, 1, 0))
    #gravity_pointer = vp.arrow(shaftwidth = 0.1, color=vp.vector(0, 1, 1))
    #lift_pointer = vp.arrow(shaftwidth = 0.1, color=vp.vector(1, 0, 1))
    #drag_pointer = vp.arrow(shaftwidth = 0.1, color=vp.vector(0, 0, 1))

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

    print("Placing references...")
    for step in range(2, steps):
        if not step % 10:
            a = 4
            c = 0.3
            b = a - c
            ref = vp.extrusion(path=[vp.vector(0, 0, 0), vp.vector(c, 0, 0)], shape=[sqr_out, sqr_in1, sqr_in2], axis=vp.vector(1, 0, 0), up=vp.vector(0, 1, 0), pos=vp.vector(position[step][0], position[step][1], position[step][2]))
            ref.up = vp.vector(0, 0, 1)
            ref.rotate(angle=orientation[step][0], axis=vp.cross(ref.axis, ref.up))
            ref.rotate(angle=orientation[step][1], axis=ref.up)
            ref.rotate(angle=orientation[step][2], axis=ref.axis)

    print("Beginning graphics loop")
    render.caption = "Running"
    render.visible = True
    t = 0
    for step in range(steps):
        t += 1/sample_rate
        render.caption = "Simulated time %.2f" %  t
        x = position[step][0]
        y = position[step][1]
        z = position[step][2]

        rocket.normal = vp.cross(rocket.axis, rocket.up)

        pitch, yaw, roll = orientation[step][0], orientation[step][1], orientation[step][2]

        #COM_x = x + COM[step] * np.cos(yaw) * np.cos(pitch)
        #COM_y = y + COM[step] * np.sin(yaw) * np.cos(pitch)
        #COM_z = z + COM[step] * np.sin(pitch)

        #COP_x = x + COP[step] * np.cos(yaw) * np.cos(pitch)
        #COP_y = y + COP[step] * np.sin(yaw) * np.cos(pitch)
        #COP_z = z + COP[step] * np.sin(pitch)

        #thrust_x = x + (-l) * np.cos(yaw) * np.cos(pitch)
        #thrust_y = y + (-l) * np.sin(yaw) * np.cos(pitch)
        #thrust_z = z + (-l) * np.sin(pitch)

        #thrust_mag = np.linalg.norm(thrust[step]) * force_scale
        #thrust_ax_x = thrust_mag * np.cos(yaw) * np.cos(pitch)
        #thrust_ax_y = thrust_mag * np.sin(yaw) * np.cos(pitch)
        #thrust_ax_z = thrust_mag * np.sin(pitch)

        #lift_mag = lift[step][0] * force_scale
        #lift_ax_x = 0
        #lift_ax_y = 0
        #lift_ax_z = lift_mag

        #drag_mag = np.linalg.norm(drag[step]) * force_scale
        #drag_ax_x = drag_mag
        #drag_ax_y = 0
        #drag_ax_z = 0

        #gravity_mag = np.linalg.norm(gravity[step]) * force_scale

        #thrust_pointer.pos = vp.vector(thrust_x, thrust_y, thrust_z)
        #thrust_pointer.axis = vp.vector(thrust_ax_x, thrust_ax_y, thrust_ax_z)

        #gravity_pointer.pos = vp.vector(COM_x, COM_y, COM_z)
        #gravity_pointer.axis = vp.vector(0, 0, -gravity_mag)

        #lift_pointer.pos = vp.vector(COP_x, COP_y, COP_z)
        #lift_pointer.axis = vp.vector(lift_ax_x, lift_ax_y, lift_ax_z)

        #drag_pointer.pos = vp.vector(COP_x, COP_y, COP_z)
        #drag_pointer.axis = vp.vector(drag_ax_x,drag_ax_y, drag_ax_z)

        #COM_sphere.pos = vp.vector(COM_x, COM_y, COM_z)
        #COP_sphere.pos = vp.vector(COP_x, COP_y, COP_z)

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
