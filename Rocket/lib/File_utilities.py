"""
This is a module that contains functions for file reading.

Last edit: 16.11.2018
"""
from numpy import linspace, reshape, flip, loadtxt, sin, pi, unique, stack


def find_parameter(file, parameter):
    File = open(file, 'r')
    arr = ["", ""]
    while arr[0] != parameter.lower():
        base = File.readline()
        if base == '':
            print("ERROR: Could not find parameter '" + parameter + "' in '" + file + "'.")
            return False
        base = base.replace(" ", "")
        base = base.replace("\n", "")
        arr = base.split("=")
    File.close()
    return arr[1]


def unwrap_report1(file, T, alpha_max, delta_v, v0):
    f = 1/T
    report = loadtxt(file, skiprows=7).transpose()[1:]
    c = len(report.transpose())
    c -= c%T

    def make2d(report, c, T, idx):
        arr = report[idx]
        arr = arr[:c]
        arr = reshape(arr, (-1, T))
        for i in range(len(arr)):
            if i%2 == 0:
                arr[i] = flip(arr[i])
        return arr.transpose()

    drag = make2d(report, c, T, 0)
    lift = make2d(report, c, T, 1)
    moment = make2d(report, c, T, 2)

    # Speed axis seems to be working fine regardless of T
    # Alpha axis acts weird when given other integers than 5
    alpha_axis = alpha_max*sin(linspace(0, T, T)*pi*f/2)
    speed_axis = v0 + (linspace(0, c/T - 1, c//T)*delta_v)

    return alpha_axis, speed_axis, drag, lift, moment


def unwrap_CFD_report(file):
    """
    Reads CFD file with CFD values (drag, lift or moment about COM) 
        *** Assuming 1 airspeed per column, and 1 AoA per row
        *** Assuming that there are only floats in the file.

    :param file: The CFD file [File]
    :return: CFD_data [2D array] 
    """
    CFD = loadtxt(file, skiprows=1, dtype=float)

    return CFD


def unwrap_openRocket(file):
        # Assuming file is in a proper format. See file in Tests/V9 folder
        f = loadtxt(file, skiprows=1, delimiter=',')

        # Longitudinal MOI is Iyy and Izz
        time, mass, propMass, I_longitudinal, Ixx, COM = f[:,0], f[:,1], f[:,2], f[:,3], f[:,4], f[:,5]

        return time, mass, propMass, COM, Ixx, I_longitudinal