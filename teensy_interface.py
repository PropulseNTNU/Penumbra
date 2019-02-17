import serial
import serial.tools.list_ports
import time

def showPortList():
    if len(serial.tools.list_ports.comports()) == 0:
        print("No ports detected. Try using another USB port or reconnect it")
    for item in serial.tools.list_ports.comports():
        print(item)


def initSerial(port, baudrate, timeout):
    """
    Initialize the serial connection with the microcontroller.

    :param port: string; The serial port connected to the microcontroller.

    :param baudrate: int; The baud rate in bits per second. Needs to match the baud rate that the microcontroller uses.

    :param timeout: Float; The maximum time we wait to read X number of bytes or a whole line.

    :return: serial.Serial; returns the initialized serial object or None if it failed.
    """
    try:
        ser = serial.Serial(port, baudrate, timeout=timeout)
        print("Successfully initialized connection")
        return ser
    except (ValueError, serial.SerialException) as error:
        print("Teensy Interface: failed to initialize serial connection with the device")
        print("You tried to initialize port:", port)
        print("These are the available ports: ")
        showPortList()
        print(error)
        return None

def readControlSignal(ser, prefix='', size=0):
    """
    Read the control signal. This function alows other data than the control signal to be on the serial connection.
    Specify the prefix if you know there will be other data in addition to the control signal on the serial link.

    :param serial.Serial: The serial object already initialized.

    :param str prefix: The control signal will have a prefix if other data is on the link.
    Specify the prefix so we know what data to read.

    :param int size: How many bytes we should read before. This number should be more than the
    number of bytes  writen to the serial link in one loop by the air brakes program.

    :return int: returns the control signal or None if it cant find one.
    """
    if prefix != '':
        try:
            data = ser.read(size).decode("utf-8").split("\n")
            for i in range(len(data)-1, 0, -1):
                if data[i].startswith("c_s") and data[i].endswith("\r"):
                    return float(data[i].replace("c_s", '').replace("\r", ''))
            print("Did not find the control signal")

        except AttributeError as error:
            print("The serial connection is not initialized. Run the initSerial function first")
            print(error)
            return None
    else:
        try:
            return int.from_bytes(ser.read(), byteorder="big")

        except AttributeError as error:
            print("The serial connection is not initialized. Run the initSerial function first")
            print(error)
            return None

def sendHeightAndVelocity(ser, height, velocity):
    try:
        ser.write((("h"+ str(height) + "v" + str(velocity))).encode("utf-8"))
        return True
    except serial.SerialTimeoutException  as error:
        print("The write process timed out.")
        print(error)
        return False


def close(ser):
    ser.close()
