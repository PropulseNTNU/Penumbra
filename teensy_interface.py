import serial
import serial.tools.list_ports
import time
import numpy as np

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

def readFloatData(ser, prefix='', lines = 1):
    """
    Read data that has the float datatype. This function alows other data than the control signal to be on the serial connection.
    Specify the prefix if you know there will be other data in addition to the data you want on the serial link.

    :param serial.Serial: The serial object already initialized.

    :param str prefix: The data needs to have a prefix if other data is on the link. 
    Specify the prefix so we know what data to read.  
    
    :return: returns ta float or None if it fails to find one.
    """
    if prefix != '':
        try:
            #print("Trying to read: ", prefix)
            lines_read = 0
            while lines_read < lines:
                lines_read += 1
                data = ser.readline().decode("utf-8")
                if prefix in data:
                    try:
                        retdata = float(data.replace(prefix, '').replace("\r", ''))
                        return retdata
                    except:
                        pass
            print("Read " + str(lines_read) + " lines without finding the data prefix.")
            return None
        except AttributeError as error:
            print("The serial connection is not initialized. Run the initSerial function first")
            print(error)
            return None
    else:
        try:
            return ser.readline().decode("utf-8")
                    
        except AttributeError as error:
            print("The serial connection is not initialized. Run the initSerial function first")
            print(error)
            return None

def readData(ser, prefix='', lines = 1):
    """
    Read data that has the float datatype. This function alows other data than the control signal to be on the serial connection.
    Specify the prefix if you know there will be other data in addition to the data you want on the serial link.

    :param serial.Serial: The serial object already initialized.

    :param str prefix: The data needs to have a prefix if other data is on the link. 
    Specify the prefix so we know what data to read.  
    
    :return: returns data or None if it fails to find one.
    """
    if prefix != '':
        try:
            #print("Trying to read: ", prefix)
            lines_read = 0
            while lines_read < lines:
                lines_read += 1
                data = ser.readline().decode("utf-8")
                if prefix in data:
                    try:
                        retdata = data
                        return retdata
                    except:
                        pass
            print("Read " + str(lines_read) + " lines without finding the data prefix.")
            return None
        except AttributeError as error:
            print("The serial connection is not initialized. Run the initSerial function first")
            print(error)
            return None
    else:
        try:
            return ser.readline().decode("utf-8")
                    
        except AttributeError as error:
            print("The serial connection is not initialized. Run the initSerial function first")
            print(error)
            return None
            

def sendData(ser, data):
    try:
        string = "<"
        for i in range(len(data)):
            if isinstance(data[i], float) or isinstance(data[i], np.float64):
                string += str(round(data[i], 3)) + ","
            else:
                string += str(data[i]) + ","
        string = string[:-1] + ">"
        string = string.encode('utf-8')
        ser.write(string)
        #print("The printed string sent to teensy: ", string)
        return True
    except serial.SerialTimeoutException as error:
        print("The write process timed out.")
        print(error)
        return False
 
    
def close(ser):
    ser.close()
