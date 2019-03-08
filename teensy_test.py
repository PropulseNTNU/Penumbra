import teensy_interface as ti
import time

ser = ti.initSerial("COM12", 9600, 1)

while True:
    controlDump = ti.readControlSignal(ser, "c_s", 50)
    print(controlDump)
    time.sleep(1/2)
