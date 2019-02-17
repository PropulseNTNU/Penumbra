import serial

ser = serial.Serial('COM11')
for i in range(10000):
    print(ser.readline().decode("utf-8").strip().split(','))
ser.close()
