import serial
import struct

data = []
received_data = []

def handleData(data: bytes):
    temp = data.decode().strip().split(",")
    received_data.append([ float(t) for t in temp ])
    


file = open('imu_data.csv', 'r')
file_data = file.readlines()
for line in file_data[1:]:
    if line != '\n':
        line = line.strip()
        temp = line.split(',')
        temp = [ float(t) for t in temp[1:] ]
        data.append(temp)
file.close()


ser = serial.Serial('COM6', 115200)

for d in data:
    # send all 36 bytes at once
    bytes_to_send = []
    for num in d:
        ba = bytearray(struct.pack("f", num)) # use "!f" to swap endianness
        bytes_to_send += [ byte for byte in ba ]

    ser.write(bytearray(bytes_to_send))

    # wait for response from arty
    while True:
        ser_data = ser.readline()
        
        if len(ser_data) >= 1:
            handleData(ser_data)
            break

ser.close()

file = open('arty_data_lfsr.csv', 'w')
for received in received_data:
    file.write(",".join([ str(r) for r in received ]) + "\n")

file.close()