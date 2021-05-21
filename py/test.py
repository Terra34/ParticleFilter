import numpy as np
from time import time
from scipy.stats import norm
import csv
import matplotlib.pyplot as plt
import math
from numpy.random import uniform
import scipy


def get_roll_pitch(meas):
    roll_acc = np.arctan2(-meas[1], math.sqrt(meas[0]**2 + meas[2]**2))
    pitch_acc = np.arctan2(meas[0], math.sqrt(meas[1]**2 + meas[2]**2))
    return roll_acc, pitch_acc


def get_yaw_from_state(filter_state, meas):
    sina = np.sin(filter_state[0])
    sinb = np.sin(filter_state[1])
    cosa = np.cos(filter_state[0])
    XH = meas[0] * np.cos(filter_state[1]) + meas[1] * sinb * sina + meas[2] * sinb * cosa
    YH = meas[1] * cosa + meas[2] * sina
    return np.arctan2(-YH, XH)


def generate_particles(N):
    particles = np.empty((N,2))
    particles[:,0] = uniform(-np.pi, np.pi,size=N)
    particles[:,1] = uniform(-np.pi, np.pi,size=N)
    particles = particles * 180/np.pi
    return particles


def predict(particles, gyroData_t, deltaTime=0.1, std=10):
    particles[:,0] += gyroData_t[0]*deltaTime + std*np.random.randn(len(particles))
    particles[:,1] += gyroData_t[1]*deltaTime + std*np.random.randn(len(particles))
    # now we have to check if any of the particles got outside -180:180 interval
    # I assumed that gyroData*deltaTime will never be larger than 180 degrees
    for i in range(len(particles)):
        for j in range(2):
            if particles[i,j] > 180:
                particles[i,j] -= 360
            elif particles[i,j] < -180:
                particles[i,j] += 360
            else:
                continue


def update(particles, weights, state, std=10):
    weights *= 1 / (pow(particles[:,0] - state[0], 2) + pow(particles[:,1] - state[1], 2))
    weights += 1e-300


def estimate(particles, weights):
    mean = np.average(particles, weights=weights, axis=0)
    var = np.average((particles-mean)**2, weights=weights, axis=0)
    return mean,var


def resample(weights):
    K = N/sum(weights)
    wcp = 0
    j = 0
    indexes = np.zeros(N, 'i')
    cp = np.zeros(N)
    for i in range(N):
        wcp += weights[i]*K
        cp[i] = np.rint(wcp)
        fact = cp[i]-cp[i-1]
        for _ in range(int(fact)):
            indexes[j] = i
            j += 1

    return indexes


def resample_from_index(particles, weights, indexes):
    particles[:] = particles[indexes]
    weights.resize(len(particles),refcheck=False)
    weights.fill(1.0/len(weights))
    return particles


with open('imu_data.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)
    for counter, row in enumerate(reader):
        if (row != []):
            continue
print("Number of lines in this file: ", counter)

time = np.empty([int(counter/2), 1])
meas = np.empty([int(counter/2), 9])

with open('imu_data.csv', 'r') as file:
    reader = csv.reader(file)
    next(reader)
    for i, row in enumerate(reader):
        if (row != []):
            time[int((i-1)/2),0] = float(row[0])
            for j, data in enumerate(row[1:]):
                meas[int((i-1)/2), j] = float(data)

# acclData, gyroData, magData
acclData = meas[:,0:3]
gyroData = meas[:,3:6]
magData = meas[:,6:]

roll = np.empty([int(counter/2),1])
pitch = np.empty([int(counter/2),1])
for i, meas_t in enumerate(meas):
    roll[i,0], pitch[i,0] = get_roll_pitch(meas_t)

roll = roll * 180/np.pi
pitch = pitch * 180/np.pi


yaw = np.empty([int(counter/2),1])
for i, magData_t in enumerate(magData):
    yaw[i,0] = get_yaw_from_state([roll[i]*np.pi/180,pitch[i]*np.pi/180], magData_t)*180/np.pi

N = 512
particles = generate_particles(N)
weights = 1/N * np.ones(N)

roll_PF = np.zeros_like(roll)
pitch_PF = np.zeros_like(pitch)
yaw_PF = np.zeros_like(yaw)

for i, meas_t in enumerate(meas):
    predict(particles, meas_t[3:5])
    roll_t, pitch_t = get_roll_pitch(meas_t)
    roll_t *= 180/np.pi
    pitch_t *= 180/np.pi

    update(particles, weights, [roll_t, pitch_t])

    [roll_PF[i], pitch_PF[i]], _ = estimate(particles, weights)

    indexes = resample(weights)
    particles = resample_from_index(particles, weights, indexes)
    yaw_PF[i] = get_yaw_from_state([roll_PF[i]*np.pi/180, pitch_PF[i]*np.pi/180], meas_t[6:]) * 180/np.pi

out = np.zeros((98, 3), dtype=np.float32)
with open('arty_data.csv', 'r') as file:
    ff = file.readlines()
    for i in range(len(ff)):
        out[i,:] = ff[i].strip().split(',')

out2 = np.zeros((98, 3), dtype=np.float32)
with open('arty_data_lfsr.csv', 'r') as file:
    ff = file.readlines()
    for i in range(len(ff)):
        out2[i,:] = ff[i].strip().split(',')


plt.figure()
plt.subplot(311)
plt.title("Roll")
plt.plot(out[:, 0], label="Arty gauss")
plt.plot(roll_PF, label="Python")
plt.plot(out2[:, 0], label="Arty lfsr")
plt.legend()
plt.subplot(312)
plt.title("Pitch")
plt.plot(out[:, 1])
plt.plot(pitch_PF)
plt.plot(out2[:, 1])
plt.subplot(313)
plt.title("Yaw")
plt.plot(out[:, 2])
plt.plot(yaw_PF)
plt.plot(out2[:, 2])
plt.show()