#!/usr/bin/python
import matplotlib.pyplot as plt


laser_nis = []
laser_nis_flag = []
radar_nis = []
radar_nis_flag = []

with open('nis.txt') as f:
    lines = f.readlines()
    for line in lines:
        data = line.split(':')
        if len(data) != 2:
            continue

        if data[0] == 'laser_nis':
            laser_nis.append(float(data[1].strip()))
            laser_nis_flag.append(7.8)
        if data[0] == 'radar_nis':
            radar_nis.append(float(data[1].strip()))
            radar_nis_flag.append(7.8)

plt.subplot(2,1,1)
plt.title('Laser NIS')
plt.plot(range(0, len(laser_nis)), laser_nis, 'r-')
plt.plot(range(0, len(laser_nis)), laser_nis_flag, 'r')
plt.subplot(2,1,2)
plt.title('Radar NIS')
plt.plot(range(0, len(radar_nis)), radar_nis, 'b-')
plt.plot(range(0, len(radar_nis)), radar_nis_flag, 'b')
plt.show()

