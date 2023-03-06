###################################################################
#                                                                 #
#                 Produced by Luis Miguel Castillo Rápalo         #
#                 e-mail: luis.castillo@unah.hn                   #
#                           February 2023                         #
#                                                                 #
#                 Last Updated: 2 February23                      #
#                                                                 #
###################################################################


#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
#   -----  DESCRIPTION  -----
#   Name : Forward-Mole
#   Objetive: Correct DEM in urban areas according to a river shapefile in order to
#   improve hydrodynamic Modeling
#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

# --- requeried libraries --- ##
from osgeo import gdal
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin
import numpy as np
import glob
from matplotlib import pyplot
import pandas as pd
from osgeo import ogr, osr
import matplotlib.pyplot as plt


# --- Reading the required data --- #
path = "D:/Google_drive/Meu Drive/Papers/Paper - Sao_carlos/main_150.csv"  # Reading the Data
path2 = "D:/Google_drive/Meu Drive/Papers/Paper - Sao_carlos/"
channel = pd.read_csv(path)
channel = np.asarray(channel.sort_values(['lineid','row','newid'], ascending = [True, False, True]))  #Reordening the data acording LineID, Row and NewID
channel2 = channel.copy()
res = 30  # In meters
SRC = "EPSG:31983"  # Check your SRC

# --- Main loop --- #
loop, k, i, j = 0, 0, 0, 1

def hypo(i2, m):  # Calculated distance between two points
    return ((channel[i2, 8] - channel[m, 8]) ** 2 + (channel[i2, 7] - channel[m, 7]) ** 2) ** (1 / 2)

def forward(i, j, record):
    flag, flag_1, flag_k = 0, 0, 0
    while hypo(i, j) < np.sqrt(res**2 + (res*2)**2):  # Check if are cells neighborhoods
        if flag == 1:
            break
        if channel[i, 6] < channel[j, 6]:  # Check if downstream is higher than upstream cell
            i2 = j
            k, dist = [], [0]
            while i2 < len(channel):
                if flag_1 == 1:
                    flag = 1
                    break
                k.append(i2)
                if hypo(i2 - 1, i2) > np.sqrt(res**2 + (res*2)**2):  # Check if are not cells neighborhoods
                    if flag_k == 0:
                        i_f, j_f = i2, i2 + 1
                        m = i2  #
                        while m < len(channel):
                            if m == len(channel):
                                break
                            if (hypo(i2-1, m) <= np.sqrt(res**2 + (res*2)**2)) and (channel[m,4] != 1):  # Check if are cells neighborhoods
                                break
                            m += 1
                            if m == len(channel):
                                m = int(np.argwhere(channel[:, 0] > 1)[0])
                        i2 = m
                        k[-1] = i2
                        dist.append(dist[-1] + hypo(k[-2], k[-1]))
                        while (channel[i, 6] <= channel[i2, 6]) and (hypo(k[-2],k[-1])<= np.sqrt(res**2 + (res*2)**2)):  # Check if downstream is higher than upstream cell
                            i2 += 1
                            k.append(i2)
                            dist.append(dist[-1] + hypo(k[-2], k[-1]))
                        channel[k, 6] = channel[i, 6] - ((channel[i, 6] - channel[i2, 6]) * np.asarray(dist[1:])) / (dist[-1])
                        flag_k = 1
                        flag_1 = 1
                    else:
                        while hypo(k[-2], k[-1]) > np.sqrt(res**2 + (res*2)**2):
                            i_f, j_f = i2, i2 + 1
                            m = i2  #
                            while m < len(channel):
                                if m == len(channel):
                                    break
                                if (hypo(i2-1, m) <= np.sqrt(res**2 + (res*2)**2)) and (channel[m,4] != 1):  # Check if are cells neighborhoods
                                    break
                                m += 1
                                if m == len(channel):
                                    m = int(np.argwhere(channel[:, 0] > 1)[0])
                            i2 = m
                            k[-1] = i2
                            dist.append(dist[-1] + hypo(k[-2], k[-1]))
                            while (channel[i, 6] <= channel[i2, 6]) and (hypo(k[-2],k[-1])<= np.sqrt(res**2 + (res*2)**2)):  # Check if downstream is higher than upstream cell
                                i2 += 1
                                k.append(i2)
                                dist.append(dist[-1] + hypo(k[-2], k[-1]))
                            channel[k, 6] = channel[i, 6] - ((channel[i, 6] - channel[i2, 6]) * np.asarray(dist[1:])) / (dist[-1])
                            flag_k = 1
                            flag_1 = 1
                elif channel[i, 6] >= channel[i2, 6]:
                    i_f, j_f = i2, i2 + 1
                    dist.append(dist[-1] + hypo(k[-2], k[-1]))
                    channel[k, 6] = channel[i, 6] - ((channel[i, 6] - channel[i2, 6]) * np.asarray(dist[1:])) / (dist[-1])
                    i = i2
                    j = i2 + 1
                    break
                else:
                    dist.append(dist[-1] + hypo(k[-1]-1, k[-1]))
                    i2 += 1
                if i2 == len(channel):
                    i_f, j_f = i, j
                    flag = 1
                    break
        else:
            i += 1
            j += 1
            if j == len(channel):
                i -=1
                j -=1
                break
            if hypo(i, j) > np.sqrt(res**2 + (res*2)**2):  # Check if are neighborhoods
                i_f, j_f = i + 1, j + 1
    i_f, j_f = i + 1, j + 1
    return i_f, j_f

i, j, i_f, j_f, out = 0, 1, 0, 0, 0
while i < len(channel):
    if out == 1:
        break
    while j < len(channel):
        record = channel[i, 4]
        i, j = forward(i, j, record)
        if channel[i,4] == channel[-2,4]:
            out = 1
            break
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
scat1 = ax.scatter(channel[:, 7], channel[:, 8], channel[:, 6], marker='o', s=20, c='black', label="Corrected Data")
scat2 = ax.scatter(channel2[:, 7], channel2[:, 8], channel2[:, 6], marker='o', s=5, alpha=0.5, c='tab:orange', label="Original Data")
fig.legend(loc="center right")
fig.show()

# --- exporting the raster --- #
output = np.column_stack((channel[:, [7, 8, 6]])).T
df1 = pd.DataFrame(output)
df1.to_csv(path2+'main_channel_ajusted_150.csv', index=False)

df1 = df1.sort_values(by = [1, 0], ascending = [False, True])
df1.to_csv(path2+'channel_ajusted_150.xyz', index=False, header=None, sep=" ")

dem_1 = gdal.Translate(path2+"main_channel_ajusted_150.tif",
                       path2+"channel_ajusted_150.xyz",
                       outputSRS=SRC)
dem_1 = None


