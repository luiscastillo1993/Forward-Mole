###################################################################
#                                                                 #
#           Produced by Luis Miguel Castillo Rápalo               #
#                 e-mail: luis.castillo@unah.hn /                 #
#                      luis.caver17@gmail.com                     #
#                           February 2023                         #
#                                                                 #
#                      Last Updated: 5 May 2024                   #
#                                                                 #
###################################################################

#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
#   -----  DESCRIPTION  -----
#   Name : Forward-Mole
#   Objetive: Correct DEM in urban areas according to a river shapefile in order to
#   improve hydrodynamic connectivity in water bodies
#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

# --- required libraries --- ##
import os
import time
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString

import matplotlib
matplotlib.use('TkAgg',force=True)
import matplotlib.pyplot as plt

start = time.time()
#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
# FM for 1D (0) or 2D (1)
v_2D = 0  # is the FM will works in 2D version
plot_results = 0  # is the FM will plot the results in a 3D scatter on the browser
alpha = 1  # This aids to control height differences considered to scan the channel width
elements_ahead = 6  # This controls the number of rivers ahead considered to find a reference point
# --- Reading the required data --- #
path2 = "D:/Google_drive/Meu Drive/Papers/Paper - Forward Mole/GIS/RGS/FW_RGS/"
dem_name = "rgs_dem"
hydro_name = "drenagem2"
#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
# Reading the original DEM
with rasterio.open(path2+dem_name+'.tif') as src:
    dem_o = src.read(1)
    kwargs = src.meta
    bounds = src.bounds
    res = np.round(src.res[1])
    res2 = (res ** 2 + res ** 2) ** 0.5
    res_spacing = 0.001 * res
# dem_new = dem_o.copy()

# The new raster name
if v_2D == 0:
    new_name = dem_name + '_1D_FM' + '.tif'
else:
    new_name = dem_name + '_2D_FM' + '.tif'

# Load the shapefile
gdf = gpd.read_file(path2+hydro_name+'.shp')
line_bounds = gdf.total_bounds
main = np.full((0, 7), np.nan, dtype=np.float32)
lines = gdf.geometry
df = pd.DataFrame()

# Extracting the coordinates with their id element
for i in range(len(lines)):
    if lines[i] is None:
        continue
    df_temp = pd.concat([pd.DataFrame(np.array(np.linspace(i,i,len(lines[i].coords[:])))) , pd.DataFrame(np.array(lines[i].coords[:]))],axis=1)
    df = pd.concat([df, df_temp], axis=0)
df = np.array(df)

# if res > 1:
#     df[:, 1] = np.round(df[:, 1],3)
#     df[:, 2] = np.round(df[:, 2],3)
# else:
#     df[:, 1] = np.round(df[:, 1],5)
#     df[:, 2] = np.round(df[:, 2],5)

print('###### Calculating the strahler order ######')
df_temp = df.copy()
# plt.scatter(df[:,1],df[:,2])
lines_stralher = np.empty((0, 2))
stralher = 1
flag = 1
flag_outlet = 1
id_outlet = -1
while len(np.unique(df_temp[:, 0])) > 1:
    order = []
    if flag == 0:
        break

    plt.scatter(df_temp[:, 1], df_temp[:, 2])
    # plt.show()
    for i in np.unique(df_temp[:, 0]):

        if flag == 0:
            break
        # if the element is empty, continue to the next one
        if np.any(df_temp[:, 0] == i) == False:
            continue
        # if the outlet element was identified, continue to the next one
        if i == id_outlet:
            continue

        # looking in df_temp if the first or last pair or coords of the i element are present
        cond = df_temp[:, 1:] == df_temp[df_temp[:, 0] == i, 1:][0]
        cond2 = df_temp[:, 1:] == df_temp[df_temp[:, 0] == i, 1:][-1]
        cond_s = np.size(np.where((cond[:, 0] == True) & (cond[:, 1] == True)), 1)
        cond_s2 = np.size(np.where((cond2[:, 0] == True) & (cond2[:, 1] == True)), 1)
        # if they are present just once at least, means that i element is strahler 1
        if cond_s == 1 or cond_s2 == 1:
            # We find which cond_s is the one that does not have neighbors, then, we work with the other extreme because
            # is connected with other elements. Here, we catch if the i element is over
            # because it meets other elements with Strahler 1 and 2.
            if cond_s == 1:
                cond = cond2
            elif cond_s2 == 1:
                cond = cond
            cond_s = np.size(np.where((cond[:, 0] == True) & (cond[:, 1] == True)), 1)

            if cond_s > 2:
                order = np.append(order, int(i))
            elif cond_s == 2:
                # j is the new i to merge the next i element, we start over from df_temp
                temp_order = []
                j = i
                # temp_order = np.append(temp_order,j)
                while cond_s == 2:
                    # we cumulate the j values because they have the same stralher
                    temp_order = np.append(temp_order, j)
                    temp = df_temp[np.where((cond[:,0]==True) & (cond[:,1]==True)),0]
                    # k is the new j to merge the next j element
                    k = temp[np.where(temp != j)]
                    # to avoid the endless loop when the next element is inverted
                    if (df_temp[df_temp[:, 0] == k, 1:][-1] != df_temp[df_temp[:, 0] == j, 1:][-1]).all():
                        # remove the last one because it's a bug from the endless loop
                        cond = df_temp[:, 1:] == df_temp[df_temp[:, 0] == k, 1:][-1]
                    else:
                        cond = df_temp[:, 1:] == df_temp[df_temp[:, 0] == k, 1:][0]

                    cond_s = np.size(np.where((cond[:, 0] == True) & (cond[:, 1] == True)), 1)
                    j = k
                    if cond_s == 1:
                        # this is the outlet
                        flag = 0
                        break
                # to add the last one element detected
                temp_order = np.append(temp_order, k)
                # to delete repeated values that happened inside the while loop
                temp_order = np.unique(temp_order)
                order = np.append(order,temp_order)
    # Identify the outlet element, this to exclude it to be deleted in every iteration
    if flag_outlet == 1:
        idx = np.argwhere(df[:, 0] == order[:, np.newaxis])
        idx[:, 0] = df[np.any((df[:, 0] == order[:, np.newaxis]).T, axis=1), 0]
        x_y = np.array([np.round((df[idx[:, 1], 1] - bounds[0]) / res), np.round((bounds[3]-df[idx[:, 1], 2]) / res)]).astype(int).T
        z = np.array([dem_o[x_y[:, 1], x_y[:, 0]], idx[:, 0]])
        id_outlet = z[1, np.argwhere(z[0, :] == np.min(z[0, :]))[0]]
        flag_outlet = 0
        order = order[np.where(order != id_outlet)]
    # Dropping the identified lines
    idx = np.any(np.equal(df_temp[:,0],order[:,np.newaxis]),axis=0)
    df_temp = df_temp[np.invert(idx)]
    # For all the elements, we assign a stralher order
    lines_stralher = np.vstack([lines_stralher,np.stack((order,np.linspace(stralher,stralher,len(order))),axis=1)])
    stralher = stralher+1
    # plt.scatter(df_temp[:,1],df_temp[:,2])

print('###### Updating the shapefile with row_id and strahler order ######')
gdf["row_id"] = range(len(gdf))
gdf["strahler"] = np.empty((len(gdf["row_id"]),1))
for i in gdf["row_id"]:
    if len(lines_stralher[np.argwhere(i == lines_stralher[:, 0]), 1]) == 0:
        continue
    else:
        gdf.loc[i, "strahler"] = lines_stralher[np.argwhere(i == lines_stralher[:, 0]), 1][0][0]

df = np.append(df, np.zeros((len(df), 1)), axis=1)
for i in lines_stralher[:, 0]:
    df[np.argwhere(df[:, 0] == i), 3] = lines_stralher[np.argwhere(lines_stralher[:, 0] == i), 1]

print('###### Check if there are inverted elements. If yes, we correct them ######')
idx_pass = []
while len(np.unique(df[:,0])) > 1:
    # idx list to wipe up elements to simplify the routing on each loop
    idx = []
    for geom, strahler, row_id in zip(gdf.geometry, gdf['strahler'], gdf['row_id']):
        if geom and not geom.is_empty and geom.is_valid and geom.geom_type == 'LineString':
            # To no repeat already processed elements
            if np.any(idx_pass == row_id):
                continue
            # Extracting nodes info as how many neighbor are
            xy = np.array([geom.coords[0][0],geom.coords[0][1]])
            xy2 = np.array([geom.coords[-1][0],geom.coords[-1][1]])
            xy_c = len(np.argwhere((df[:,1:2]==xy[0]) & (df[:, 2:3]==xy[1]))[:, 0])
            xy2_c = len(np.argwhere((df[:,1:2]==xy2[0]) & (df[:, 2:3]==xy2[1]))[:, 0])
            # The outlet element is the only one that should be inverted
            if row_id == id_outlet:
                if (xy_c == 1 and np.any(df[np.argwhere(df[:,1:2]==xy2),3]>1)):
                    # Reverse geometry
                    x,y = gdf.iloc[row_id].geometry.xy
                    gdf.loc[row_id,'geometry'] = LineString(list(zip(np.flip(x),np.flip(y))))
                continue
            # We check if the element have one extreme alone, if yes, we check if inverted
            if np.any(xy_c == 1 or xy2_c == 1):
                # If the first element is alone and the other one is connected with neighbor with strahler higher than 1  so the direction it's ok
                if (strahler == 1):
                    if (xy_c == 1 and np.any(df[np.argwhere(df[:,1:2]==xy2),3]>1)):
                        idx = np.append(idx,row_id)
                    else:
                        # Reverse geometry
                        x, y = gdf.iloc[row_id].geometry.xy
                        gdf.loc[row_id, 'geometry'] = LineString(list(zip(np.flip(x), np.flip(y))))
                        idx = np.append(idx, row_id)
                # If the first element is alone and the other one have more than 1 neighbor so the direction it's ok
                elif (strahler > 1):
                    if (xy_c == 1) and (len(np.unique((np.argwhere(df[:,1:3]==xy2))[:,0]))>1):
                        idx = np.append(idx, row_id)
                    else:
                        # Reverse geometry
                        x, y = gdf.iloc[row_id].geometry.xy
                        gdf.loc[row_id, 'geometry'] = LineString(list(zip(np.flip(x), np.flip(y))))
                        idx = np.append(idx, row_id)
    idx_pass = np.append(idx_pass,idx)
    idx = np.any(np.equal(df[:,0], idx[:,np.newaxis]),axis=0)
    df = df[np.invert(idx)]

print('###### Iterate over the shapefile features ######')
for geom, strahler, row_id in zip(gdf.geometry, gdf['strahler'], gdf['row_id']):
    if geom and not geom.is_empty and geom.is_valid and geom.geom_type == 'LineString':
        # Convert line to points
        k = 1
        for i in range(len(geom.coords)-1):
            num_points = int(
                ((geom.coords[i][0] - geom.coords[i+1][0]) ** 2 + (geom.coords[i][1] - geom.coords[i+1][1]) ** 2) ** 0.5 / res_spacing) + 1
            coordinates_array = np.linspace(geom.coords[i], geom.coords[i+1], num_points)
            coordinates_array = np.array([np.round((coordinates_array[:,0]-(bounds[0]-res/2))/res)*res+bounds[0]-res/2,
                                 np.round((coordinates_array[:,1]-(bounds[1]-res/2))/res)*res+bounds[1]-res/2]).T
            coordinates_array = coordinates_array[np.sort(np.unique(coordinates_array,axis=0, return_index=True)[1])]
            # Dropping the last row to avoid repeat the first and last point between segments
            coordinates_array = coordinates_array[:-1]
            # Populate the matrix with line points and additional data

            col = np.array((coordinates_array[:, 0] - bounds[0]) / res, dtype=int)
            row = np.array((bounds[3] - coordinates_array[:, 1]) / res, dtype=int)
            temp = np.vstack([col, row, np.repeat(strahler, len(col)), np.repeat(row_id, len(col)),
                              coordinates_array[:, 0], coordinates_array[:, 1], range(k, k + len(col), 1)]).T
            main = np.vstack([main, temp])
            k = k + len(col)

print('###### Overall forward mole loop ######')
for stra in range(1, int(np.max(main[:,2]))+1):
    qq = 0
    if stra == np.max(main[:,2]):
        order = np.unique(main[np.argwhere(main[:, 2] == stra), 3]).astype(int)
        order = order[::-1]
    else:
        order = np.unique(main[np.argwhere(main[:, 2] == stra), 3]).astype(int)
    for j in order:
        qq += 1
        print(' ####################################################')
        print(' ##### Processing element from strahler order -' + str(stra) + '- #####')
        print(' ##### Processing element No. ' + str(j) + ' ### ' +
              str(len(np.unique(main[np.argwhere(main[:, 2] == stra), 3]).astype(int)) - qq) + ' Elements left #####')
        # Segments to analyze
        temp_main = main[np.argwhere(main[:, 3] == j)[:, 0], :]

        if len(temp_main) == 0:  # if some line_id does not exist, the code will jump the line
            continue

        fl, lines_forward = j, 0
        # Creating a 3x3 empty windows
        matrix = np.full((3, 3, 3), np.nan, dtype=np.float32)
        t2 = np.zeros((8, 2))

        while lines_forward < elements_ahead:  # number of lines to see forward
            # Taking the last pair of coordinates (row and col)
            x_coor = temp_main[-1, 0]
            y_coor = temp_main[-1, 1]

            # Window to identify neighboring coordinates
            t2[:, 0] = [x_coor-1, x_coor, x_coor+1, x_coor-1, x_coor+1, x_coor-1, x_coor, x_coor+1]
            t2[:, 1] = [y_coor-1, y_coor-1, y_coor-1, y_coor, y_coor, y_coor+1, y_coor+1, y_coor+1]
            # Finding the rows in main where required data is stored for our 3x3 window
            t3 = np.argwhere([(main[:, 0] == t2[value, 0]) & (main[:, 1] == t2[value, 1]) for value in range(len(t2))])[:, 1]
            # Filling the 3x3 window with rowd_id values
            matrix[0,(main[t3, 0] - min(t2[:, 0])).astype(int), (main[t3, 1] - min(t2[:, 1])).astype(int)] = main[t3, 3]
            # Filling the 3x3 window with k order line
            matrix[1,(main[t3, 0] - min(t2[:, 0])).astype(int), (main[t3, 1] - min(t2[:, 1])).astype(int)] = main[t3, 6]
            # Filling the 3x3 window with stralher line
            matrix[2,(main[t3, 0] - min(t2[:, 0])).astype(int), (main[t3, 1] - min(t2[:, 1])).astype(int)] = main[t3, 2]
            # Dropping not allowed neighbors
            # Dropping Nan values
            # Transforming t2 values of rows and columns for a window of 3x3 size
            t4 = t2.copy()
            t4[:,0], t4[:,1] = t4[:,0]-min(t4[:,0]), t4[:,1]-min(t4[:,1])
            t4 = t4.astype(int)
            condition = np.isnan(matrix[0, t4[:, 0], t4[:, 1]])

            # Windows of the row_id values
            condition2 = np.array([value in main[:,3] for value in matrix[0, t4[:, 0], t4[:, 1]]])

            # Dropping values that are lower than the last value K and belongs to the same row_id
            condition2[np.argwhere(condition2==True)[((matrix[1,t4[condition2][:,0],t4[condition2][:,1]]<temp_main[-1,6])&
                                                     (matrix[0,t4[condition2][:,0],t4[condition2][:,1]]==temp_main[-1,3]))]] = False

            # Dropping values that are lower than the last value K and with lower stralher order
            t4 = np.where(condition.reshape(-1, 1), np.nan, t4)
            t4 = np.where(condition2.reshape(-1, 1), t4, np.nan)
            t4 = t4[~np.isnan(t4).any(axis=1)].astype(int)

            # sorting t2 indexes
            t4 = t4[np.lexsort((matrix[1,t4[:,0],t4[:,1]], matrix[0,t4[:,0],t4[:,1]]))]

            # returning original values of t4 for row and columns as t2
            t4[:,0], t4[:,1] = t4[:,0]+min(t2[:,0]), t4[:,1]+min(t2[:,1])

            # if t4 is empty, means that could be the outlet because there is no more parts downstream of it
            if len(t4)==0:
                break
            # extracting all info from main for the filtered cells
            temp2 = main[np.argwhere([(main[:,0]==t4[value,0])&(main[:,1]==t4[value,1]) for value in range(len(t4))])[:,1],:]

            if len(temp2) == 0: # this check if temp 2 is empty, this happend when we are working with the outlet line (last one)
                break
            if np.any(temp2[:, 3] != temp2[0, 3]):
                # if there is row_id different, we sort according x[2] and x[6] which are stralher number and k order, respectively.
                if np.all(temp2[:,2]==temp2[0,2]):
                    # Finding the group with the largest K order, this means that is the end of one channel, we looking for
                    # the beginning of the downstream channel
                    df = pd.DataFrame(temp2)
                    # Ordening the pandas group by the mean mentioned before.
                    temp2 = pd.concat([df[df[3] == df.groupby(3)[6].mean().idxmax()],df[df[3] != df.groupby(3)[6].mean().idxmax()]]).to_numpy()
                else:
                    df = pd.DataFrame(temp2)
                    grouped_df = df.groupby(3)[6].mean()
                    sorted_group_indices = np.argsort(df.groupby(3)[6].mean().values)[::-1]
                    # Sort the DataFrame based on the sorted group indices
                    temp2 = df.merge(grouped_df, left_on=3, right_index=True).sort_values(by='6_y', ascending=False).to_numpy()[:, :-1]
            if temp2[-1, 2] >= temp_main[-1, 2]:
                # Gathering data from main according the last row in temp2, this is to connect with the next channel

                t = main[np.argwhere(main[:, 3] == temp2[-1, 3])[:, 0], :]
                if len(t) <= 2:
                    temp_main = np.vstack([temp_main, t])
                else:
                    try:
                        # Dropping all rows before our interest reference in t
                        rc = temp2[np.argwhere((temp2[:,2]==max(temp2[:,2]))&(temp2[:,6]==min(temp2[:,6])))[0],0:2]
                        # Finding the index in t and erasing rows
                        t = t[int(np.argwhere((t[:, 0] == rc[0, 0]) & (t[:, 1] == rc[0, 1]))[0]):,:]
                    except:
                        # Dropping all rows before our interest reference in t
                        rc = temp2[np.argwhere((temp2[:,2]==max(temp2[:, 2]))&(temp2[:, 6]==min(temp2[:,6])))[1],0:2]
                        # Finding the index in t and erasing rows
                        t = t[int(np.argwhere((t[:, 0] == rc[0, 0]) & (t[:, 1] == rc[0, 1]))[0]):, :]
                    # joining with the temp_main
                    temp_main = np.vstack([temp_main, t])
            else:
                temp_main = np.vstack([temp_main, temp2])
            lines_forward += 1
            fl = temp_main[-1,6].astype(int)

        # Calculate distances between points according to the coordinates
        temp = temp_main
        temp = np.hstack([temp, np.zeros((temp.shape[0], 1))])
        for i in range(1, len(temp)):
            temp[i, -1] = np.sqrt((temp[i, 4]-temp[i-1, 4])**2 + (temp[i, 5]-temp[i-1, 5])**2)

        flag_r = 0 # to indicate that right comes first
        flag_l = 0
        flag_2D_off = 0
        flag_2nd_loop_out = 0
        flag_2D_repeat = 0
        flag_backup = 0
        flag_out_element = 0
        flag_first = 0
        flag_once = 0
        k_backup = 0
        flag_empty = 0
        flag_comment = 0

        t = -1
        while t < 10000:
            flag_2nd_loop_out = 0

            if flag_out_element == 1:
                break

            if t == -1:  # the first iteration for the center line
                temp_w = temp.copy()
                if flag_empty == 1:
                    flag_2D_repeat = 0
                # if v_2D == 0:  # for the 1D case, just one loop
                #     t = 10000
            else:
                if flag_2D_repeat == 1 and flag_comment == 0:
                    print('...fixing the right side')
                    flag_comment = 1

                if t > np.size(width_right, 0)-1 and flag_r == 0:
                    print('...fixing the left side')
                    flag_r = 1
                    if flag_once == 0:
                        flag_first = 1
                        flag_once = 1  # only once per segment
                    t = 0  # now is left turn

                if flag_r == 0:
                    temp_w = temp.copy()
                    temp_w[:, 0:2] = temp_w[:, 0:2] + width_right[t, :]
                    if flag_once == 1:
                        flag_empty = 1
                    t = t + 1
                else:
                    temp_w = temp.copy()
                    if t > np.size(width_left, 0)-1:
                        temp_w = temp.copy()  #  here we return to the central line
                        flag_backup = 1
                        flag_r = 0
                        flag_2D_repeat = 0  # to activate the width scan
                        if flag_first == 1:
                            flag_first = 0
                        t = -1
                    else:
                        temp_w[:, 0:2] = temp_w[:, 0:2] + width_left[t, :]
                        t = t + 1

            # Main Forward Mole script
            i = j

            if flag_2D_repeat == 1 and flag_first == 1:
                k = 0
            elif flag_2D_repeat == 0 and flag_empty == 1:
                k = k_0
                flag_empty = 0
            elif flag_2D_repeat == 1 and flag_empty == 1 and flag_once == 1:
                k = k_0
                flag_empty = 1
            else:
                k = 0

            if flag_backup == 1:
                i = i_backup
                k = k_backup
                flag_backup = 0

            while i == j and flag_2nd_loop_out == 0:
                # Main Forward Mole Logic
                if np.size(temp_w, 0) == k + 1:  # to avoid index out of bounds
                    break

                # this if identify a vertical pit
                if(dem_o[int(temp_w[k,1]),int(temp_w[k,0])] < dem_o[int(temp_w[k+1, 1]),int(temp_w[k+1, 0])]):

                    if flag_2D_repeat < 1:
                        print('-------- Pit found at ' + str(k*res) + ' meters')

                    # if it is empty, means that there is no cells for solution and this only happens at the outlet
                    if len(np.argwhere((dem_o[int(temp_w[k, 1]), int(temp_w[k, 0])] >= dem_o[temp_w[:, 1].astype(int), temp_w[:, 0].astype(int)]) &
                                                   (dem_o[temp_w[:,1].astype(int),temp_w[:,0].astype(int)] < dem_o[int(temp_w[k,1]),int(temp_w[k,0])])))==0:

                        flag_out_element = 1
                        break

                    #######################################################################################################
                    # identify the cells that belongs to the error and solution by the FM logic
                    refs = np.argwhere((dem_o[int(temp_w[k, 1]), int(temp_w[k, 0])] > dem_o[temp_w[:, 1].astype(int), temp_w[:, 0].astype(int)]) &
                                                   (dem_o[temp_w[:, 1].astype(int),temp_w[:, 0].astype(int)] < dem_o[int(temp_w[k, 1]),int(temp_w[k, 0])]))
                    if len(refs)<=1:
                        k = k + 1
                        continue
                    elif len(np.argwhere(k < refs)) == 0:
                        k = k + 1
                        continue
                    refs = np.arange(k, refs[np.argwhere(k<refs)[0][0].astype(int)], 1)

                    # and temp_w[refs[-1], 3] != i
                    if t != -1:
                        if refs[-1] >= refs_backup * 1:
                            # if j ==5 and flag_r == 1 and t == 8 and k > 700:
                            #     print('t = ' + str(t) + 'and k = ' + str(k))
                            if len(np.argwhere(refs < refs_backup)) == 0:  # no solution
                                flag_2nd_loop_out = 1
                                continue
                            refs = refs[:np.argwhere(refs == refs_backup)[-1][0].astype(int)+1]
                            # flag_2nd_loop_out = 1
                            if (t > np.size(width_left, 0)-1) and flag_r == 1:
                                flag_r = 0
                                flag_backup = 1
                                t = -1  #  to get back to the central line

                            dem_o[int(temp_w[int(refs[-1]), 1]), int(temp_w[int(refs[-1]), 0])] = dem_backup

                    if len(refs) == 0:
                        k = k + 1
                        if k > refs_backup and t != -1:
                            flag_2nd_loop_out = 1
                        continue

                    # if flag_2D_repeat == 1 and refs[-1] >= refs_backup*1.1:
                    #     flag_2nd_loop_out = 1
                    #     flag_empty = 1
                    #     continue

                    # for the identified cells, we apply linear interpolation.
                    dem_o[temp_w[refs[1:-1],1].astype(int),temp_w[refs[1:-1],0].astype(int)] = \
                        ((dem_o[int(temp_w[refs[0],1]),int(temp_w[refs[0],0])] - dem_o[int(temp_w[refs[-1], 1]),int(temp_w[refs[-1], 0])]) *
                         (np.cumsum(temp_w[refs[1:-1], 7]) - np.sum(temp_w[refs[1:], 7]))) / (
                                       temp_w[refs[0], 7] - np.sum(temp_w[refs[:], 7])) + dem_o[int(temp_w[refs[-1], 1]), int(temp_w[refs[-1], 0])]

                    #######################################################################################################

                    if j == 39 and flag_once == 0:
                        x_ = np.array((0))
                        y_ = np.array((0))
                        # plt.scatter(temp_w[refs, 0], temp_w[refs, 1])
                        # plt.show()
                        # plt.close()
                        to_plot = dem_o
                        to_plot[dem_o < 0] = np.nan
                    if j == 39 and j < 41:
                        if k >0:
                            plt.scatter(temp_w[refs, 0], temp_w[refs, 1], s=2)
                            # plt.show()
                    # if j >41:
                    #     zzz = 1
                    #     plt.imshow(to_plot, cmap='gray')
                    #     plt.scatter(x_, y_, s=2)

                    ######################################################################################################
                    if v_2D == 1 and flag_2D_repeat == 0:
                        # secs = np.array((k_backup, int(np.ceil(np.average((k_backup, k)))), k))
                        secs = np.array((int(np.ceil(np.average((k_backup, k)))), int(np.ceil(np.average((k_backup, k))))))

                        # dif = dem_o[int(temp_w[k + 1, 1]), int(temp_w[k + 1, 0])] - dem_o[int(temp_w[k, 1]), int(temp_w[k, 0])]
                        # dif = np.average(dem_o[(temp_w[refs, 1]).astype(int),(temp_w[refs, 0]).astype(int)]) - dem_o[int(temp_w[k, 1]), int(temp_w[k, 0])]
                        # dif = dif*0.5

                        dif = np.max(dem_o[(temp_w[refs, 1]).astype(int),(temp_w[refs, 0]).astype(int)]) - dem_o[int(temp_w[k, 1]), int(temp_w[k, 0])]
                        delta = np.array((dif*alpha, dif*alpha))
                        wa_l, wa_r, wa_u, wa_d = [], [], [], []
                        for jj in secs:
                            # this makes a vector through the dem from the points as reference, it takes the distance
                            # from the reference to the no condition of height
                            for ii in delta:
                                dem_ref = dem_o[int(temp_w[jj, 1]), int(temp_w[jj, 0])] + ii

                                L = dem_ref > dem_o[int(temp_w[jj, 1]), range(int(temp_w[jj, 0] - 1), -1, -1)]
                                if L.all() == True:
                                    L = -1
                                else:
                                    L = (np.where(~L)[0][0] + 1) * res
                                R = dem_ref > dem_o[int(temp_w[jj, 1]), range(int(temp_w[jj, 0] + 1), np.shape(dem_o)[1])]
                                if R.all() == True:
                                    R = -1
                                else:
                                    R = (np.where(~R)[0][0] + 1) * res
                                U = dem_ref > dem_o[range(int(temp_w[jj, 1] - 1), -1, -1), int(temp_w[jj, 0])]
                                if U.all() == True:
                                    U = -1
                                else:
                                    U = (np.where(~U)[0][0] + 1) * res
                                D = dem_ref > dem_o[range(int(temp_w[jj, 1] + 1), np.shape(dem_o)[0]), int(temp_w[jj, 0])]
                                if D.all() == True:
                                    D = -1
                                else:
                                    D = (np.where(~D)[0][0] + 1) * res
                                # for diagonal, we need to build the coordinates separately
                                temp_l = np.arange(int(temp_w[jj, 0]), -1, -1)
                                temp_u = np.arange(int(temp_w[jj, 1]), -1, -1)
                                temp_r = np.arange(int(temp_w[jj, 0]), np.shape(dem_o)[1])
                                temp_d = np.arange(int(temp_w[jj, 1]), np.shape(dem_o)[0])
                                if len(temp_u) >= len(temp_l):
                                    temp_lu = np.column_stack((temp_u[0:len(temp_l)], temp_l))
                                else:
                                    temp_lu = np.column_stack((temp_u, temp_l[0:len(temp_u)]))
                                if len(temp_d) >= len(temp_r):
                                    temp_rd = np.column_stack((temp_d[0:len(temp_r)], temp_r))
                                else:
                                    temp_rd = np.column_stack((temp_d, temp_r[0:len(temp_d)]))
                                if len(temp_d) >= len(temp_l):
                                    temp_ld = np.column_stack((temp_d[0:len(temp_l)], temp_l))
                                else:
                                    temp_ld = np.column_stack((temp_d, temp_l[0:len(temp_d)]))
                                if len(temp_u) >= len(temp_r):
                                    temp_ru = np.column_stack((temp_u[0:len(temp_r)], temp_r))
                                else:
                                    temp_ru = np.column_stack((temp_u, temp_r[0:len(temp_u)]))

                                RU = dem_ref >= dem_o[temp_ru[:, 0], temp_ru[:, 1]]
                                if RU.all() == True:
                                    RU = -1
                                else:
                                    RU = (np.where(~RU)[0][0]-0)*(np.sqrt(res**2 + res**2))
                                LD = dem_ref >= dem_o[temp_ld[:, 0], temp_ld[:, 1]]
                                if LD.all() == True:
                                    LD = -1
                                else:
                                    LD = (np.where(~LD)[0][0]-0)*(np.sqrt(res**2 + res**2))
                                RD = dem_ref >= dem_o[temp_rd[:, 0], temp_rd[:, 1]]
                                if RD.all() == True:
                                    RD = -1
                                else:
                                    RD = (np.where(~RD)[0][0]-0)*(np.sqrt(res**2 + res**2))
                                LU = dem_ref >= dem_o[temp_lu[:, 0], temp_lu[:, 1]]
                                if LU.all() == True:
                                    LU = -1
                                else:
                                    LU = (np.where(~LU)[0][0]-0)*(np.sqrt(res**2 + res**2))

                                final_vector = np.array((L, LU, U, RU, R, RD, D, LD))
                                final_vector[final_vector < 0] = 0

                                # if ii == dif*0.1:
                                #     final_vector_10 = np.array((L, LU, U, RU, R, RD, D, LD))
                                #     final_vector_10[final_vector_10 < 0] = 0
                                # elif ii == dif / 2:
                                #     final_vector_50 = np.array((L, LU, U, RU, R, RD, D, LD))
                                #     final_vector_50[final_vector_50 < 0] = 0
                                # elif ii == dif:
                                #     final_vector_100 = np.array((L, LU, U, RU, R, RD, D, LD))
                                #     final_vector_100[final_vector_100 < 0] = 0

                            # final_vector = np.mean((final_vector_10, final_vector_50, final_vector_100), axis=0)

                            final_sum_vector = np.array((final_vector[3] + final_vector[7],
                                                         final_vector[1] + final_vector[5],
                                                         final_vector[0] + final_vector[4],
                                                         final_vector[2] + final_vector[6]))
                            final_matrix = np.array(((final_vector[7], final_vector[3]),
                                                     (final_vector[1], final_vector[5]),
                                                     (final_vector[0], final_vector[4]),
                                                     (final_vector[2], final_vector[6])))

                            # to avoid not undetermined fraction n/0
                            final_sum_vector[final_sum_vector == 0] = 0.1

                            ratio_1 = np.array(((final_sum_vector[2]) / (final_sum_vector[0]),
                                                (final_sum_vector[3]) / (final_sum_vector[0])))
                            ratio_11 = (final_sum_vector[1]) / (final_sum_vector[0])
                            ratio_2 = np.array(((final_sum_vector[2]) / (final_sum_vector[1]),
                                                (final_sum_vector[3]) / (final_sum_vector[1])))
                            ratio_22 = (final_sum_vector[0]) / (final_sum_vector[1])
                            ratio_3 = np.array(((final_sum_vector[0]) / (final_sum_vector[2]),
                                                (final_sum_vector[1]) / (final_sum_vector[2])))
                            ratio_33 = (final_sum_vector[3]) / (final_sum_vector[2])
                            ratio_4 = np.array(((final_sum_vector[0]) / (final_sum_vector[3]),
                                                (final_sum_vector[1]) / (final_sum_vector[3])))
                            ratio_44 = (final_sum_vector[2]) / (final_sum_vector[3])

                            # calculating ratios between diagonals and orthogonal lines
                            ratios = np.array((ratio_1, ratio_2, ratio_3, ratio_4))
                            ratios_2 = np.array((ratio_11, ratio_22, ratio_33, ratio_44))

                            index = np.argmin((final_vector[3] + final_vector[7], final_vector[1] + final_vector[5],
                                               final_vector[0] + final_vector[4], final_vector[2] + final_vector[6]))
                            # we compare if it is the minor element
                            if sum((ratios[index, :] < ratios_2[index]).astype(int)) == 2:
                                index_2 = np.argmin(ratios[index]) + int(index > 1) * (-2) + 2
                            elif sum((ratios[index, :] > ratios_2[index]).astype(int)) != 2:
                                index = index + int(index == 0) + index * (-1) + int(index > 1) * 3 + int(index > 2) * (-1)
                                index_2 = np.argmin(ratios[index]) + int(index > 1) * (-2) + 2
                            else:
                                index_2 = np.argmin(ratios[index]) + int(index > 1) * (-2) + 2
                            # we calculate the coordinates from the average points of width
                            if index < 2:
                                w_l = (np.sqrt(2) / 2) * final_matrix[index, 0]
                                w_r = (np.sqrt(2) / 2) * final_matrix[index, 1]
                                w_u = (np.sqrt(2) / 2) * final_matrix[index, 0]
                                w_d = (np.sqrt(2) / 2) * final_matrix[index, 1]
                                if index_2 == 2:
                                    w_l2 = final_matrix[index_2, 0]
                                    w_r2 = final_matrix[index_2, 1]
                                    w_u2, w_d2 = 0, 0
                                else:
                                    w_u2 = final_matrix[index_2, 0]
                                    w_d2 = final_matrix[index_2, 1]
                                    w_l2, w_r2 = 0, 0
                            else:
                                w_l = (np.sqrt(2) / 2) * final_matrix[index_2, 0]
                                w_r = (np.sqrt(2) / 2) * final_matrix[index_2, 1]
                                w_u = (np.sqrt(2) / 2) * final_matrix[index_2, 0]
                                w_d = (np.sqrt(2) / 2) * final_matrix[index_2, 1]
                                if index == 2:
                                    w_l2 = final_matrix[index, 0]
                                    w_r2 = final_matrix[index, 1]
                                    w_u2, w_d2 = 0, 0
                                else:
                                    w_u2 = final_matrix[index, 0]
                                    w_d2 = final_matrix[index, 1]
                                    w_l2, w_r2 = 0, 0

                            # two points which could be in any of the 4 quarters
                            w_l, w_r = -np.ceil(np.average((w_l, w_l2) / res)), np.ceil(np.average((w_r, w_r2) / res))
                            w_u, w_d = -np.ceil(np.average((w_u, w_u2) / res)), np.ceil(np.average((w_d, w_d2) / res))
                            wa_l.append(w_l), wa_r.append(w_r), wa_u.append(w_u), wa_d.append(w_d)

                        wa_l, wa_r, wa_u, wa_d = np.average(wa_l), np.average(wa_r), np.average(wa_u), np.average(wa_d)
                        w_l, w_r, w_u, w_d = np.abs(wa_l), np.abs(wa_r), np.abs(wa_u), np.abs(wa_d)
                        # two vectors of coordinates are made to represent the movement
                        if (index == 0) or (index_2 == 0):
                            if w_l == 0 and w_d == 0:
                                width_left = np.array(((0, 0), (0, 0)))
                            elif w_l == 0:
                                width_left = np.array(
                                    (np.ceil(np.arange(1, w_d + 1, res)), np.zeros(np.ceil(w_d / res).astype(int)))).astype(
                                    int)
                            elif w_d == 0:
                                width_left = np.array((np.zeros(np.ceil(w_l / res).astype(int)),
                                                       np.ceil(-1 * np.arange(1, w_l + 1, res)))).astype(int)
                            else:
                                if w_l < w_d:
                                    width_left = np.array((np.ceil(np.arange(1, w_d + 1, res)),
                                                       np.ceil(-1 * np.arange(1, w_d + 1, res) * (w_l / w_d)))).astype(int)
                                else:
                                    width_left = np.array((np.ceil(np.arange(1, w_l + 1, res) * (w_d / w_l)),
                                                          np.ceil(-1 * np.arange(1, w_l + 1, res)))).astype(int)

                            if w_r == 0 and w_u == 0:
                                width_right = np.array(((0, 0), (0, 0)))
                            elif w_r == 0:
                                width_right = np.array((-1 * np.ceil(np.arange(1, w_u + 1, res)),
                                                        np.zeros(np.ceil(w_u / res).astype(int)))).astype(int)
                            elif w_u == 0:
                                width_right = np.array(
                                    (np.zeros(np.ceil(w_r / res).astype(int)), np.ceil(np.arange(1, w_r + 1, res)))).astype(
                                    int)
                            else:
                                if w_r < w_u:
                                    width_right = np.array((-1 * np.ceil(np.arange(1, w_u + 1, res)),
                                                        np.ceil(np.arange(1, w_u + 1, res) * (w_r / w_u)))).astype(int)
                                else:
                                    width_right = np.array((-1 * np.ceil(np.arange(1, w_r + 1, res) * (w_u / w_r)),
                                                            np.ceil(np.arange(1, w_r + 1, res)))).astype(int)

                        elif (index == 1) or (index_2 == 1):
                            if w_l == 0 and w_u == 0:
                                width_left = np.array(((0, 0), (0, 0)))
                            elif w_l == 0:
                                width_left = np.array((-1 * np.ceil(np.arange(1, w_u + 1, res)),
                                                       np.zeros(np.ceil(w_u / res).astype(int)))).astype(int)
                            elif w_u == 0:
                                width_left = np.array((np.zeros(np.ceil(w_l / res).astype(int)),
                                                       np.ceil(-1 * np.arange(1, w_l + 1, res)))).astype(int)
                            else:
                                if w_l < w_u:
                                    width_left = np.array((-1 * np.ceil(np.arange(1, w_u + 1, res)),
                                                       np.ceil(-1 * np.arange(1, w_u + 1, res) * (w_l / w_u)))).astype(int)
                                else:
                                    width_left = np.array((-1 * np.ceil(np.arange(1, w_l + 1, res)* (w_u / w_l)),
                                                       np.ceil(-1 * np.arange(1, w_l + 1, res)))).astype(int)

                            if w_r == 0 and w_d == 0:
                                width_right = np.array(((0, 0), (0, 0)))
                            elif w_r == 0:
                                width_right = np.array(
                                    (np.ceil(np.arange(1, w_d + 1, res)), np.zeros(np.ceil(w_d / res).astype(int)))).astype(
                                    int)
                            elif w_d == 0:
                                width_right = np.array(
                                    (np.zeros(np.ceil(w_r / res).astype(int)), np.ceil(np.arange(1, w_r + 1, res)))).astype(
                                    int)
                            else:
                                if w_r < w_d:
                                    width_right = np.array((np.ceil(np.arange(1, w_d + 1, res)),
                                                        np.ceil(np.arange(1, w_d + 1, res) * (w_r / w_d)))).astype(int)
                                else:
                                    width_right = np.array((np.ceil(np.arange(1, w_r + 1, res)*(w_d / w_r)),
                                                        np.ceil(np.arange(1, w_r + 1, res)))).astype(int)

                        if np.size(width_right, 0) != 1:
                            width_right = width_right.T
                        if np.size(width_right, 1) != 2:
                            width_right = width_right.T
                        width_right[:, [0, 1]] = width_right[:, [1, 0]]
                        if np.size(width_left, 0) != 1:
                            width_left = width_left.T
                        if np.size(width_left, 1) != 2:
                            width_left = width_left.T
                        width_left[:, [0, 1]] = width_left[:, [1, 0]]

                        if len(refs)*res < np.sqrt(np.max(abs(width_right[:, 0]))**2 + np.max(abs(width_right[:, 1]))**2)*res:
                            width_right = width_right[0:len(refs), :]
                        if len(refs)*res < np.sqrt(np.max(abs(width_left[:, 0]))**2 + np.max(abs(width_left[:, 1]))**2)*res:
                            width_left = width_left[0:len(refs), :]


                        flag_2D_off = 1
                        flag_2D_repeat = 1

                    ######################################################################################################

                    if flag_2D_off == 1:

                        k = refs[-1]
                        k_0 = refs[0]
                        # i = temp_w[k, 3]
                        flag_2nd_loop_out = 1   # to get out from the center line and work the river width
                        flag_2D_off = 0  # to avoid lose k and i backup in the next iteration when river width is fixed
                        flag_comment = 0
                        i_backup = i
                        k_backup = k
                        dem_backup = dem_o[int(temp_w[k, 1]),int(temp_w[k, 0])]
                        refs_backup = refs[-1]
                        t = 0  # to get out from the center line
                        if flag_once == 0:
                            flag_first = 1
                    else:
                        if flag_2D_repeat == 1 and refs[-1] >= refs_backup:
                            flag_2nd_loop_out = 1
                            flag_empty = 1
                            if flag_r == 1 and flag_l == 1:
                                t = -1  # to return t the central line
                                flag_empty = 1
                                flag_r = 0
                        else:
                            k = refs[-1]
                            i = temp_w[k, 3]
                else:
                    if flag_2D_repeat == 1 and k >= refs_backup:
                        flag_2nd_loop_out = 1
                    else:
                        k += 1
                        i = temp_w[k, 3]

                if flag_r == 1 and t > np.size(width_left, 0)-1:
                    i = temp_w[k_backup, 3]
                    flag_first = 0

                if i != j:
                    flag_out_element = 1


# saving the new raster
with rasterio.open(os.path.join(path2, new_name), 'w', **kwargs) as dst:
    dst.write_band(1, dem_o.astype(rasterio.float32))

end = time.time()
print('Modified DEM successfully exported')
print('Time elapsed: ' + str((end-start)/60) + ' minutes')

# # to plot in plotly to compare results with original points
# df = {'x':main[:, 4], 'y': main[:, 5], 'z': dem_o[main[:, 1].astype(int), main[:, 0].astype(int)],'color':str('Original Points'),'size':10}
# df2 = {'x':main[:, 4], 'y': main[:, 5], 'z': dem_new[main[:, 1].astype(int), main[:, 0].astype(int)],'color':str('New Points'),'size':1}
# df = pd.DataFrame(df)
# df2 = pd.DataFrame(df2)
# df3 = pd.concat([df,df2])
#
# if plot_results == 1:
#     # import plotly.graph_objects as go
#     import plotly.io as pio
#     import plotly.express as px
#     # from scipy.interpolate import griddata
#     pio.renderers.default = 'browser'
#
#     fig = px.scatter_3d(df3, x = 'x', y='y', z='z',color='color',size='size',color_discrete_sequence=px.colors.qualitative.D3)
#     # Remove outline or border from points
#     fig.update_traces(marker=dict(line=dict(width=0)))
#     fig.show()
