import datetime as dt
import xarray as xr
from netCDF4 import Dataset
import numpy as np 
import time 
import pandas as pd 
from tqdm import tqdm
from datetime import datetime, timedelta
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from shapely.geometry import mapping, box
import rioxarray
import fsspec, re, aiohttp, requests
from tqdm import tqdm
import glob

# def contar_archivos_npy(directorio):
#     archivos_npy = glob.glob(f'{directorio}/*.npy', recursive=True)
#     return len(archivos_npy)

# year = '2023'

# for m in range(1, 12):
#     month = str(m+1).zfill(2)

#     directorio = 'E:/Heladas/matrices_tablas/{}/{}/'.format(year, month)
#     total_archivos_npy = contar_archivos_npy(directorio)
#     # mat_todo = np.zeros((58,41,total_archivos_npy))
#     mat_todo = []
#     dates_todo = []
#     index_count = 0
#     for j in range(0,  total_archivos_npy):
#         day = str(j+1).zfill(2)
#         array_t = np.load('E:/Heladas/matrices_tablas/{}/{}/mat_1dia_{}.npy'.format(year, month, day))
#         df_t = pd.read_csv('E:/Heladas/matrices_tablas/{}/{}/mat_1dia_info_{}.csv'.format(year, month, day), header=None)
#         tabla_t = df_t[0].astype(str).values.tolist()

#         for i in range(0,  array_t.shape[-1]):
#             mat_todo.append(array_t[:,:,i])
#             dates_todo.append(str(tabla_t[i]))
#             pass
#         pass

#     mat_todo_dim = np.zeros((58,41,len(mat_todo)))

#     for i in range(0, len(mat_todo)):
#         mat_todo_dim[:,:,i] = mat_todo[i]

#     heladas = mat_todo_dim <= 2.3
#     heladas_01 = np.where(heladas, 1, 0)

#     mat_sum = np.zeros((heladas_01.shape[0],heladas_01.shape[1]))

#     tabla = pd.DataFrame(columns=["Fila", "Columna", "Intensidad", "Fecha_inicial", "Fecha_final"])


#     for fil in tqdm(range(0, heladas_01.shape[0])):
#         for col in range(0, heladas_01.shape[1]):
#             for dim in range(0, heladas_01.shape[-1]):

#                 if heladas_01[fil, col, dim] == 1:
#                     mat_sum[fil, col] += 1

#                     if mat_sum[fil, col] >= 2:

#                         if mat_sum[fil, col] == 2:
#                             ind_fc = tabla.index[(tabla['Fila'] == fil) & (tabla['Columna'] == col)]

#                             if ind_fc.empty or tabla.loc[ind_fc, 'Fecha_inicial'].values[0] != dates_todo[dim-1]:
#                                 fila_nueva = [fil, col, 2, dates_todo[dim-1], dates_todo[dim]]
#                                 tabla.loc[len(tabla)] = fila_nueva
#                             else:
#                                 tabla.loc[ind_fc, 'Intensidad'] =  tabla.loc[ind_fc, 'Intensidad'].values[0]+1
#                                 tabla.loc[ind_fc, 'Fecha_final'] = dates_todo[dim]
                        
#                         else:
#                             ind_fc = tabla.index[(tabla['Fila'] == fil) & (tabla['Columna'] == col) & (tabla['Fecha_inicial'] == dates_todo[int(dim-(mat_sum[fil, col]-1))])]
#                             tabla.loc[ind_fc, 'Intensidad'] =  tabla.loc[ind_fc, 'Intensidad'].values[0]+1
#                             tabla.loc[ind_fc, 'Fecha_final'] = dates_todo[dim]
#                     else:
#                         pass
#                 else:
#                     mat_sum[fil, col] = 0

#     print(tabla)

#     tabla.to_pickle('E:/Heladas/matrices_tablas/tabla_heladas/tabla_{}_{}.pkl'.format(month, year))




#################################################################################

# year = '2023'
# month = '01'
# tabla = pd.read_pickle('E:/Heladas/matrices_tablas/tabla_heladas/tabla_{}_{}.pkl'.format(month, year))
# mat_intensidad_max = np.zeros((58,41))
# max_intensity_index = tabla.groupby(['Fila', 'Columna'])['Intensidad'].mean()
# max_intensity_df = max_intensity_index.reset_index()
# max_intensity_df.columns = ['Fila', 'Columna', 'Intensidad (máxima)']

# for i in range(0, len(max_intensity_df)):
#     fil = max_intensity_df.loc[i]['Fila']
#     col = max_intensity_df.loc[i]['Columna']
#     intens = max_intensity_df.loc[i]['Intensidad (máxima)']

#     if intens <= 6:
#         mat_intensidad_max[fil][col] = 1
#     elif intens <= 12:
#         mat_intensidad_max[fil][col] = 2
#     elif intens <= 24:
#         mat_intensidad_max[fil][col] = 3
#     elif intens <= 48:
#         mat_intensidad_max[fil][col] = 4
#     elif intens <= 96:
#         mat_intensidad_max[fil][col] = 5
#     else:
#         mat_intensidad_max[fil][col] = 6


########################################################################


def contar_archivos_npy(directorio):
    archivos_npy = glob.glob(f'{directorio}/*.npy', recursive=True)
    return len(archivos_npy)

year = '2023'

for m in range(1, 12):
    month = str(m+1).zfill(2)

    directorio = 'E:/Heladas/matrices_tablas/{}/{}/'.format(year, month)
    total_archivos_npy = contar_archivos_npy(directorio)
    # mat_todo = np.zeros((58,41,total_archivos_npy))
    mat_todo = []
    dates_todo = []
    for j in range(0,  total_archivos_npy):
        day = str(j+1).zfill(2)
        array_t = np.load('E:/Heladas/matrices_tablas/{}/{}/mat_1dia_{}.npy'.format(year, month, day))
        df_t = pd.read_csv('E:/Heladas/matrices_tablas/{}/{}/mat_1dia_info_{}.csv'.format(year, month, day), header=None)
        tabla_t = df_t[0].astype(str).values.tolist()

        for i in range(0,  array_t.shape[-1]):
            mat_todo.append(array_t[:,:,i])
            dates_todo.append(str(tabla_t[i]))
            pass
        pass

mat_todo_dim = np.zeros((58,41,len(mat_todo)))

for i in range(0, len(mat_todo)):
    mat_todo_dim[:,:,i] = mat_todo[i]

heladas = mat_todo_dim <= 2.3
heladas_01 = np.where(heladas, 1, 0)

mat_sum_heladas_tot = np.zeros((heladas_01.shape[0],heladas_01.shape[1]))

for fil in tqdm(range(0, heladas_01.shape[0])):
    for col in range(0, heladas_01.shape[1]):
        for dim in range(0, heladas_01.shape[-1]):

            if heladas_01[fil, col, dim] == 1:
                mat_sum_heladas_tot[fil, col] += 1

            else:
                pass


np.save('E:/Heladas/matrices_tablas/{}/mat_veces_heladas_{}.npy'.format(year, year), mat_sum_heladas_tot)