import iotbx.xplor.map
from iotbx import ccp4_map
#ccp4 = ccp4_map.map_reader(file_name=map)
import numpy as np
map_name = "/Users/coquelle/IBS_2022/RS_SwissFEL/Fext_grad/FoFo_dark/1ps_7uJ/Fextr_grad_MC_mqFoFo.map"
xplor = iotbx.xplor.map.reader(file_name=map_name)


def map_data(map_obj):
    m_data = map_obj.data.as_double()
    n_real = map_obj.unit_cell_grid
    if(n_real == m_data.all()):
        return map_obj.data.as_numpy_array() / map_obj.header_rms
    else:
      # XXX hideously SLOW! MOVE TO C++
      # map_new = flex.double(flex.grid(n_real), 0)
      map_new = np.empty(n_real)
      o = m_data.origin()
      f = m_data.focus()
      for i in range(o[0],f[0]):
        for j in range(o[1],f[1]):
          for k in range(o[2],f[2]):
            map_new[i%n_real[0], j%n_real[1], k%n_real[2]] = m_data[i, j, k] / map_obj.header_rms
      return map_new


grid = np.array(xplor.gridding.n, dtype=np.float32)
first = np.array(xplor.gridding.first, dtype=np.float32)
unit_cell = np.array(xplor.unit_cell.parameters(), dtype=np.float32)
data = xplor.data.as_numpy_array()

#grid = np.array(ccp4.unit_cell_grid, dtype=np.float32)
#first = np.array(ccp4.data.as_double().origin(), dtype=np.float32)
#unit_cell = np.array(ccp4.unit_cell_parameters, dtype=np.float32)
#data = map_data(ccp4)

print(grid)
print(first)
print(unit_cell)
print(data.shape)


print(data.mean())

from matplotlib import pyplot as plt

#plt.hist(data.flatten(), "auto") #, weights=counts)
#plt.show()


import scipy.stats
data = data.flatten()
#_, p = scipy.stats.normaltest(data.flatten())
#print(_)
#print(p)

z_scores = scipy.stats.zscore(data)

print(z_scores.shape)
print(z_scores.max())
print(data.size)


##################################
xplor_1 = iotbx.xplor.map.reader(file_name=map_name_1) #map_name_1 is FoFo
xplor_2 = iotbx.xplor.map.reader(file_name=map_name_2) #map_name_2 is Fextr-Fc
grid_1 = np.array(xplor_1.gridding.n, dtype=np.float32)
grid_2 = np.array(xplor_2.gridding.n, dtype=np.float32)
first_1 = np.array(xplor_1.gridding.first, dtype=np.float32)
first_2 = np.array(xplor_2.gridding.first, dtype=np.float32)
unit_cell_1 = np.array(xplor_1.unit_cell.parameters(), dtype=np.float32)
unit_cell_2 = np.array(xplor_2.unit_cell.parameters(), dtype=np.float32)
data_1 = xplor_1.data.as_numpy_array()
data_2 = xplor_2.data.as_numpy_array()

#Since maps are sigma-scaled, the value of the map and the z-score should be very similar.
#Therefore we can work directly with maps instead of mapping the z-score on the maps
#generate the mask
mask = np.where(np.abs(data_1) >= 3, 1, 0)

#Write mask to pickle as to avoid opening the FoFo map the whole time
out=open('map_mask.pickle' ,'wb')
pickle.dump(mask,out)
out.close()

#To open the mask file from the pickle
with open("map_mask.pickle","r") as mask_file:
    mask_test=pickle.load(mask_file)
    
#Use the mask to only keep the peaks in the Fextr-Fc map that are also observed in the Fo-Fo maps
if np.all(grid_1 != grid_2):
    print("Cannot apply the mask")
elif np.all(first_1 != first_2):
    print("Cannot apply the mask")
elif np.all(unit_cell_1 != unit_cell_1):
    print("Cannot apply the mask")
else:
    data_2_masked = np.multiply(data_2, mask)

        

