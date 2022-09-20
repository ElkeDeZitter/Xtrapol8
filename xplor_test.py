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