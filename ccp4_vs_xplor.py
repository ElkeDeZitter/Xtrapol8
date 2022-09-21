import iotbx.xplor.map
from iotbx import ccp4_map
import numpy as np
map_name = "/Users/coquelle/IBS_2022/RS_SwissFEL/Fext_grad/FoFo_dark/1ps_7uJ/Fextr_grad_MC_mqFoFo.map"
map_ccp4 = "/Users/coquelle/IBS_2022/RS_SwissFEL/Fext_grad/FoFo_dark/1ps_7uJ/Fextr_grad_MC_mqFoFo.ccp4"
pdb = "/Users/coquelle/IBS_2022/RS_SwissFEL/inputs/newdarkjpnew_001.pdb"
xplor = iotbx.xplor.map.reader(file_name=map_name)
from iotbx.file_reader import any_file

any = any_file(file_name=map_ccp4)
ccp4 = ccp4_map.map_reader(file_name=map_ccp4)




#grid = np.array(xplor.gridding.n, dtype=np.float32)
#first = np.array(xplor.gridding.first, dtype=np.float32)
#unit_cell = np.array(xplor.unit_cell.parameters(), dtype=np.float32)
data_xplor = xplor.data.as_numpy_array()


print(any.file_object)

grid = np.array(ccp4.unit_cell_grid, dtype=np.float32)
print(grid)
first = np.array(ccp4.data.as_double().origin(), dtype=np.float32)
print(first)
unit_cell = np.array(ccp4.unit_cell_parameters, dtype=np.float32)
print(unit_cell)
#data = map_data(ccp4)

n_real = ccp4.unit_cell_grid
print(n_real)
data = ccp4.data.as_numpy_array()
print("ccp4 shape")
print(data.shape)

diff = data - data_xplor
print(diff.min(), diff.max(), diff.mean(), diff.std())

from map_explorer import map_explorer
map_explorer(map_ccp4, pdb, 1.5, 3.0, 4.0, 'ccp4')

map_explorer(map_name, pdb, 1.5, 3.0, 4.0, 'xplor')

#plt.hist(data.flatten(), "auto") #, weights=counts)
#plt.show()


