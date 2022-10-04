import iotbx.xplor.map
from iotbx import ccp4_map
import numpy as np
map_ccp4 = "/Users/coquelle/IBS_2022/RS_SwissFEL/Fext_grad/1ps_2uJ_MC_2/Fextr_grad_MC_mqFoFo.ccp4"
#map_ccp4 = "/Users/coquelle/IBS_2022/RS_SwissFEL/Fext_grad/FoFo_dark/1ps_7uJ/Fextr_grad_MC_mqFoFo.ccp4"
pdb = "/Users/coquelle/IBS_2022/RS_SwissFEL/inputs/newdarkjpnew_001.pdb"
ccp4 = ccp4_map.map_reader(file_name=map_ccp4)







grid = np.array(ccp4.unit_cell_grid, dtype=np.float32)
print(grid)
first = np.array(ccp4.data.as_double().origin(), dtype=np.float32)
print(first)
unit_cell = np.array(ccp4.unit_cell_parameters, dtype=np.float32)
print(unit_cell)

n_real = ccp4.unit_cell_grid
print(n_real)
data = ccp4.data.as_numpy_array()
print("ccp4 shape")
print(data.shape)


from map_explorer import map_explorer
map_explorer(map_ccp4, pdb, 1.5, 3.0, 4.5, 'ccp4')

