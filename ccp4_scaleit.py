"""
authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolac Coquelle - nicolas.coquelle@esrf.fr
Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr

-------
Use scaleit from CCP4 for scaling instead of scale_cnslike
"""
from __future__ import print_function
import os, sys
from cctbx import miller
from iotbx.file_reader import any_file

def make_mtz_for_scaleit(f_obs_ref, f_obs_2):
    mtz_out  = "forscaleit.mtz"
    f_obs_ms = f_obs_ref.as_mtz_dataset(column_root_label="F_obs_ref")
    f_obs_ms.add_miller_array(f_obs_2,  column_root_label="F_obs_2")
    f_obs_ms.mtz_object().write(file_name = mtz_out)
    return mtz_out

def write_scaleit_input(mtz_in, b_scaling, high_res):
    mtz_out    = "fromscaleit.mtz"
    script_out = 'launch_scaleit.sh'
    i = open(script_out,'w')
    i.write("#!/bin/bash\n\
scaleit HKLIN %s HKLOUT %s <<eof > scaleit.log\n\
REFINE %s \n\
LABIN FP = F_obs_ref SIGFP = SIGF_obs_ref  FPH1 = F_obs_2 SIGFPH1 = SIGF_obs_2 \n\
eof" %(mtz_in, mtz_out, b_scaling))
    i.close()
    return mtz_out, script_out

def run_scaleit(f_obs_ref, f_obs_2, b_scaling, log=sys.stdout):
    if b_scaling   == 'isotropic':
        b_scaling  = 'ISOTROPIC'
    elif b_scaling == 'anisotropic':
        b_scaling  = 'ANISOTROPIC'
    else:
        b_scaling  = 'SCALE'
        
    mtz_forscaleit = make_mtz_for_scaleit(f_obs_ref, f_obs_2)
    mtz_fromscaleit, script_scaleit = write_scaleit_input(mtz_forscaleit, b_scaling, f_obs_2.d_min())
    os.system("chmod +x %s" %(script_scaleit))
    print("Running scaleit, see %s" %(script_scaleit), file=log)
    os.system("./%s" %(script_scaleit))

    reflections = any_file(mtz_fromscaleit, force_type="hkl", raise_sorry_if_errors=True)
    refl        = reflections.file_object.as_miller_arrays()
    f_obs_2 = [st for i,st in enumerate(refl) if "F_obs_2" in str(refl[i].info())][0]
    
    return f_obs_2

"""
Optional input for scaleit
NOWT\n\
converge NCYC 4 \n\
converge ABS 0.001 \n\
converge TOLR -7 \n\
RESOLUTION %.2f %.2f \n\
"""