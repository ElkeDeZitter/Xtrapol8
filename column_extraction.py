"""
authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolac Coquelle - nicolas.coquelle@esrf.fr
Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr

-------
Use truncate from CCP4 for converting intensities to structure factors
"""
from __future__ import print_function
import os,sys
from cctbx import miller, xray
from iotbx.file_reader import any_file


class Column_extraction(object):
    def __init__(self, reflections_ref, reflections_trig, low_res=None, high_res=None, log=sys.stdout):
        self.reflections_ref  = reflections_ref
        self.reflections_trig = reflections_trig
        self.low_res          = low_res
        self.high_res         = high_res
        self.log              = log
        
    def check_if_ano(self, hkl):
        """
        Loop over columns in dataset. The data set is considered anomalous if at least one column with anomalous flag is found.
        """
        #ano = False
        #for array in hkl.file_object.as_miller_arrays():
            #if 1 in ["I(+)" in label for label in array.info().labels]:
                #ano = True
                #break
                
        amplitude_types = [xray.observation_types.amplitude,
                        xray.observation_types.reconstructed_amplitude]
        intensity_types = [xray.observation_types.intensity,]
                 
        ano = False
        for array in hkl.file_object.as_miller_arrays():
            if type(array.observation_type()) in amplitude_types:
                if array.anomalous_flag():
                    ano = True
                    #break
                else:
                    ano = False
                    break
        return ano
    
    def find_column_with_labels(self, hkl, labels):
        """
        Search for column in data set with specific labels
        """
        for array in hkl.file_object.as_miller_arrays():
            if array.info().labels == labels:
                return array
    
    def get_F(self, hkl, labels = None, ano_flag = False):
        """
        Search for columns in the data set with correct anomalous flag.
        If labels is indicated, first a search for those columns take place. This is to try to use the same columns for different data sets.
        if F's are found, these are used, otherwise I's are used and converted to F using truncate
        If the finally chosen column is anomalous, the Friedel mates will be merged as we can't handle anomalous data yet.
        """
        amplitude_types = [xray.observation_types.amplitude,
                        xray.observation_types.reconstructed_amplitude]
        intensity_types = [xray.observation_types.intensity,]
        
        f_obs = i_obs = None
        run_truncate = False
        
        if labels != None:
            array = self.find_column_with_labels(hkl, labels)
            if array != None:
                if type(array.observation_type()) in amplitude_types:
                    run_truncate = False
                    #f_obs = array.map_to_asu()
                    f_obs = array
                    labels = f_obs.info().labels
                    print("Found F's: %s" %(labels), file = self.log)
                    print("Found F's: %s" %(labels))
                if type(array.observation_type()) in intensity_types:
                    run_truncate = True
                    #i_obs = array.map_to_asu()
                    i_obs = array
                    labels = i_obs.info().labels
                    print("Found I's: %s Conversion to Fs with truncate" %(labels), file = self.log)
                    print("Found I's: %s Conversion to Fs with truncate" %(labels))
                    
                if (f_obs == None and i_obs != None): #I's found, need to convert to F's
                    print("Found I's: %s Conversion to Fs with truncate" %(labels), file = self.log)
                    print("Found I's: %s Conversion to Fs with truncate" %(labels))
                    
                    #merge bijvoet mates incase of anomalous data:
                    if ano_flag:
                        i_obs = i_obs.average_bijvoet_mates()
                        
                elif f_obs != None: #F's found
                    if ano_flag:
                        f_obs = f_obs.average_bijvoet_mates()
                    
                else: #No I's and no F's found, problem with input file
                    print("Could not find F's nor I's, please check the input file.", file = self.log)
                    hkl.show_summary(self.log)
                    print("Could not find F's nor I's, please check the input file.")
                    hkl.show_summary()
                
                return f_obs, i_obs, run_truncate, labels

        for array in hkl.file_object.as_miller_arrays():
            #print(array)
            if array.anomalous_flag() == ano_flag:
                #print("ano flag", array.anomalous_flag())
                if type(array.observation_type()) in amplitude_types:
                    run_truncate = False
                    #f_obs = array.map_to_asu()
                    f_obs = array
                    labels = f_obs.info().labels
                    print("Found F's: %s" %(labels), file = self.log)
                    print("Found F's: %s" %(labels))
                    break
                if type(array.observation_type()) in intensity_types:
                    run_truncate = True
                    #i_obs = array.map_to_asu()
                    i_obs = array
                    labels = i_obs.info().labels
                    #print("Found I's: %s" %(labels), file = self.log)
                    #print("Found I's: %s" %(labels))
        
        if (f_obs == None and i_obs != None): #I's found, need to convert to F's
            print("Found I's: %s Conversion to Fs with truncate" %(labels), file = self.log)
            print("Found I's: %s Conversion to Fs with truncate" %(labels))
            
            #merge bijvoet mates incase of anomalous data:
            if ano_flag:
                i_obs = i_obs.average_bijvoet_mates()
                
        elif f_obs != None: #F's found
            if ano_flag:
                f_obs = f_obs.average_bijvoet_mates()
            
        else: #No I's and no F's found, problem with input file
            print("Could not find F's nor I's, please check the input file.", file = self.log)
            hkl.show_summary(self.log)
            print("Could not find F's nor I's, please check the input file.")
            hkl.show_summary()
    
        return f_obs, i_obs, run_truncate, labels

    def run_pointless(self, reflections, reflections_ref, prefix):
        """
        Run pointless with a reference data set. This is to make sure that the reflections have been processed with
        the same indexing possibility. This should normally be done for the triggered data set.
        Afterwards run truncate as usual to convert the intensities to structure factors
        """
        mtz_out, script_pointless = self.write_pointless_input(reflections, reflections_ref, prefix)
        os.system("chmod +x %s" % (script_pointless))
        print("Running pointless to avoid indexing issues", file=self.log)
        print("Running pointless to avoid indexing issues, see %s" % (script_pointless))
        pointless_success = os.system("./%s" % (script_pointless))
        if pointless_success == 0:
            pointless_success = True
        else:
            pointless_success = False
        return mtz_out, pointless_success

        # reflections_pointless = any_file(mtz_out, force_type="hkl", raise_sorry_if_errors=True)
        # return reflections_pointless

    def write_pointless_input(self, reflections, reflections_ref, prefix):
        """
        Prepaare mtz files and script for pointless
        """
        mtz_for_pointless = "%s_forpointless.mtz" % (prefix)
        if type(reflections.observation_type()) == xray.observation_types.intensity:
            #Use the make_mtz_for_truncate function, because it does exactly what we need for pointless
            mtz_for_pointless, inlabels = self.make_mtz_for_truncate(reflections, mtz_for_pointless, "I")
        else:
            mtz_for_pointless, inlabels = self.make_mtz_for_truncate(reflections, mtz_for_pointless, "F")

        mtz_reference_for_pointless = "%s_reference_forpointless.mtz" % (prefix)
        if type(reflections_ref.observation_type()) == xray.observation_types.intensity:
            #Use the make_mtz_for_truncate function, because it does exactly what we need for pointless
            #The reference data set should already be Fs by now
            mtz_reference_for_pointless, inlabels_ref = self.make_mtz_for_truncate(reflections_ref, mtz_reference_for_pointless, "I")
        else:
            mtz_reference_for_pointless, inlabels_ref = self.make_mtz_for_truncate(reflections_ref,
                                                                                   mtz_reference_for_pointless, "F")

        if 'I(+)' in inlabels:
            labin_line = ' '
            for label in inlabels:
                labin_line += 'labin %s ' % (label)
        elif 'I' in inlabels:
            labin_line = 'labin I'
        elif 'F(+)' in inlabels:
            labin_line = ' '
            for label in inlabels:
                labin_line += 'labin %s ' % (label)
        elif 'F' in inlabels:
            labin_line = 'labin F'
        else:
            labin_line = " "

        if 'I(+)' in inlabels_ref:
            labref_line = ' '
            for label in inlabels_ref:
                labref_line += 'labref %s ' % (label)
        elif 'I' in inlabels_ref:
            labref_line = 'labref I'
        elif 'F(+)' in inlabels_ref:
            labref_line = ' '
            for label in inlabels_ref:
                labref_line += 'labref %s ' % (label)
        elif 'F' in inlabels_ref:
            labref_line = 'labref F'
        else:
            labref_line = " "

        mtz_from_pointless = "%s_frompointless.mtz" % (prefix)

        script_out = 'launch_pointless_%s.sh' % (prefix)
        i = open(script_out, 'w')
        i.write("#!/bin/bash\n\
pointless HKLIN %s HKLOUT %s HKLREF %s <<eof_pointless > pointless_%s.log\n\
%s\n\
%s\n\
eof_pointless" %(mtz_for_pointless, mtz_from_pointless, mtz_reference_for_pointless, prefix, labin_line, labref_line))

        i.close()
        return mtz_from_pointless, script_out
    
    def run_truncate(self, reflections, prefix, high_res, low_res):
        """
        Run truncate to convert Is to Fs
        Take anomalous signal into account if present in the two datasets considered
        Afterwards, extract the reflections
        """
        mtz_out, script_truncate = self.write_truncate_input(reflections, prefix, high_res, low_res)
        os.system("chmod +x %s" %(script_truncate))
        print("Running truncate, see %s" %(script_truncate))
        os.system("./%s" %(script_truncate))
        reflections_off = any_file(mtz_out, force_type="hkl", raise_sorry_if_errors=True)
        return reflections_off
        
    def make_mtz_for_truncate(self, reflections, outname, column_root_label="I"):
        """
        Write the reflections to an mtz-file. In principle this should be equal for the reference/triggered data set, although with default column labels.
        """
        reflections_ms = reflections.as_mtz_dataset(column_root_label="%s"%(column_root_label))
        reflections_ms.mtz_object().write(file_name = outname)
        labels = reflections_ms.column_labels()[3:]
        return outname, labels

    def write_truncate_input(self, reflections, prefix, high_res, low_res):
        """
        Prepare everything mtz and script for truncate
        """
        mtz_for_truncate  = "%s_fortruncate.mtz"%(prefix)
        assert type(reflections.observation_type()) == xray.observation_types.intensity
        mtz_for_truncate, inlabels = self.make_mtz_for_truncate(reflections, mtz_for_truncate, "I")
    
        #Normally, the labels labels should be of the form I(+), SIGI(+), I(-), SIGI(-)
        
        if 'I(+)' in inlabels:
            #This is not going to work because truncate needs in addition also IMEAN and SIGIMEAN columns
            labin_line = ' '
            for label in inlabels:
                labin_line+='%s=%s '%(label, label)
            labout_line = " F(+)=F(+) SIGF(+)=SIGF(+) F(-)=F(-) SIGF(-)=SIGF(-)"
        elif 'I' in inlabels:
            labin_line = ' IMEAN=I SIGIMEAN=SIGI'
            labout_line = " F=F SIGF=SIGF"
        else: #with this truncate won't run proporly
            labin_line = " "
            labout_line = " "

        if self.ano:
            ano = "YES"
        else:
            ano = "NO"

        mtz_from_truncate = "%s_fromtruncate.mtz"%(prefix)
        
        script_out = 'launch_truncate_%s.sh'%(prefix)
        i = open(script_out,'w')
        i.write("#!/bin/bash\n\
truncate HKLIN %s HKLOUT %s <<eof_truncate > truncate_%s.log\n\
truncate YES\n\
anomalous %s\n\
resolution %.2f %.2f\n\
plot OFF\n\
header BRIEF BATCH\n\
labin %s\n\
labout %s \n\
eof_truncate"
    %(mtz_for_truncate, mtz_from_truncate, prefix, ano, high_res, low_res, labin_line, labout_line))

        i.close()
        return mtz_from_truncate, script_out
    
    def resolution_cutoff(self, obs):
        """
        Cut data at low and high resolution, only if the data extend beyond the limit.
        """
        dmax, dmin = obs.d_max_min()
        if (self.high_res != None and self.high_res > dmin):
            dmin = self.high_res
        if (self.low_res != None and self.low_res < dmax):
            dmax = self.low_res
        return dmax, dmin

    def extract_columns(self):
        """
        Extract columns:
        1) check if the data should be kept anomalous or not. This is only the case if both data sets are anomalous.
        If only one of them is anomalous, then it will be converted to non-anomalous!
        2) extract Fs if present, else extract Is and convert to Fs using truncate.
        """
        
        if self.check_if_ano(self.reflections_ref)+self.check_if_ano(self.reflections_trig) == 2:
            print("Found anomalous data in both datasets. Xtrapol8 is not yet ready to handle anomalous data.\nData will be converted to non-anomalous structure factors", file=self.log)
            print("Found anomalous data in both datasets. Xtrapol8 is not yet ready to handle anomalous data.\nData will be converted to non-anomalous structure factors")
            self.ano = False #This should become True if we can deal with anomalous data
            ano_ref = ano_trig = True
        elif self.check_if_ano(self.reflections_ref):
            print("One of the input files contains anonalous data but the other not.\nData will be converted to non-anomalous structure factors", file=self.log)
            print("One of the input files contains anonalous data but the other not.\nData will be converted to non-anomalous structure factors")
            self.ano = ano_trig = False
            ano_ref = True
        elif self.check_if_ano(self.reflections_trig):
            print("One of the input files contains anonalous data but the other not.\nData will be converted to non-anomalous structure factors", file=self.log)
            print("One of the input files contains anonalous data but the other not.\nData will be converted to non-anomalous structure factors")
            self.ano = ano_ref = False
            ano_trig = True
        else:
            print("Both data sets contain non-anomalous data", file=self.log)
            print("Both data sets contain non-anomalous data")
            self.ano = ano_ref = ano_trig = False

        print("Reference")
        print("Reference", file=self.log)
        f_obs_ref, i_obs_ref, run_truncate, labels_ref = self.get_F(self.reflections_ref, ano_flag=ano_ref)
        if run_truncate: #Data is I
            dmax, dmin = self.resolution_cutoff(i_obs_ref)
            # Run truncate to convert I to F
            self.reflections_ref = self.run_truncate(i_obs_ref, "reference", dmax, dmin)
            f_obs_ref,_,_,_ = self.get_F(self.reflections_ref, ano_flag=False)
            #f_obs_ref,_,_,_ = self.get_F(self.reflections_ref, ano_flag=ano_ref) #when we can handle anomalous data properly
        else: #Data is F
            dmax, dmin = self.resolution_cutoff(f_obs_ref)
        #This should not be done anymore because bijvoet mates already merged in get_F
        #if (self.check_if_ano(self.reflections_ref) == True and self.ano == False): 
            ##merge Friedel pairs in case reflections_ref is ano but reflections_trig is not.
            #print("Merge Friedel pairs for reference data set",file=self.log)
            #print("Merge Friedel pairs for reference data set")
            #f_obs_ref = f_obs_ref.average_bijvoet_mates()
        f_obs_ref = f_obs_ref.resolution_filter(dmax, dmin)

        print("------")
        print("------", file=self.log)
        print("Triggered")
        print("Triggered", file=self.log)
        f_obs_2, i_obs_2, run_truncate,_ = self.get_F(self.reflections_trig, labels_ref, ano_flag=ano_trig)
        if run_truncate: #Data is I
            dmax, dmin = self.resolution_cutoff(i_obs_2)
            #Run pointless to avoid indexing issues
            mtz_pointless, pointless_success = self.run_pointless(i_obs_2, f_obs_ref, "triggered")
            if (os.path.isfile(mtz_pointless) and pointless_success): #Pointless run correctly
                reflections_pointless = any_file(mtz_pointless, force_type="hkl", raise_sorry_if_errors=True)
                _, i_obs_2, _, _ = self.get_F(reflections_pointless,ano_flag=False)
                #_, i_obs_2, _, _ = self.get_F(reflections_pointless,ano_flag=ano_trig) #when we can handle anomalous data properly
            else: #Pointless run incorrectly
                print("Pointless failed, data not reindexed.")
                print("Pointless failed, data not reindexed.", file=self.log)
            #Run truncate to convert I to F
            self.reflections_trig = self.run_truncate(i_obs_2, "triggered", dmax, dmin)
            f_obs_2,_,_,_= self.get_F(self.reflections_trig, ano_flag=False)
            #f_obs_2,_,_,_= self.get_F(self.reflections_trig, ano_flag=ano_trig) #when we can handle anomalous data properly
        else: #Data is F
            # Run pointless to avoid indexing issues
            dmax, dmin = self.resolution_cutoff(f_obs_2)
            mtz_pointless, pointless_success = self.run_pointless(f_obs_2, f_obs_ref, "triggered")
            if (os.path.isfile(mtz_pointless) and pointless_success): #Pointless run correctly
                reflections_pointless = any_file(mtz_pointless, force_type="hkl", raise_sorry_if_errors=True)
                f_obs_2, _, _, _ = self.get_F(reflections_pointless,ano_flag=False)
                #f_obs_2, _, _, _ = self.get_F(reflections_pointless,ano_flag=ano_trig) #when we can handle anomalous data properly
            else: #Pointless run incorrectly
                print("Pointless failed, data not reindexed.")
                print("Pointless failed, data not reindexed.", file=self.log)
        ##This should not be done anymore because bijvoet mates already merged in get_F
        #if (self.check_if_ano(self.reflections_trig) == True and self.ano == False):
            ##merge Friedel pairs in case reflections_trig is ano but reflections_ref is not
            #print("Merge Friedel pairs for other data sets", file=self.log)
            #print("Merge Friedel pairs for other data sets")
            #f_obs_2 = f_obs_2.average_bijvoet_mates()
        f_obs_2 = f_obs_2.resolution_filter(dmax, dmin)

        print("------")
        print("------", file=self.log)

        #print("f_obs_ref.indices().size()",f_obs_ref.indices().size())
        return f_obs_ref, f_obs_2

class Extrapolated_column_extraction(object):
    def __init__(self, mtz_in, column_labels = "I", log=sys.stdout):
        self.mtz_in        = mtz_in
        self.column_labels = column_labels
        self.log           = log
        
    def write_truncate_input(self):
        """
        Prepare script for truncate
        """    
        
        labin_line = ' IMEAN=%s SIGIMEAN=SIG%s' %(self.column_labels, self.column_labels)
        labout_line = ' F=F SIGF=SIGF'

        mtz_from_truncate = "%s_fromtruncate.mtz"%(self.column_labels.lower())
        
        script_out = 'launch_truncate_%s.sh'%(self.column_labels)
        i = open(script_out,'w')
        i.write("#!/bin/bash\n\
truncate HKLIN %s HKLOUT %s <<eof_truncate > truncate_%s.log \n\
truncate YES\n\
anomalous NO\n\
plot OFF\n\
header BRIEF BATCH\n\
labin %s\n\
labout %s \n\
eof_truncate"
    %(self.mtz_in, mtz_from_truncate, self.column_labels.lower(), labin_line, labout_line))

        i.close()
        return mtz_from_truncate, script_out
    
    def find_column_with_Flabels(self, hkl):
        """
        Search for column in data set with specific labels
        """
        for array in hkl.file_object.as_miller_arrays():
            if array.info().labels == ['F', 'SIGF']:
                return array
    
    def run_truncate(self):
        """
        Run truncate to convert Is to Fs
        Afterwards, extract the reflections
        """
        mtz_out, script_truncate = self.write_truncate_input()
        os.system("chmod +x %s" %(script_truncate))
        print("Running truncate, see %s, to convert Is to Fs" %(script_truncate))
        print("Running truncate, see %s, to convert Is to Fs" %(script_truncate), file=self.log)
        os.system("./%s" %(script_truncate))
        reflections_extrapolate = any_file(mtz_out, force_type="hkl", raise_sorry_if_errors=False)
        return reflections_extrapolate
    
    def reflection_file_converter_arguments(self):
        
        mtz_from_converter = "%s_fromconverter.mtz"%(self.column_labels)
        
        inlabels  = '%s,SIG%s' %(self.column_labels, self.column_labels)
        outlabels = 'F'
        
        argument_line = '%s --massage_intensities --mtz=%s --write_mtz_amplitudes --mtz_root_label="%s" --label="%s"' %(self.mtz_in, mtz_from_converter, outlabels, inlabels)
        
        return mtz_from_converter, argument_line
    
    def run_reflection_file_converter(self):
        """
        Run phenix.reflection_file_converter with --massage_intensities option to convert negative reflections to positive
        Afterwards, extract the reflections
        """
        mtz_out, arguments_converter = self.reflection_file_converter_arguments()
        os.system("phenix.reflection_file_converter %s" %(arguments_converter))
        reflections_extrapolate = any_file(mtz_out, force_type="hkl", raise_sorry_if_errors=False)
        return reflections_extrapolate        
        
        
    def get_Fs_from_truncate(self):
        try:
            reflections_extrapolate = self.run_truncate()
            reflections_F = self.find_column_with_Flabels(reflections_extrapolate)
        except AssertionError:
            reflections_F = None
        return reflections_F
            
    
    def get_Fs_from_reflection_file_converter(self):
        try:
            reflections_extrapolate = self.run_reflection_file_converter()
            reflections_F = self.find_column_with_Flabels(reflections_extrapolate)
        except AssertionError:
            reflections_F = None        
        return reflections_F
