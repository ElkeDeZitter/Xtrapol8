# -*- coding: utf-8 -*-
"""
authors and contact information
-------
Elke De Zitter - elke.de-zitter@ibs.fr
Nicolas Coquelle - nicolas.coquelle@esrf.fr
Thomas Barends - Thomas.Barends@mpimf-heidelberg.mpg.de
Jacques Philippe Colletier - jacques-Philippe.colletier@ibs.fr

-------

license information
-------
Copyright (c) 2021 Elke De Zitter, Nicolas Coquelle, Thomas Barends and Jacques-Philippe Colletier
see https://github.com/ElkeDeZitter/Xtrapol8/blob/main/LICENSE

-------
"""

from __future__ import division, print_function
import sys
import numpy as np
import math
import scipy.stats
from matplotlib import pyplot as plt
from iotbx.pdb import hierarchy
from cctbx.array_family import flex
from mmtbx.secondary_structure.find_ss_from_ca import find_secondary_structure

class Map_explorer_analysis(object):
    """
    This class is for the analysis of the output of mapexplorer: sum the integrated peakvolume for each residue and locate residues with highest peaks based on a Z-score analysis (asumes normal distribution when omitting central part of the distribution which is not taken into account + plot peakintegration in function of residue number and secondary structure.
    """
    def __init__(self,
                 peakintegration_file=None,
                 ligands=[],
                 log=sys.stdout):
        self.log = log
        if peakintegration_file != None:
            with open(peakintegration_file) as fle:
                self.peaklines = fle.read().split('\n')
            self.peaklines = [line for line in self.peaklines if len(line.split()) > 1 and ("Resn" not in line)]
        self.check_peaks()
        self.ligands = ligands
        # self.waters = self.check_water(self.peaklines)

    def check_water(self, table):
        """
        search for waters in a table, e.g. self.peaklines. We should not need this anymore
        """
        waters_found = False
        for lne in table:
            if lne.split()[0] == 'HOH':
                waters_found = True
        return waters_found

    def check_peaks(self):
        try:
            if len(self.peaklines) == 0:
                raise IndexError
        except IndexError:
            print("No peaks identified. Choose better map_explorer parameters.", file=self.log)
            print("No peaks identified. Choose better map_explorer parameters.")
            self.peaklines = None
        except AttributeError:
            print("No peakintegration file provided.", file=self.log)
            print("No peakintegration file provided.")
            self.peaklines = None

    def get_sorted_peaklines_ligands(self):
        # peaklines_ligands = [lne for lne in self.peaklines for ligand in self.ligands if ligand in lne]
        peaklines_ligands = [lne for lne in self.peaklines for ligand in self.ligands if (ligand) in lne]
        cnp = [["-".join([line.split()[2], line.split()[0], line.split()[1]]), line.split()[-2]] for line in
               peaklines_ligands]
        cnp = map(lambda x: [x[0], float(x[-1])], cnp)
        cnp = sorted(cnp, key=lambda x: x[0])
        pos_peaklist = [x for x in cnp if x[-1] > 0]
        pos_peaklist = self.peaksum_per_residue(pos_peaklist)
        neg_peaklist = [x for x in cnp if x[-1] < 0]
        neg_peaklist = self.peaksum_per_residue(neg_peaklist)
        return pos_peaklist, neg_peaklist

    def peaklist_corrected(self):
        # peaklist_corrected(self, threshold):
        assert self.peaklines != None
        # peaks = np.array([line.split()[-3] for line in self.peaklines], dtype=np.float32)
        # self.peaks_corr = np.where(peaks<0, peaks+threshold, peaks-threshold)
        peaks = np.array([line.split()[-2] for line in self.peaklines], dtype=np.float32)
        self.peaks_corr = np.where(peaks < 0, peaks - (np.max(np.where(peaks > 0, -np.infty, peaks))), peaks - (np.min(
            np.where(peaks < 0, np.infty,
                     peaks))))  # substract central part which is not taken into account because of threshold

    def residlist_top(self, Z=2):
        """
        Function to get the highest integrated peak areas.
        As only peaks are listed above a certain threshold, we need to subtract the middle part of the distribution to obtain a normal distribution.
        If no peakintion file is provided, then an empty file is returned.
        If the distribution is not Gausian, then all peaks are returned.
        If no peaks have a z-score higher higher than Z, then all peaks are returned.
        """
        outname = "residlist_Zscore%.2f.txt" % (Z)
        residlst = open(outname, "w")
        residlst.write("%4s %4s %4s %3s\n" % ("Resn", "Resv", "Chain", "Alt"))

        # write empty file if no peaks (empty peaklist or no peaklist)
        if self.peaklines == None:
            residlst.close()
            return outname

        # correct distrubution to obtain normal distribution
        self.peaklist_corrected()
        # test for normal distribution
        try:
            _, p = scipy.stats.normaltest(self.peaks_corr)
            if p >= 0.05:
                print(
                    "Peaks not normal distributed (p-value for normality test = %.3f). All peaks will be written to Z-score file. This is however NOT correct!" % (
                        p), file=self.log)
                print(
                    "Peaks not normal distributed (p-value for normality test = %.3f). All peaks will be written to Z-score file. This is however NOT correct!" % (
                        p))
                for line in self.peaklines:
                    residlst.write("%s\n" % (line))
                residlst.close()
                return outname
        except ValueError:
            print(
                "Test for normality could not be performed. Probably not enough peaks. Reconsider the map_explorer parameters",
                file=self.log)
            print(
                "Test for normality could not be performed. Probably not enough peaks. Reconsider the map_explorer parameters")
            for line in self.peaklines:
                residlst.write("%s\n" % (line))
            residlst.close()
            return outname

        # Get highest integrated peak areas
        zscores = scipy.stats.zscore(self.peaks_corr)
        peak_indices = np.where(abs(zscores) > Z)[0]

        if peak_indices.shape[0] > 0:
            for x in peak_indices:
                residlst.write("%s\n" % (self.peaklines[x][:19]))
        else:
            print("No peaks found above Z-score.", file=self.log)
            print("No peaks found above Z-score.")
            for line in self.peaklines:
                residlst.write("%s\n" % (line))
        residlst.close()
        return outname

    def remove_ligands_from_peaklines(self):
        assert self.ligands != None
        assert self.peaklines != None
        self.peaklines_noligands = [line for line in self.peaklines if
                                    len(line.split()) > 1 and line.split()[0] not in self.ligands and line.split()[
                                        0] != 'HOH']

    def get_highest_peak_per_residue(self):
        # to be compatible with secondary_structure plot, waters and ligands have to be removed
        self.remove_ligands_from_peaklines()
        cnp = [[line.split()[2], line.split()[1], line.split()[-2]] for line in self.peaklines_noligands]
        cnp = map(lambda x: [x[0], int(x[1]), float(x[2])], cnp)
        cnp.sort()
        # only keep hihgest peak per residue
        x = 0
        while x < len(cnp) - 1:
            try:
                if ','.join(str(cnp[x][0:2])) == ','.join(str(cnp[x + 1][0:2])):
                    if np.abs(cnp[x][2]) <= np.abs(cnp[x + 1][2]):
                        cnp.remove(cnp[x])
                    else:
                        cnp.remove(cnp[x + 1])
                else:
                    x += 1
            except IndexError:  # due to removed items the list has a smaller length than len(cnp)-2
                break
        return cnp

    def get_sorted_peaklist(self):
        self.remove_ligands_from_peaklines()  # to be compatible with secondary_structure plot, waters and ligands have to be removed
        cnp = [[line.split()[2], line.split()[1], line.split()[-2]] for line in self.peaklines_noligands]
        cnp = map(lambda x: [x[0], int(x[1]), float(x[2])], cnp)
        cnp.sort()
        pos_peaklist = [x for x in cnp if x[2] > 0]
        pos_peaklist = self.peaksum_per_residue(pos_peaklist)
        neg_peaklist = [x for x in cnp if x[2] < 0]
        neg_peaklist = self.peaksum_per_residue(neg_peaklist)
        return pos_peaklist, neg_peaklist

    def peaksum_per_residue(self, peaklist):
        x = 0
        while x < len(peaklist) - 1:
            try:
                if ','.join(str(peaklist[x][0:-1])) == ','.join(str(peaklist[x + 1][0:-1])):
                    peaklist[x][-1] += peaklist[x + 1][-1]
                    peaklist.remove(peaklist[x + 1])
                else:
                    x += 1
            except IndexError:
                print(x)
        return peaklist

    def get_plot_limits(self, *lists):
        """
        find maximum postive or negative value to get total range for plotting, with equal range in positive and negative values.'
        """
        mx = 0
        for lst in lists:
            try:
                m = max([abs(x[-1]) for x in lst])
                if m > mx:
                    mx = m
            except ValueError:
                pass
        return -mx, mx

    def get_ss(self, pdb_file):
        """
        This function does not allow intrachain TER cards as can be added by Refmac. So take care of your TER cards beforehand!!!
        TODO: write function to remove TER cards from PDB (or smarter: search how cctbx-people deal with it). A spossibilty is p.write(pdb_hier.hierarchy.as_pdb_string(crystal_symmetry=pdb_hier.input.crystal_symmetry(),output_break_records=False)) as in functions above but then we still need a way to check if there are intrachain TER cards.
        """
        plt.close()
        pdb_hier = hierarchy.input(file_name=pdb_file)
        hier = pdb_hier.hierarchy
        # hierarchy=iotbx.pdb.input(pdb_file).construct_hierarchy()

        # 1.add the HETATM residues to the list of ligands
        # 2.remove ligands because they might mess up with numbers, independent of definition being HETATM or ATOM
        # and remove atom groups with HETATM
        # 3. remove the empty residue groups and chains
        self.hetatm_residues = list()
        # create a list for the HETATM residues
        for m in hier.models():
            for chain in hier.chains():
                for res_group in chain.residue_groups():
                    for atom_group in res_group.atom_groups():
                        for atom in atom_group.atoms():
                            # for each atom of the pdb file:

                            # 1.add the HETATM residues to the list of ligands
                            if atom.hetero:
                                # in the pdb file, if the line of the atom starts with HETATM
                                self.hetatm_residues.append(atom_group.resname)
                                # add the residue name to the list hetatm_residues
                                atom_group.remove_atom(atom)
                                #removing the HETATM from the atom group
                        self.ligands = self.hetatm_residues + self.ligands
                        # add the list of HETATM to the list of self.ligands
                        self.ligands = list(set(self.ligands))
                        # get rid of the repetition of residue name of the list of ligands

                        # 2.remove ligands and atom groups with HETATM (empty at this point)
                        for lig in self.ligands:
                            if atom_group.resname == lig or atom_group.atoms_size() == 0:
                                # print(chain.id, atom_group.resname)
                                # fonctionne pas rs_group pas terable if atom_group in res_group:
                                try:
                                    res_group.remove_atom_group(atom_group)
                                except:
                                    pass #the atom_group was already removed from the residue group

                    # 3.remove residue group and chain if only composed of HETATM
                    if res_group.atom_groups_size() == 0:
                        chain.remove_residue_group(res_group)
                        #remove the residue group
                if chain.residue_groups_size() == 0:
                    m.remove_chain(chain)
                    #remove the chain

        # for chain in hier.chains():
        # for res_group in chain.residue_groups():
        # for atom_group in res_group.atom_groups():
        # if atom_group.resname=='HOH':

        # for chain in hier.chains():
        # for res_group in chain.residue_groups():
        # for atom_group in res_group.atom_groups():
        # if atom_group.resname=='Na':

        # get the secondary structure
        fss = find_secondary_structure(hierarchy=hier)
        results = fss.get_results()

        # make lists with for each ss the chain ID and start and stop residue
        helix_lst = []
        for helix in results.helices:
            helix_range = helix.start_chain_id, helix.get_start_resseq_as_int(), helix.get_end_resseq_as_int()
            helix_lst.append(list(helix_range))
        helix_lst.sort()
        strand_lst = []
        for sheet in results.sheets:
            for strand in sheet.strands:
                strand_range = strand.start_chain_id, strand.get_start_resseq_as_int(), strand.get_end_resseq_as_int()
                strand_lst.append(list(strand_range))
        strand_lst.sort()

        # Get residues with negative and positive peaks
        if self.peaklines != None:
            pos_peaklist, neg_peaklist = self.get_sorted_peaklist()
        else:
            pos_peaklist = neg_peaklist = [['NotARealList', 0, 0]]

        # Get ligands with positive and negative peaks
        if len(self.ligands) > 0:
            if self.peaklines != None:
                pos_peaklist_ligands, neg_peaklist_ligands = self.get_sorted_peaklines_ligands()
                ligands = list(set([peaklist[x][0] for peaklist in [pos_peaklist_ligands, neg_peaklist_ligands] for x in
                                    range(len(peaklist))]))
                num_lig = len(ligands)
            else:
                pos_peaklist_ligands = neg_peaklist_ligands = [['NotARealList', 0, 0]]
                num_lig = 0
                ligands = []
        else:
            pos_peaklist_ligands = neg_peaklist_ligands = [['NotARealList', 0, 0]]
            ligands = []
            num_lig = 0

        # get the number of chains
        chains = [chain.id for chain in hier.chains()]  # if chain.is_protein()] deleted for DNA
        num_chains = len(set(chains))
        # get the numbers and rows to plot
        if num_chains == 1:
            n_cols = 1
        else:
            n_cols = 2
        n_rows_prot = int(math.ceil(num_chains / n_cols))
        # n_rows_lig  = int(math.ceil(num_lig/n_cols))
        if num_lig == 0:
            n_rows_lig = 0
        else:
            n_rows_lig = 1
        n_rows = n_rows_prot + n_rows_lig

        # initiate plot
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(5 * n_rows, 5 * n_cols),
                                squeeze=False)  # , constrained_layout=True)

        # get minimum and maximum value in order to have the same y-axis on each plot
        mn, mx = self.get_plot_limits(pos_peaklist, neg_peaklist, pos_peaklist_ligands, neg_peaklist_ligands)
        if mn == mx:
            mx += 0.001
            mn -= 0.001

        col = -1
        row = -1
        old_chain_id = ''
        # loop of the protein chains
        for chain in hier.chains():
            # if chain.is_protein():  deleted to be used for DNA too, HETATM already deleted, only protein chain left
            ID = chain.id
            # Check if we have already worked with this chain (can happen because of badly placed TER cards
            if ID != old_chain_id:
                # go the next position in the plot
                if col == 0:
                    col = 1
                else:
                    col = 0
                    row += 1
                # get the helix and strands sublists for the correct chain ID
                # unclear why the next 2 lines had an indentation
                helix_ranges = [range(lne[1], lne[2] + 1, 1) for lne in helix_lst if ID in lne[0]]
                strand_ranges = [range(lne[1], lne[2] + 1, 1) for lne in strand_lst if ID in lne[0]]
                try:
                    # get number of amino acids, won't work if multiple conformers
                    AA = flex.double([chain.residues()[0].resseq_as_int(), chain.residues()[-1].resseq_as_int()])
                except AssertionError:
                    # get number of amino acids when multiple conformers
                    AA = flex.double([chain.conformers()[0].residues()[0].resseq_as_int(),
                                      chain.conformers()[0].residues()[-1].resseq_as_int()])
                # flatten the helix and strand list and convert to flex array
                helix_ranges_flat = [val for sublist in helix_ranges for val in sublist]
                helices = flex.int(helix_ranges_flat)
                strand_ranges_flat = [val for sublist in strand_ranges for val in sublist]
                strands = flex.int(strand_ranges_flat)

                # create flex arrays with the length of the helices and strands, needed for plotting later
                h = flex.int(helices.size(), 0)
                s = flex.int(strands.size(), 0)
                a = flex.int(2, 0)

                # get all the positive peaks and negative peaks and the associated residues
                pos_peaks_chain = [x[2] for x in pos_peaklist if x[0] == ID]
                pos_peaks = flex.double(pos_peaks_chain)
                pos_peak_resids_chain = [x[1] for x in pos_peaklist if x[0] == ID]
                pos_peak_resids = flex.int(pos_peak_resids_chain)
                neg_peaks_chain = [x[2] for x in neg_peaklist if x[0] == ID]
                neg_peaks = flex.double(neg_peaks_chain)
                neg_peak_resids_chain = [x[1] for x in neg_peaklist if x[0] == ID]
                neg_peak_resids = flex.int(neg_peak_resids_chain)

                # Check if everything went well during the generation of the lists
                assert pos_peak_resids.size() == pos_peaks.size()
                assert AA.size() == a.size()

                # The way we need to plot depends on the number and rows, therefore 4 possible cases
                # plot the number of amino acids as a black line
                # plot the helices as magenta circles
                # plot the strands as blue triangles
                # plot the positive peaks with a green bar
                # plot the negative peaks with a red bar
                # if (n_rows > 1 and n_cols > 1):
                axs[row, col].plot(AA, a, color='black', zorder=0, linewidth=0.75)
                axs[row, col].scatter(helices, h, color='magenta', marker='o')
                axs[row, col].scatter(strands, s, color='blue', marker='>')

                axs[row, col].bar(pos_peak_resids, pos_peaks, color='green', width=1.0)
                axs[row, col].bar(neg_peak_resids, neg_peaks, color='red', width=1.0)

                axs[row, col].set_ylim(mn, mx)
                axs[row, col].set_title('Chain %s' % ID, fontsize='medium', fontweight="bold")
                axs[row, col].set_xlabel('Residues')
                axs[row, col].set_ylabel('Summed integrated\npeak volume')

                old_chain_id = ID
            # If we already worked with this chain ID
            else:
                try:
                    # get number of amino acids, won't work if multiple conformers
                    AA = flex.double([chain.residues()[0].resseq_as_int(), chain.residues()[-1].resseq_as_int()])
                except AssertionError:
                    # get number of amino acids when multiple conformers
                    AA = flex.double([chain.conformers()[0].residues()[0].resseq_as_int(),
                                      chain.conformers()[0].residues()[-1].resseq_as_int()])
                # continue plotting the black line
                #if (n_rows > 1 and n_cols > 1):
                    #axs[row, col].plot(AA, a, color='black', zorder=0)
                #elif (n_rows > 1 and n_cols == 1):
                    #axs[row].plot(AA, a, color='black', zorder=0)
                #elif (n_rows == 1 and n_cols > 1):
                    #axs[col].plot(AA, a, color='black', zorder=0)
                #else:
                    #axs.plot(AA, a, color='black', zorder=0)
                axs[row, col].plot(AA, a, color='black', zorder=0)

        # do the same thing for the ligands and water molecules
        if num_lig > 0:
            ligands.sort()
            row += 1
            col = 0

            # get the positive and negative peaks in properly in lists
            n = np.asarray(neg_peaklist_ligands)
            nt = n.transpose()
            p = np.asarray(pos_peaklist_ligands)
            pt = p.transpose()

            # The way we need to plot depends on the number and rows, therefore 4 possible cases
            # plot a black line (for reference and consistency)
            # plot the positive peaks with a green bar
            # plot the negative peaks with a red bar
            # if n_cols > 1: #n_rows will always be larger > 1 (except if something went wrong and no protein present ==> should be checked for nucleic acids!!
            axs[row, col].plot(ligands, flex.int(num_lig, 0), color='black', zorder=0, linewidth=0.75)
            if pt.shape[0] > 0:
                axs[row, col].bar(pt[0], pt[1].astype(np.float), color='green', width=0.25)
            if nt.shape[0] > 0:
                axs[row, col].bar(nt[0], nt[1].astype(np.float), color='red', width=0.25)
            axs[row, col].set_ylim(mn, mx)
            axs[row, col].set_ylabel('Summed integrated\npeak volume')
            axs[row, col].set_title("Ligands", fontsize='medium', fontweight="bold")
            plt.setp(axs[row, col].xaxis.get_majorticklabels(), rotation=90, fontsize='small')

        fig.tight_layout()
        plt.savefig('summed_difference_peaks.pdf', dpi=300, transparent=True)
        plt.savefig('summed_difference_peaks.png', dpi=300)
        plt.close()
