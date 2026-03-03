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

import pickle
import sys

import numpy as np
from matplotlib import pyplot as plt


class plotalpha:
    """
    Calculation of the alpha and occupancy based on the integrated difference maps.
    """

    def __init__(
        self,
        occupancies,
        extrapolation_results,
        reference,
        outsuffix="qFextr",
        log=sys.stdout,
    ):
        self.occupancies = np.array(occupancies)
        self.alphas = np.reciprocal(occupancies)
        self.pearsonCC, self.pos, self.neg, self.sum = zip(*extrapolation_results)
        self.reference = reference
        self.outsuffix = outsuffix
        self.log = log

    def get_residlist(self):
        """
        Extract the residues from the residue list
        """
        if self.resids_lst != None:
            with open(self.resids_lst) as rs_lst:
                residlst = rs_lst.readlines()
                residlst = [lne for lne in residlst if len(lne) > 0]
            if len(residlst) == 0:
                self.resids_lst = None
            residlst = np.array(residlst)
            residlst_unique = list(set(residlst))
            residlst_unique.remove(
                [lne for lne in residlst_unique if "Resn Resv Chain Alt" in lne][0]
            )
            residlst = np.array(residlst_unique)
            self.residlst = residlst

    def normalize_array(self, array):
        if np.max(array) == 0:
            raise ZeroDivisionError
            # other option: return np.zeros_like(array)
        else:
            return array / np.max(array)

    def alphadetermination(self):
        """
        Alphadetermination with signal enhancement. Uses selected peaks only
        """
        try:
            self.results = self.normalize_array(self.sum)
        except ZeroDivisionError:
            print("no peaks found in at least one of the maps")
            self.results = np.zeros_like(self.sum)
        gooda = np.where(self.results == np.max(self.results))
        goodcc = np.where(self.pearsonCC == np.nanmax(self.pearsonCC))
        if np.max(self.alphas[gooda]) != 0:
            a = np.max(self.alphas[gooda])
            o = np.max(self.occupancies[gooda])
            acc = np.max(self.alphas[goodcc])
            occ = np.max(self.occupancies[goodcc])
            self.alphafound = True
        elif np.max(self.alphas[goodcc]) != 0:
            print(
                "No good alpha value could be found based on map_explorer. The alpha value based on the Pearson correlation between FoFo and mFextr-DF will be used instead. Consider changing map_explorer parameters and run again."
            )
            print(
                "No good alpha value could be found based on map_explorer. The alpha value based on the Pearson correlation between FoFo and mFextr-DFc will be used instead. Consider changing map_explorer parameters and run again.",
                file=self.log,
            )
            a = 1
            o = 1
            acc = np.max(self.alphas[goodcc])
            occ = np.max(self.occupancies[goodcc])

            self.alphafound = True
        else:
            print(
                "No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again.",
                file=self.log,
            )
            print(
                "No good alpha value could be found. Best alpha will be set to 1, occupancy will be set to 1 and subsequently corrected to the highest value in the input list. Consider changing map_explorer parameters and run again."
            )
            a = 1
            o = 1
            acc = 1
            occ = 1
            self.alphafound = False
        return a, o, acc, occ

    def plot_alpha_determination(self, outname):
        # write pickle for GUI
        out = open("alpha_occupancy_determination_%s.pickle" % (self.outsuffix), "wb")
        stats = [
            self.alphas,
            self.occupancies,
            self.pos,
            self.neg,
            self.sum,
            self.reference,
            self.pearsonCC,
            self.alpha,
            self.occ,
            self.alpha_CC,
            self.occ_CC,
        ]

        if self.alphafound:
            stats.append(True)
        else:
            stats.append(False)
        pickle.dump(stats, out)
        out.close()

        fig, axes = plt.subplots(2, 2, figsize=(10, 10))

        if self.reference[0] == 0:
            pos_features = np.zeros(len(self.pos))
        else:
            pos_features = np.asarray(self.pos) / self.reference[0]

        if self.reference[1] == 0:
            neg_features = np.zeros(len(self.neg))
        else:
            neg_features = np.asarray(self.neg) / self.reference[1]

        if self.reference[2] == 0:
            all_features = np.zeros(len(self.sum))
        else:
            all_features = np.asarray(self.sum) / self.reference[2]

        axes[0, 0].plot(
            self.alphas, pos_features, "o", color="green", label="Positive features"
        )
        axes[0, 0].plot(
            self.alphas,
            neg_features,
            "s",
            markersize=5,
            color="red",
            label="Negative features",
        )
        axes[0, 0].plot(self.alphas, all_features, "^", color="k", label="All features")
        axes[0, 0].set_xlim([np.min(self.alphas) * 0.95, np.max(self.alphas) * 1.05])
        axes[0, 0].set_xlabel("Alpha value = 1/occupancy")
        axes[0, 0].set_ylabel("Normalized difference map signal")
        axes[0, 0].legend(
            loc="lower right",
            bbox_to_anchor=(0.94, -0.05, 0.45, 0.5),
            fontsize="x-small",
            framealpha=0.5,
        )

        axes[0, 1].plot(
            self.occupancies,
            pos_features,
            "o",
            color="green",
            label="Positive features",
        )
        axes[0, 1].plot(
            self.occupancies,
            neg_features,
            "s",
            markersize=5,
            color="red",
            label="Negative features",
        )
        axes[0, 1].plot(
            self.occupancies, all_features, "^", color="k", label="All features"
        )
        axes[0, 1].set_xlim(
            [np.min(self.occupancies) * 0.95, np.max(self.occupancies) * 1.05]
        )
        axes[0, 1].set_xlabel("Triggered state occupancy")
        axes[0, 1].set_ylabel("Normalized difference map signal")

        axes[1, 0].plot(
            self.alphas, self.pearsonCC, "X", color="blue", label="PearsonCC"
        )
        axes[1, 0].set_xlim([np.min(self.alphas) * 0.95, np.max(self.alphas) * 1.05])
        axes[1, 0].set_xlabel("Alpha value = 1/occupancy")
        axes[1, 0].set_ylabel("Pearson CC")
        axes[1, 0].legend(
            loc="lower right",
            bbox_to_anchor=(0.84, -0.05, 0.45, 0.5),
            fontsize="x-small",
            framealpha=0.5,
        )

        axes[1, 1].plot(
            self.occupancies, self.pearsonCC, "X", color="blue", label="PearsonCC"
        )
        axes[1, 1].set_ylabel("Pearson CC")
        axes[1, 1].set_xlabel("Triggered state occupancy")
        axes[1, 1].set_xlim(
            [np.min(self.occupancies) * 0.95, np.max(self.occupancies) * 1.05]
        )

        if self.alphafound:
            axes[0, 0].set_title(
                "Alpha determination", fontsize="medium", fontweight="bold"
            )
            axes[0, 1].set_title(
                "Occupancy determination", fontsize="medium", fontweight="bold"
            )

            if self.reference[2] == 0:
                ymax = 1
            else:
                ymax = np.max(self.sum / self.reference[2])
            axes[0, 0].vlines(
                self.alpha, ymin=0, ymax=ymax, color="magenta", linestyles="dashed"
            )
            axes[0, 1].vlines(
                self.occ, ymin=0, ymax=ymax, color="magenta", linestyles="dashed"
            )
            axes[1, 0].vlines(
                self.alpha_CC,
                ymin=0,
                ymax=np.max(self.pearsonCC),
                color="magenta",
                linestyles="dashed",
            )
            axes[1, 1].vlines(
                self.occ_CC,
                ymin=0,
                ymax=np.max(self.pearsonCC),
                color="magenta",
                linestyles="dashed",
            )

        else:
            axes[0, 0].set_title(
                "Alpha determination IMPOSSIBLE", fontsize="medium", fontweight="bold"
            )
            axes[0, 0].text(
                np.min(self.alphas),
                0.5,
                "no peaks found in at least one of the maps\n for the selected residues",
            )
            axes[0, 1].set_title(
                "Occupancy determination IMPOSSIBLE",
                fontsize="medium",
                fontweight="bold",
            )
            axes[0, 1].text(
                np.min(self.occupancies),
                0.5,
                "no peaks found in at least one of the maps\n for the selected residues",
            )
        plt.subplots_adjust(hspace=0.25, wspace=0.5, left=0.09, right=0.88, top=0.95)
        plt.savefig(
            "%s.pdf" % (outname),
            dpi=300,
            transparent=True,
            bbox_inches="tight",
            pad_inches=0,
        )
        plt.savefig("%s.png" % (outname), dpi=300, bbox_inches="tight", pad_inches=0)

        plt.close()

    def estimate_alpha(self):
        print("DIFEFRENCE MAP ANALYSIS")
        print("DIFEFRENCE MAP ANALYSIS", file=self.log)
        print(
            "mFext-DFc maps          | Sum of selected integrated peaks | Pearson CC",
            file=self.log,
        )
        print("Occupancy       alpha         absolute   normalized", file=self.log)
        print("mFext-DFc maps          | Sum of selected integrated peaks | Pearson CC")
        print("Occupancy       alpha         absolute   normalized")

        for i in range(self.alphas.shape[0]):
            if self.reference[2] == 0:
                all_features = 0
            else:
                all_features = self.sum[i] / float(self.reference[2])
            print(
                "{:>5.3f} {:>15.3f} {:>15.3f} {:>10.3f} {:>20.3f}".format(
                    self.occupancies[i],
                    self.alphas[i],
                    self.sum[i],
                    all_features,
                    self.pearsonCC[i],
                ),
                file=self.log,
            )
            print(
                "{:>5.3f} {:>15.3f} {:>15.3f} {:>10.3f} {:>20.3f}".format(
                    self.occupancies[i],
                    self.alphas[i],
                    self.sum[i],
                    all_features,
                    self.pearsonCC[i],
                )
            )

        print("Fobs-Fobs map {:>23.3f}".format(self.reference[2]), file=self.log)
        print("Fobs-Fobs map {:>23.3f}".format(self.reference[2]))

        self.alpha, self.occ, self.alpha_CC, self.occ_CC = self.alphadetermination()
        print(
            "Difference map maximization: Best alpha value is %.3f, meaning an occupancy of %.3f"
            % (self.alpha, self.occ),
            file=self.log,
        )
        print(
            "Difference map maximization: Best alpha value is %.3f, meaning an occupancy of %.3f"
            % (self.alpha, self.occ)
        )

        print(
            "Difference map PearsonCC: Best alpha value is %.3f, meaning an occupancy of %.3f"
            % (self.alpha_CC, self.occ_CC),
            file=self.log,
        )
        print(
            "Difference map PearsonCC: Best alpha value is %.3f, meaning an occupancy of %.3f"
            % (self.alpha_CC, self.occ_CC)
        )

        self.plot_alpha_determination(
            "alpha_occupancy_determination_%s" % (self.outsuffix)
        )

        return self.alpha, self.occ, self.alpha_CC, self.occ_CC
