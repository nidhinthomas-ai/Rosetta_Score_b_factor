#!/usr/bin/env python
# coding: utf-8

import os
import sys
# import numpy as np

# import py3Dmol
# import matplotlib as mpl
# import matplotlib.pyplot as plt

import argparse
import collections


def get_parser():
    """"Return a command line parser for this script."""
    parser = argparse.ArgumentParser(
        description="This script extracts the Rosetta Energy scores per residue from "
        "Rosetta generated pdb files and adds them to the PDB file as b-factor.")
    parser.add_argument(
        "-p",
        "--pdb",
        dest="pdb_file",
        required=True,
        help="Input PDB file from which the b-factor is extracted.")
    parser.add_argument(
        "-b",
        "--b_factors",
        dest="b_factors_file",
        required=True,
        help="Output b-factor file with the format <Rosetta residue number> <Rosetta energy score>")
    parser.add_argument(
        "-o",
        "--output_file",
        dest="output_file",
        required=True,
        help="The output PDB file produced by this script.")

    return parser

def read_PDB_Rosetta(pdbfilename, split_string="pose"):
    
    with open (pdbfilename, "r") as pdbfile:
        
        score = False
        
        score_list = []
        
        for line in pdbfile:
            
            if line.split(" ")[0] == split_string:
                
                score = True
                
            if score == True:
                
                score_list += [str(line)]
                
    with open ("{0}.bfactor".format(pdbfilename.split('.')[0]), "w") as bfactor_score:
        
        for line in score_list[1:-2]:
            
                bfactor_score.write((line.split(" ")[0]).split("_")[-1] +" "+line.split(" ")[-1])    
                
    with open ("{0}.score".format(pdbfilename.split('.')[0]), "w") as scorefile:
        
        for line in score_list:
            
            scorefile.write(line)

def write_pdb_Rosetta (pdb_filename, bfactor_filename, pdb_output_filename):
    
    with open (pdb_filename, "r") as pdb_file:
        
        pdb_lines = []
        
        for line in pdb_file:
            
            if line[0:6] == "ATOM  ":
                
                pdb_lines += [str(line)]
                
    bfactor_dict = {}
                
    with open (bfactor_filename, "r") as bfactor_file:
        
        for line in bfactor_file:
                
            bfactor_dict[int(line.split(" ")[0])] = float((line.strip()).split(" ")[1])
    
    with open (pdb_output_filename, "w") as pdb_output_file:
    
        for lineno, line in enumerate(pdb_lines):
            
            if lineno == 0:
            
                first_resnr = int(line[23:26])
            
            resnr = int(line[23:26])
            
            bfactor_iter = resnr - first_resnr + 1
        
            pdb_output_file.write("{prefix}{bfac:6.2f}{suffix}".format(prefix=line[:60], bfac=bfactor_dict[bfactor_iter], suffix=line[66:]))
    
def main(args):

    read_PDB_Rosetta(args.pdb_file, split_string="pose")

    write_pdb_Rosetta (args.pdb_file, args.b_factors_file, args.output_file)

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    main(args)

# If you want to plot the pdb structure along with the color bar, use the script below. 

# fig, ax = plt.subplots(figsize=(6, 1))
# fig.subplots_adjust(bottom=0.5)

# cmap = mpl.cm.jet_r
# norm = mpl.colors.Normalize(vmin=-10, vmax=10)

# cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,norm=norm,orientation='horizontal')
# cb1.set_label('Rosetta Energy',fontsize=20,fontweight='bold')
# cb1.ax.tick_params(labelsize=16) 

# view1 = py3Dmol.view("./rosetta_score_bfactor.pdb", height=600,width=1000)
# view1.setStyle({'cartoon': {'color':'white'}})
# view1.addSurface(py3Dmol.VDW,{'opacity':0.7,'colorscheme':{'prop':'b','gradient':'roygb','min':-10,'max':10}})





