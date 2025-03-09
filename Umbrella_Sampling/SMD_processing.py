#!/usr/bin/env python

import pytraj as pt
import os 
from math import *
import argparse 
import numpy as np 

parser = argparse.ArgumentParser()
parser.add_argument("-sp","--spacing",help="US windown spacing in Angstroms",required=True)
parser.add_argument("-p","--topology",help="topology file",required=True)
parser.add_argument("-c","--trajectory",help="trajectory DCD file",required=True)
parser.add_argument("-s","--pdb",help="trajectory PDB file",required=True)
args = parser.parse_args()

def main():
	# Load SMD trajectory file
	print("\t- Loading trajectory and topology files...")
	traj = pt.load(args.trajectory,args.topology)

	# Calculating COM of Fc and FcReceptor
	fc = pt.center_of_mass(traj,mask=":1-116@CA")
	fcr = pt.center_of_mass(traj,mask=":117-166@CA")

	# Splitting generated frames from SMD simulation
	def split_pdb(input_pdb, output_prefix):
		with open(input_pdb, 'r') as infile:
			model_index = 0
			outfile = None
			
			for line in infile:
				if line.startswith('MODEL'):
					if outfile:
						outfile.close()
					model_index += 1
					outfile = open(f'{output_prefix}_{model_index}.pdb', 'w')
				
				if outfile:
					outfile.write(line)
				
				if line.startswith('ENDMDL'):
					outfile.close()
					outfile = None

			if outfile:
				outfile.close()
	print("\t- Splitting PDB trajectory file...")
	split_pdb(args.pdb, 'frame')

	# Calculating COM distance between Fc and FcReceptor
	print("\t- Calculating COM Distance and Extracting US Frames...")
	com_distance = []
	for i in range(len(fc)):
		com_distance.append(sqrt((fc[i][0] - fcr[i][0])**2
		 + (fc[i][1] - fcr[i][1])**2 
		 + (fc[i][2] - fcr[i][2])**2))

	# Selecting windows for US simulation
	start = 0
	frames = []
	while start < len(com_distance):
		frames.append(start+1) # To match the PDB file numbering
		# Find the next frame that meets the spacing criterion
		for i in range(start + 1, len(com_distance)):
			if abs(com_distance[i] - com_distance[start]) >= 0.5:
				start = i 
				break
		else:
			break

	# Adding Box vector coordinates to each frame 
	print("\t- Processing PDB Files...")
	with open(args.pdb,"r") as f:
		lines = f.readlines()
		for i in lines:
			if "CRYST1" in i:
				cryst = i
				break

	files = [x for x in os.listdir() if x.startswith("frame_") and int(x.split(".")[0].split("_")[-1]) in frames]
	
	for i in files:
		with open(i,"r") as f:
			lines = f.readlines()
			lines.insert(0,cryst)
		with open(i,"w") as f:
			for i in lines:
				f.write(i)
				
	# Delete unwanted frames
	def CreateFolders(i:int):
		os.mkdir("f_{}/".format(i))
		os.system("cp frame_{}.pdb f_{}/".format(i,i))
		os.system("cp plumed.dat f_{}".format(i))
		os.chdir("f_{}/".format(i))
		os.system("mv frame_{}.pdb protein.pdb".format(i)) 
		os.system(""" sed 's|AT=|AT={}|g' plumed.dat > plumed_file.dat""".format((com_distance[i-1]) / 10)) # The index value in COM_distance starts from 0 and frames start from 1; converting angstrom to nm (1 nm = 10 angstroms)
		os.chdir("../")
	
	for i in os.listdir():
		if i in files:
			num = int(i.split("_")[1].split(".")[0])
			CreateFolders(num)
		elif i.startswith("frame") and i not in frames:
			os.system("rm {}".format(i))

	os.system("rm *.pdb")
	print("Number of Umbrella windows:",len([x for x in os.listdir("./") if os.path.isdir(x)]))
	
if __name__ == "__main__":
	main()
