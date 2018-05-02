#!/usr/bin/python2.7 
import time
import sys, os
from os.path import basename
from MDAnalysis import *
import numpy as np
start_time = time.time()
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def homespinner(state):
        if state:
                print "processing...\\",
                bs = '\b'
                syms = ['\\', '|', '/', '-']
                for sym in syms:
                        sys.stdout.write("\b%s" % sym)
                        sys.stdout.flush()
                        time.sleep(.5)
        else:
                sys.stdout.write("\b\n\n")
                sys.stdout.flush()

#routine to convert pdb to pdb with center of mass.
def PDB2COM (pdbname,pdbout):
	import MDAnalysis
	#load pdbfile
	u = MDAnalysis.Universe(pdbname) 
	#select protein molecule
	protein = u.select_atoms('protein')
	#serial COM definition.
	serial=0
	#output variable
	lines=[]
	#iterate over residues.
	for residuo in protein.residues:
		#get resname
		resname=residuo.resname
		#sidechain definition.
		sidechain='protein and not (backbone or name H*)'
		# if GLY keep CA atom coordinates.
		if resname != "GLY":
			selection=protein.select_atoms(('resid ' + str(residuo.resid) + ' and ' + sidechain), updating=True)
		# else select sidechain
		else:
			selection=protein.select_atoms(('protein and resid ' + str(residuo.resid) + ' and name CA'), updating=True)
		CM=selection.center_of_mass()
		#increment COM counter
		serial+=1
		#Make pdb line output.
		line="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM ",serial,"COM","",resname,residuo.segid,residuo.resnum,"",CM[0],CM[1],CM[2],0.00,0.00,"C")
		lines.append(str(line))
	#open output file
        out = open(pdbout, 'w')
	# writefile.
	out.writelines(["%s\n" % item  for item in lines])
	out.close()

# Return a dictionary of residue pair Ai,Aj  and distance cutoff.
def parse_distance_file(inputfile):
	dist_dictionary = {}
	with open(inputfile) as f:
		for line in f:
			#separete files by columns.
			column = line.split()
			#separate first column by residue1 and residue2
			par= column[0].strip(":").split("_")
			par_format =          par[0] + "_" + par[1]
			par_format_inverse =  par[1] + "_" + par[0]
			#add values to dictionary.
			dist_dictionary[par_format] = column[1]
			dist_dictionary[par_format_inverse] = column[1]
	return dist_dictionary



#routine to generan 4 matrices from pdbFile:
	#coordinate matrix  Coordinates(nxn)
	#Distance   matrix  Distances(nxn)
	#cutoff     matrix  cutoff(nxn)
#Return
	#contact_map matrix ContactMap(nxn) contact formed 1 otherwise 0.

def make_matrices(pdbFile):
	import MDAnalysis
	from MDAnalysis.analysis import contacts
	selection_string ='name CA' 
	u = MDAnalysis.Universe(pdbFile)
	ref = MDAnalysis.Universe(pdbFile)
        #select protein molecule
        protein       = u.select_atoms(selection_string)
	#ref_selection = ref.select_atoms(selection_string)
	#coms = contacts.Contacts(u, selection=(selection_string),
        #                    refgroup=(ref_selection), radius=6.0)
	#coms.run()

	

def is_any_closer(r, r0, dist=2.5):
        return np.any(r < dist)

def Contact_finder(sa,sb,selection_atom):
	if selection_atom=='CA':
		print("Entro en CA")
	elif selection_atom=='COM':
		print(' Selected %s as an atom probe to search native contacts' % (selection_atom))
		# Parse distance dictionary from file.
		pair_distances_dicc=parse_distance_file(dist_file)
		# write PDB files with COM of each sidechain.
		#PDB2COM (sa,'stateA.pdb')
		#PDB2COM (sb,'stateB.pdb')
		#print(pair_distances_dicc['A_S'])
		make_matrices('1AKE_solva.pdb')				
	elif selection_atom=='all':
		print("Entro en all")

def dynamic_analizer(topfile, trajfile, topformat, trjformat, dist_file,selection_atom ):
        # cargo la dinamica.
        print('Loading Molecular Dynamic simulations....\n')
        traj=Universe(topfile,trajfile,topology_format=topformat,format=trjformat)
        X = np.empty(shape=[0, 3])
        #creo una lista de residuos
        residues=traj.select_atoms("all").residues.resids
        #print(residues)
        # Defino el numero de frames.
        NumbOfFrames=len(traj.trajectory)
        # Progress bar
	homespinner(True)
        #for ts in traj.trajectory:
	#	a=1
	from MDAnalysis.analysis import contacts
	q1q2 = contacts.q1q2(traj, 'name CA', radius=8)
	q1q2.run()

	f, ax = plt.subplots(1, 2, figsize=plt.figaspect(0.5))
	ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 1], label='q1')
	ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 2], label='q2')
	ax[0].legend(loc='best')
	ax[0].set(xlabel='Frame', ylabel='Fraction Q')
	ax[1].plot(q1q2.timeseries[:, 1], q1q2.timeseries[:, 2], '.-')
	ax[1].set(xlabel='Q1', ylabel='Q2',
           title='2D Native Contacts Analysis.')
	f.savefig('q1q2.pdf')

	#Finish Progress bar.
	homespinner(False)
	print("\n\nDone\n\nThe q1,q2 plots were saved in the current directory...\n")
















#MAIN PROGRAM.
if __name__ == '__main__':
	import sys
	import getopt


	def usage():
		"Print helpful, accurate usage statement to stdout."
		print "Usage: Contact_Map.py [OPTIONS]"
		print
		print "    Description of command..."
		print "         -P         topology_filename    (.pdb/.prmtop format)"
		print "         -X         trajectory_filename  (.nc/.dcd format)"
		print "         -S         selection  (CA,COM,all)"
		print "         -a         pdb file corresponding to state A"
		print "         -b         pdb file corresponding to state B"
		print "    Optional parameters:"
		print "        [-p    topology_format]   default=prmtop"
		print "        [-x    trajectory_format] default=nc (netcdf)"
		print "        [-d    distance file] default=distances.yml"
		print "							Develop by Elias Lopez (2018) contact:eliaslopez at qb.fcen.uba.ar\n"


	# process command arguments
	try:
		opt_list, args = getopt.getopt(sys.argv[1:], 'P:X:p:x:d:s:a:b:')
	except getopt.GetoptError, msg:
		print 'Contact_Map.py: %s' %msg
		usage()
		sys.exit(2)

	# initialize required parameters
	#-l: ligand
	topfile =  None
	trajfile = None
	topformat = 'prmtop'
	trjformat = 'nc'
	dist_file='distances.yml'
	selection_atom = 'CA'
	sa=None
	sb=None

	#'l:vo:d:A:CKU:B:R:MFI:Zgs'
	for o, a in opt_list:
		if o in ('-p', '--pf'):
			topformat = a
		if o in ('-P', '--p'):
			topfile = a
		if o in ('-X', '--x'):
			trajfile = a
		if o in ('-x', '--x'):
			trjformat = a
		if o in ('-d', '--d'):
			dist_file = a
		if o in ('-s', '--s'):
                        selection_atom = a
		if o in ('-a', '--a'):
                        sa = a
		if o in ('-b', '--b'):
                        sb = a
		if o in ('-h', '--'):
			usage()
			sys.exit()
	if not  (sa or sb):
		print 'solventclust.py: pdbs files must be specified.'
		usage()
		sys.exit()
	else:
		# ejecuto el programa propiamente dicho.
		print(topfile, trajfile, topformat, trjformat, dist_file,selection_atom,sa,sb)
		print("\n\nCalculating contacts ...\n\n")
	        Contact_finder(sa,sb,selection_atom)	
		#dynamic_analizer(topfile, trajfile, topformat, trjformat, dist_file, selection_atom)
		print("\n\nIt took--- %10.3f seconds ---" % (time.time() - start_time))

