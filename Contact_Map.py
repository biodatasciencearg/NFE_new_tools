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

def add_bond(r1,r2):
        return("label add Atoms 0/%i \nlabel add Atoms 0/%i \nlabel add Bonds 0/%i 0/%i \n"% (r1-1,r2-1,r1-1,r2-1))



#list of tuples.
def bond_file_vmd(outfile,ltuples):
	with open(outfile, 'w') as f:
		for tup in ltuples:
			f.write(add_bond(tup[0],tup[1]))
	f.close


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
		line="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format("ATOM ",serial,"CA","",resname,residuo.segid,residuo.resnum,"",CM[0],CM[1],CM[2],0.00,0.00,"C")
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



# procedure to make a square matrix (NxN) from a matrix of dimension N.
def square_matrix_tuple(array_example,typea):
        import numpy as np
        # array length
        alength=len(array_example)
        #
        squareArray = np.empty((alength,alength), dtype=typea)
        for i, x in enumerate(array_example):
                for j, y in enumerate(array_example):
                        squareArray[i][j] = (x, y)
        return squareArray


# routine to write a contact map file from a lists of parameters.
def write_contact_file (initial_contacts_resid, distance_contacts_resid, resname_contacts_resid, outfile, selection_string, cutoff_array):
	# tuple list to save the used tuples.
	tuple_list =[]
	with open(outfile + ".qlist", 'w') as f:
		f.write("resid1 resid2 distance resname1 resname2 sel sel cutoff \n")
		# Iterate over lists.
		for t, d, tr, c in zip(initial_contacts_resid, distance_contacts_resid, resname_contacts_resid, cutoff_array):
			# remove diagonal elements.
			if t[0]!=t[1]:
				# remove duplicate tuples, sorting the resids and defining a unique tuple.
				if t[0]<t[1]:
					tuple_correct = (t[0], t[1])
				else:
					tuple_correct = (t[1], t[0])
				# just keep unique tuples. 
				if not tuple_correct in tuple_list:
					tuple_list.append(tuple_correct)
					# formating output.
					line= " %-10i %-10i %-10.4f %-10s %-10s %-10s %-10s %10.4f\n" % (tuple_correct[0], tuple_correct[1], d, tr[0], tr[1], selection_string, selection_string,c)
					f.write(line)
	#print bond file to load pdb with vmd.
	bond_file_vmd(outfile +"_qlist.vmd",tuple_list)
	f.close
#
def get_custom_matrix(dist_file,resnameArray,size):
	pair_distances_dicc=parse_distance_file(dist_file)
	squareArray = np.empty((size,size), dtype=float)	
	#print(pair_distances_dicc['A_S'])	
	for i, t in enumerate(resnameArray):
		for j, e in enumerate(t):
			key=str(e[0]) + "_" + str(e[1])
			#print(i,j,key,pair_distances_dicc[key])
			squareArray[i][j] = float(pair_distances_dicc[key])
	return squareArray	

# 3-letter code list to 1-letter code.
def convert3to1(aminoList):
	from MDAnalysis.lib.util import convert_aa_code
	out = []
	for r in aminoList:
		out.append(convert_aa_code(r))
	return out	
		

# get native contacts file from a string selection and a pdb file.
def get_contacts_file(pdbIN, selection_string, selection, dist_file):
	import MDAnalysis
	from MDAnalysis.analysis import contacts
	u = MDAnalysis.Universe(pdbIN  + "_COM.pdb")
	#select protein molecule
	protein       = u.select_atoms(selection_string)
	#  resid square array
	residArray   = square_matrix_tuple(protein.resids, [('r1', 'i'),('r2', 'i')])
	# resname square array
	resnameArray = square_matrix_tuple(convert3to1(protein.resnames), [('r1', 'S10'),('r2', 'S10')])
	#
	cutoff=get_custom_matrix(dist_file,resnameArray,len(protein.resnames))
	# contacts routine 
        ca1 = contacts.Contacts(u, selection=(selection_string,selection_string), refgroup=(protein, protein), radius=cutoff)
	# initial contact list 
	initial_contacts_resid  = residArray[ca1.initial_contacts]
	# initial distance list
	distance_contacts_resid = ca1.r0[0][ca1.initial_contacts] 
	# initial  resname list
	resname_contacts_resid  = resnameArray[ca1.initial_contacts]
	cutoff_array = cutoff[ca1.initial_contacts]
	# write the native contacts found.
	write_contact_file(initial_contacts_resid, distance_contacts_resid, resname_contacts_resid, pdbIN, selection, cutoff_array)


def get_contacts(pdbFile, selection, dist_file='distances.yml'):
	if (selection is "COM"):
		selection_string = "protein"
	elif (selection is "CA"):
		selection_string = "name CA"
	get_contacts_file(pdbFile, selection_string, selection, dist_file)
	

def Contact_finder(sa,sb,selection_atom,dist_file):
	if selection_atom=='CA':
		print(' Selected %s as an atom probe to search native contacts' % (selection_atom))
                get_contacts(sa, 'CA')
	elif selection_atom=='COM':
		print(' Selected %s as an atom probe to search native contacts' % (selection_atom))
		# write PDB files with COM of each sidechain.
		pdb_com_name_sa = os.path.splitext(sa)[0] 
		PDB2COM(sa,pdb_com_name_sa + "_COM.pdb")
		#PDB2COM (sb,'stateB.pdb')
		
		get_contacts(pdb_com_name_sa, 'COM')				

def dynamic_analizer(topfile, trajfile, topformat, trjformat, dist_file,selection_atom ):
        # cargo la dinamica.
        print('Loading Molecular Dynamic simulations....\n')
        traj=Universe(topfile,trajfile,topology_format=topformat,format=trjformat)
        #X = np.empty(shape=[0, 3])
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
		print 'Contact_map.py: pdbs files must be specified.'
		usage()
		sys.exit()
	else:
		# ejecuto el programa propiamente dicho.
		print(topfile, trajfile, topformat, trjformat, dist_file,selection_atom,sa,sb)
		print("\n\nCalculating contacts ...\n\n")
	        Contact_finder(sa,sb,selection_atom,dist_file)
		if topfile!=None:
			dynamic_analizer(topfile, trajfile, topformat, trjformat, dist_file, selection_atom)
		print("\n\nIt took--- %10.3f seconds ---" % (time.time() - start_time))

