#!/user/bin/python

#################################################################################################################
##
## This libary for reading in db2 file
## 
#################################################################################################################
## Writen by Trent Balius in the Shoichet Lab, UCSF in 2015
#################################################################################################################


import math, sys
import os.path
import cmath
from math import sqrt

import mol2
#import mol2_debu:qg as mol2

#################################################################################################################
#################################################################################################################
# data structure to store information about each residue with the docked ligand.
class db2Mol:
    def __init__(self,header,atom_list,bond_list,coord_list,seg_list,conformer_list):
        self.header         = str(header)
        #self.name           = str(name)
        self.atom_list      = atom_list
        self.bond_list      = bond_list
        self.coord_list     = coord_list
        self.seg_list       = seg_list
        self.conformer_list = conformer_list

class db2atom: # A
    def __init__(self,Q,type,name,num):
        #self.X = float(X)
        #self.Y = float(Y)
        #self.Z = float(Z)
        self.Q = float(Q)
        self.heavy_atom = False
        self.type = type
        self.name = name
        self.num  = int(num)
	#self.resnum  = int(resnum)
	#self.resname = resname
class db2bond: # B
     def __init__(self,a1_num,a2_num,num,type):
        self.a1_num = int(a1_num)
        self.a2_num = int(a2_num)
        self.num = int(num)
        self.type = type
class db2coords: # X
    def __init__(self,num,atomnum,segnum,X,Y,Z):
        self.num = int(num)
        self.atomnum = int(atomnum)
        self.segnum = int(segnum)
        self.X = float(X)
        self.Y = float(Y)
        self.Z = float(Z)
class db2segment: # C 
    def __init__(self,num,start,stop):
        self.num = int(num)
        self.start = int(start)
        self.stop = int(stop)
class db2conformer: # S
    def __init__(self,num,seglist):
        self.num = int(num)
        self.seglist = seglist
     


#################################################################################################################
#################################################################################################################
#def read_Mol2_filehandel(filehandel,startline):
#    lines  =  filehandel.readlines()
#def read_Moldb2_lines(lines,startline):
def read_Moldb2_file(file):
    # reads in data from multi-Mol2 file.

#T ## namexxxx (implicitly assumed to be the standard 7)
#M zincname protname #atoms #bonds #xyz #confs #sets #rigid #Mlines #clusters
#M charge polar_solv apolar_solv total_solv surface_area
#M smiles
#M longname
#[M arbitrary information preserved for writing out]
#A stuff about each atom, 1 per line 
#B stuff about each bond, 1 per line
#X coordnum atomnum confnum x y z 
#R rigidnum color x y z
#C confnum coordstart coordend
#S setnum #lines #confs_total broken hydrogens omega_energy
#S setnum linenum #confs confs [until full column]
#D clusternum setstart setend matchstart matchend #additionalmatching
#D matchnum color x y z
#E 


    # reads in data from multi-Mol2 file.

    file1 = open(file,'r')
    lines  =  file1.readlines()
    file1.close()

    mol_list = []

    header = ''
    for line in lines:
         linesplit = line.split() #split on white space
         if(line[0] == "M"): # ATOM 
             header = header + line[1:-1]
             atomlist  = []
             bondlist  = []
             coordlist = []
             seglist   = []
             conflist  = []

         elif(line[0] == "A"): # ATOM 
            #print line
            #print line[0]
            atomnum    = linesplit[1]
            atomname   = linesplit[2]
            atomtype   = linesplit[3]
            atomcharge = linesplit[6]
            tempatom = db2atom(atomcharge,atomtype,atomname,atomnum)
            atomlist.append(tempatom)
            #print atomnum, atomname, atomtype, "q = ", atomcharge

         elif(line[0] == "B"): # BOND 
            #print line
            #print line[0]
            bondnum  = linesplit[1]
            atom1 = linesplit[2]
            atom2 = linesplit[3]
            bondtype = linesplit[4]
            #print "atom1,atom2,bondnum,bondtype = [" , atom1,atom2,bondnum,bondtype, "]"
            tempbond = db2bond(atom1,atom2,bondnum,bondtype)
            bondlist.append(tempbond)
            #print bondnum, atom1,atom2, bondtype
         elif(line[0] == "X"): # COORDS
            #exit()
            #print line
            #print line[0]
            coordnum = linesplit[1]
            atomnum  = linesplit[2]
            segnum   = linesplit[3]
            X        = linesplit[4]
            Y        = linesplit[5]
            Z        = linesplit[6]
            temp_coord = db2coords(coordnum,atomnum,segnum,X,Y,Z)
            coordlist.append(temp_coord)
            #print coordnum,X,Y,Z
         #elif(line[0] == "R"): # Rigid
         #   print line
         elif(line[0] == "C"): # Segment 
            #print line
            #print line[0]
            confnum    = linesplit[1]
            coordstart = linesplit[2]
            coordstop  = linesplit[3]
            #print confnum, coordstart, coordstop 
            tempseg = db2segment(confnum, coordstart, coordstop)
            seglist.append(tempseg)
            numold = 1
            fristflag = True
         elif(line[0] == "S"): # set -- Conformer 
            #print line
            num = int(linesplit[1])
            num2 = int(linesplit[2])
            #print numold, num
            if (fristflag):
                fristflag = False
                segnum_list = []
            elif (numold != num): # we know when it is a new conf when this number changes. 
                #print "new conformation" 
                tempconf = db2conformer(num,segnum_list)
                conflist.append(tempconf)
                segnum_list = []
                # This fist line does not contain the segment information
                # The second, and higher lines have more information
            else: # there may be multiple lines for enumerating sagments for one conformer. 
                #print "continueing, size of segnum_list = " + str(len(segnum_list))
                numofseg = linesplit[3]
                #print numofseg, len(linesplit) 
                for i in range(4,len(linesplit)):
                    segnum_list.append(int(linesplit[i]))
            numold = num
         elif(line[0] == "E"): # ATOM 
             #if (len(segnum_list) > 0): # this is to put the last conformation in the the list
             tempconf = db2conformer(num,segnum_list)
             conflist.append(tempconf)

             print ("atomnum =", len(atomlist))
             print ("bondnum =", len(bondlist))
             print ("coordnum =", len(coordlist))
             print ("segnum =", len(seglist))
             print ("confnum =", len(conflist))
             tempmol = db2Mol(header, atomlist, bondlist, coordlist, seglist, conflist)  # this is an ensomble of conformation 
             header = ''
             mol_list.append(tempmol)
             #exit()
         else:
             print ("Warrning: " + line[0] + " is not found in the if statments. ")
             #exit()

    return mol_list


#################################################################################################################
#################################################################################################################

def convert_db2_to_mol2(db2mols):

    allmol2s = []
    # loop over each molecule
    for mol in db2mols: 
         # loop over each conformer in the molcule
         mol2mols = []
         for conf in mol.conformer_list:
              #print conf.seglist
              # the conformer is defined as a set of segement of the molecule
              mol2atomlist = []
              residue_list = {}
              for segint in conf.seglist:
                  segment =  mol.seg_list[segint-1]
                  print (segment.num, segment.start, segment.stop)
                  # the segement point to a bunch of coordenates, we know what atom the coordenate coresponds to. 
                  #print segment.start,segment.stop, range(segment.start,segment.stop)
                  for coordint in range(segment.start,segment.stop+1):
                      coord = mol.coord_list[coordint-1]
                      print (coord.num, coord.atomnum, coord.segnum,coord.X,coord.Y,coord.Z)
                      tempatom = mol.atom_list[coord.atomnum-1]
                      #X,Y,Z,Q,type,name,num,resnum,resname):
                      res_num = 1
                      resname = "lig"
                      mol2atom = mol2.atom(coord.X,coord.Y,coord.Z,tempatom.Q,tempatom.type,tempatom.name,tempatom.num,res_num,resname)
                      if residue_list.has_key(res_num):
                         residue_list[res_num].append(mol2atom)
                      else:
                         residue_list[res_num] = [mol2atom]
                      #residue_list[res_num] = [tempatom]
                      mol2atomlist.append(mol2atom)
              mol2bondlist = []
              for bond in mol.bond_list: 
                  #,a1_num,a2_num,num,type
                  mol2bond = mol2.bond(bond.a1_num,bond.a2_num,bond.num,bond.type)
                  mol2bondlist.append(mol2bond)
              mol2mol = mol2.Mol("","",mol2atomlist,mol2bondlist,residue_list)
              mol2mols.append(mol2mol)
         allmol2s.append(mol2mols)
         #exit()
         #return allmol2s
    return allmol2s
#################################################################################################################
#################################################################################################################
def main():
    if len(sys.argv) != 3: # if no input
        print ("Syntax: python db2_to_mol2.py db2file.db2 outprefix ")
        return

    filename  = sys.argv[1];
    outfileprefix  = sys.argv[2];

    print ("filename: " + filename )

    db2mols = read_Moldb2_file(filename)
    allmol2s = convert_db2_to_mol2(db2mols)
    filecount = 0
    for mol2mols in allmol2s:
       outfilename = outfileprefix + "." + str(filecount) + ".mol2"
       print ("writing "+ outfilename )
       count = 0
       for mol2mol in mol2mols:
          if count == 0: 
             mol2.write_mol2(mol2mol,outfilename)
          else:
             mol2.append_mol2(mol2mol,outfilename)
          count = count + 1
       filecount = filecount + 1
    return;
#################################################################################################################
#################################################################################################################
main()

