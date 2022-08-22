#!/usr/bin/env python

#################################################################################################################
##
## This libary for reading in db2 file
## 
#################################################################################################################
## Writen by Trent Balius in the Shoichet Lab, UCSF in 2015
#################################################################################################################

import argparse
import cmath
import logging
import math
import os
import sys

import mol2
#import mol2_debug as mol2


#################################################################################################################
#################################################################################################################
# data structure to store information about each residue with the docked ligand.
class db2Mol:
    def __init__(self,header,atom_list,bond_list,coord_list,seg_list,conformer_list):
        self.header         = header
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
class db2segment: # 
    def __init__(self,num,start,stop):
        self.num = int(num)
        self.start = int(start)
        self.stop = int(stop)
class db2conformer: # C
    def __init__(self,num,seglist):
        self.num = int(num)
        self.seglist = seglist
     


#################################################################################################################
#################################################################################################################
#def read_Mol2_filehandel(filehandel,startline):
#    lines  =  filehandel.readlines()
#def read_Moldb2_lines(lines,startline):
def read_Moldb2_file(source):
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
#S setnum #lines #confs_total broken hydrogens omega_energy
#S setnum linenum #confs confs [until full column]
#D clusternum setstart setend matchstart matchend #additionalmatching
#D matchnum color x y z
#E 


    # reads in data from multi-Mol2 file.
    with mol2.maybe_open(source, 'r') as lines:
        headers = []
        for line in lines:
             linesplit = line.split() #split on white space
             if(line[0] == "M"): # ATOM 
                 headers.append(line[1:-1].split())
                 atomlist  = []
                 bondlist  = []
                 coordlist = []
                 seglist   = []
                 conflist  = []

             elif(line[0] == "A"): # ATOM 
                atomnum    = linesplit[1]
                atomname   = linesplit[2]
                atomtype   = linesplit[3]
                atomcharge = linesplit[6]
                tempatom = db2atom(atomcharge,atomtype,atomname,atomnum)
                atomlist.append(tempatom)

             elif(line[0] == "B"): # BOND 
                bondnum  = linesplit[1]
                atom1 = linesplit[2]
                atom2 = linesplit[3]
                bondtype = linesplit[4]
                tempbond = db2bond(atom1,atom2,bondnum,bondtype)
                bondlist.append(tempbond)
             elif(line[0] == "X"): # COORDS
                coordnum = linesplit[1]
                atomnum  = linesplit[2]
                segnum   = linesplit[3]
                X        = linesplit[4]
                Y        = linesplit[5]
                Z        = linesplit[6]
                temp_coord = db2coords(coordnum,atomnum,segnum,X,Y,Z)
                coordlist.append(temp_coord)
             #elif(line[0] == "R"): # Rigid
             #   print line
             elif(line[0] == "C"): # Segment 
                confnum    = linesplit[1]
                coordstart = linesplit[2]
                coordstop  = linesplit[3]
                tempseg = db2segment(confnum, coordstart, coordstop)
                seglist.append(tempseg)
                numold = 1
                fristflag = True
             elif(line[0] == "S"): # set -- Conformer 
                num = int(linesplit[1])
                num2 = int(linesplit[2])
                if (fristflag):
                    fristflag = False
                    segnum_list = []
                elif (numold != num): # we know when it is a new conf when this number changes. 
                    tempconf = db2conformer(num,segnum_list)
                    conflist.append(tempconf)
                    segnum_list = []
                    # This fist line does not contain the segment information
                    # The second, and higher lines have more information
                else: # there may be multiple lines for enumerating sagments for one conformer. 
                    numofseg = linesplit[3]
                    for i in range(4,len(linesplit)):
                        segnum_list.append(int(linesplit[i]))
                numold = num
             elif(line[0] == "E"): # ATOM 
                 #if (len(segnum_list) > 0): # this is to put the last conformation in the the list
                 tempconf = db2conformer(num,segnum_list)
                 conflist.append(tempconf)

                 logging.debug("atomnum =", len(atomlist))
                 logging.debug("bondnum =", len(bondlist))
                 logging.debug("coordnum =", len(coordlist))
                 logging.debug("segnum =", len(seglist))
                 logging.debug("confnum =", len(conflist))
                 tempmol = db2Mol(headers, atomlist, bondlist, coordlist, seglist, conflist)  # this is an ensomble of conformation 
                 yield tempmol
             else:
                 logging.debug(line[0] + " is not found in the if statments.")


#################################################################################################################
#################################################################################################################

def convert_db2_to_mol2(db2mols, chimera_headers=False):
    # loop over each molecule
    for mol_idx, mol in enumerate(db2mols, start=1):
         # loop over each conformer in the molcule
         mol2mols = []
         mol_name = str(mol_idx)
         headers = []

         for header_idx, header in enumerate(mol.header):
            num_tokens = len(header)
            if header_idx == 0:
                if num_tokens > 0:
                    mol_name = header[0]
                    headers.append(('Name', mol_name))
                if num_tokens > 1:
                    headers.append(('Protonation', header[1]))
            if header_idx == 2 and num_tokens > 0:
                headers.append(('Smiles', header[0]))
            if header_idx == 3 and num_tokens > 0:
                headers.append(('Long Name', header[0]))

         header = ''
         if chimera_headers:
            header = '\n'.join(['##########{name}: {value}'.format(name=name.rjust(22), value=value.rjust(21))
                                for name, value in headers])
            header += '\n'
             
         for conf in mol.conformer_list:
              # the conformer is defined as a set of segement of the molecule
              mol2atomlist = []
              residue_list = {}
              for segint in conf.seglist:
                  segment =  mol.seg_list[segint-1]
                  logging.debug("\t".join(map(str, [segment.num, segment.start, segment.stop])))
                  # the segement point to a bunch of coordenates, we know what atom the coordenate coresponds to. 
                  for coordint in range(segment.start,segment.stop+1):
                      coord = mol.coord_list[coordint-1]
                      logging.debug("\t".join(map(str, [coord.num, coord.atomnum, coord.segnum,coord.X,coord.Y,coord.Z])))
                      tempatom = mol.atom_list[coord.atomnum-1]
                      #X,Y,Z,Q,type,name,num,resnum,resname):
                      res_num = 1
                      resname = "lig"
                      mol2atom = mol2.atom(coord.X,coord.Y,coord.Z,tempatom.Q,tempatom.type,tempatom.name,tempatom.num,res_num,resname)
                      if residue_list.has_key(res_num):
                         residue_list[res_num].append(mol2atom)
                      else:
                         residue_list[res_num] = [mol2atom]
                      mol2atomlist.append(mol2atom)
              mol2bondlist = []
              for bond in mol.bond_list: 
                  #,a1_num,a2_num,num,type
                  mol2bond = mol2.bond(bond.a1_num,bond.a2_num,bond.num,bond.type)
                  mol2bondlist.append(mol2bond)
              mol2mol = mol2.Mol(header,mol_name,mol2atomlist,mol2bondlist,residue_list)
              mol2mols.append(mol2mol)
         yield mol2mols
         #return allmol2s
#################################################################################################################
#################################################################################################################
def main(args=sys.argv, stdout=sys.stdout):

    parser = argparse.ArgumentParser()
    parser.add_argument('source', nargs='?', type=str, default=None,
                        help='Source db2 file to read [default: -]')
    parser.add_argument('destination', nargs='?', type=str, default=None, 
                        help="Prefix or mol2 file to write [default: - if multi or basename]")
    parser.add_argument('-m', '--multi-mol2', action='store_true', default=False, dest='multi', 
                        help='Write all heirarchies in a single multi mol2 file')
    parser.add_argument('-c', '--chimera-headers', action='store_true', default=False, dest='chimera', 
                        help='Include Chimera mol2 headers (for viewdock)')
    parser.add_argument('--debug', action='store_true', default=False)
    args = args[1:] if args is sys.argv else args

    params = parser.parse_args(args)

    level = logging.DEBUG if params.debug else logging.INFO
    logging.basicConfig(level=level)

    if params.source is None or params.source == '-':
        params.source = sys.stdin
    
    target_tpl = None
    if params.destination is None:
        if params.multi or params.source is sys.stdin or params.source == '-':
            target_tpl = sys.stdout
            logging.info("filename: -")
        else:
            base, ext = os.path.splitext(os.path.basename(params.source))
            if ext == '.gz':
                base, ext = os.path.splitext(base)
            params.destination = base

    if target_tpl is None:
        if params.multi:
            target_tpl = '{}.mol2'.format(params.destination)
        else:
            target_tpl = '{}{{num}}.mol2'.format(params.destination)
        logging.info("filename: {}".format(target_tpl.format(num='#')))

    db2mols = read_Moldb2_file(params.source)
    allmol2s = convert_db2_to_mol2(db2mols, chimera_headers=params.chimera)

    if params.multi:
        target_tpl = mol2.maybe_open(target_tpl, 'a', context=False)

    for heir_idx, mol2mols in enumerate(allmol2s, start=1):
        outfile = target_tpl

        if isinstance(outfile, basestring):
            outfile = target_tpl.format(num=heir_idx)
            logging.info("writing #{} -> {}".format(heir_idx, outfile))
        else:
            logging.info('writing #{}'.format(heir_idx))

        for mol_idx, mol2mol in enumerate(mol2mols, start=1):
            mol2.write_mol2(mol2mol, outfile)

#################################################################################################################
#################################################################################################################

if __name__ == '__main__':
    sys.exit(main(args=sys.argv, stdout=sys.stdout))

