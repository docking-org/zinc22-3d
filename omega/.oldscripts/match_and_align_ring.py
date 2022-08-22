#!python3 

# this script will find a the matching atoms between ligand and ring.
# it will read in the ligand conformations generated with omega. 
# we need to check if they are already aligned.  most will be aligned.  
# if they are not aligned, then calculate the matrix to move ligand from its position to the orgin and then rotate it to a comon frame
#   
import sys
import copy
import math
import mol2_tebsm as mol2



def mag(v): 
   mag_v = 0
   for e in v: 
       mag_v = mag_v + e**2.0
   mag_v = math.sqrt(mag_v)
   return mag_v

def normal(v):
    mag_v = mag(v)
    return [v[0]/mag_v,v[1]/mag_v,v[2]/mag_v]

def dot_product(v1,v2):
    
   return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])
   

def cross_product(v1,v2):
    c_x =  v1[1]*v2[2] - v1[2]*v2[1]
    c_y = -v1[0]*v2[2] + v1[2]*v2[0]
    c_z =  v1[0]*v2[1] - v1[1]*v2[0]
    #print v1
    #print v2
    #print c_x,c_y,c_z
    return [c_x,c_y,c_z]

def cross_product_normal(v1,v2):
    n = cross_product(v1,v2)
    mag_n = mag(n)
    return [n[0]/mag_n,n[1]/mag_n,n[2]/mag_n]

def calc_angle(v1, v2):
# this function calculates the cos t = v1 dot v2 / mag|v1| times mag|v2|.  cos t is returned. 
   dot = dot_product(v1,v2)
   v1mag = mag(v1)
   v2mag = mag(v2)
   cos_theata = dot/(v1mag*v2mag)
   if cos_theata > 1.0:
       cos_theata = 1.0
   return cos_theata


def rotate(atoms,M):
# use a matrix
    for atom in atoms:
        X = atom.X
        Y = atom.Y
        Z = atom.Z
        atom.X = X *M[0][0] + Y * M[0][1] + Z *M[0][2]
        atom.Y = X *M[1][0] + Y * M[1][1] + Z *M[1][2]
        atom.Z = X *M[2][0] + Y * M[2][1] + Z *M[2][2]
        #atom.X = atom.X *M[0][0] + atom.Y * M[1][0] + atom.Z *M[2][0]
        #atom.Y = atom.X *M[0][1] + atom.Y * M[1][1] + atom.Z *M[2][1]
        #atom.Z = atom.X *M[0][2] + atom.Y * M[1][2] + atom.Z *M[2][2]
#def rotate(v,M):
#    
#    #v[0] = v[0] *M[0][0] + v[1] * M[0][1] + v[2] *M[0][2]       
#    #v[1] = v[0] *M[1][0] + v[1] * M[1][1] + v[2] *M[1][2]       
#    #v[2] = v[0] *M[2][0] + v[1] * M[2][1] + v[2] *M[2][2]       
#    return v

def translate(atoms,x,y,z):

    for atom in atoms:
        atom.X = atom.X + x
        atom.Y = atom.Y + y
        atom.Z = atom.Z + z

    

def get_ring_atoms(ring,mol):
    match_dict = {}
    for i,ratom in enumerate(ring.atom_list):
        if ratom.type == 'H': 
           continue
        for j,matom in enumerate(mol.atom_list):
            if matom.type == 'H': 
               continue
            if (math.fabs(ratom.X-matom.X)<0.00001): 
               if (math.fabs(ratom.Y-matom.Y)<0.00001): 
                  if (math.fabs(ratom.Z-matom.Z)<0.00001): 
                      match_dict[j]=i
                      print(j,'->',i)

    return match_dict


def print_atom(atom):
    print("%s, %f, %f, %f\n"%(atom.name,atom.X, atom.Y, atom.Z))

def same_atoms(atomlist1,atomlist2):
     # if on coordenate is different then they are not aligned.  
     for i in range(len(atomlist1)):
         atom1 = atomlist1[i]
         atom2 = atomlist2[i]
         #print (math.fabs(atom1.X-atom2.X))
         if (math.fabs(atom1.X-atom2.X)>0.00001): 
             return False
         if (math.fabs(atom1.Y-atom2.Y)>0.00001): 
             return False
         if (math.fabs(atom1.Z-atom2.Z)>0.00001): 
                return False
     return True

             

def basis_set_rot_matrix(three_atom_list):

    #step one.# make basis set which contains . 
    # project onto the z,y plan.

    atom1 = three_atom_list[0]
    atom2 = three_atom_list[1]
    atom3 = three_atom_list[2]

    v1= [atom1.X-atom2.X,atom1.Y-atom2.Y,atom1.Z-atom2.Z]
    v2= [atom3.X-atom2.X,atom3.Y-atom2.Y,atom3.Z-atom2.Z]
    n_v1 = normal(v1)
    n_v2 = normal(v2)
    u = cross_product_normal(v1,v2)
    w = cross_product_normal(v1,u)

    mag_v1 = mag(v1)
    n = [v1[0]/mag_v1,v1[1]/mag_v1,v1[2]/mag_v1]

    print ("u = ",u)
    print ("w = ",w)
    print ("n = ",n)

    M = [[w[0],n[0],u[0]], [w[1],n[1],u[1]], [w[2],n[2],u[2]]]   
    inverseM = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range(3):
       for j in range(3):
           inverseM[j][i] = M[i][j]

    print (M)
    print (inverseM)
    return M, inverseM

def atom_in_bond(a,b):
    if b.a1_num-1 == a or b.a2_num-1 == a: 
       return True
    return False


def bonded_atom_in_atoms(b,a,alist):
    if a in alist:
       return True
    elif b.a1_num-1 == a and b.a2_num-1 in alist:
       return True
    elif b.a2_num-1 == a and b.a1_num-1 in alist:
       return True
    return False

def pick_three_atoms(atoms,mol):
    # this is important for not planer rigid groups. 
    atoms_n = []
    #atoms_n.append(atoms[0])
    for atom in atoms:
        print ("atom = ",atom)
        add_atom = True
        if len(atoms_n)>=3: 
           print("atoms:",atoms_n)
           break
        for bond in mol.bond_list:
            if atom_in_bond(atom,bond) and bonded_atom_in_atoms(bond,atom,atoms_n):
               add_atom = False
               print("bond:",bond.a1_num,bond.a2_num)
               #print("atoms:",atoms_n)
               break
        if add_atom:
               atoms_n.append(atom)
               
            #else: 
            #   print(bond.a1_num,bond.a2_num)
    if len(atoms_n)<3:
       atoms_n = atoms
    return atoms_n

def align(moving_mol,fixed_mol,atomid,sign):

    if len(atomid) <= 4: 
        print ("ERROR. atomid == %s"%atomid)
        exit()

    if not (sign == 1.0 or sign == -1.0): 
        print ("ERROR. sign is not -1 or 1")
        exit()

    M_mirror = [[1,0,0],[0,1,0],[0,0,sign]]
        
    atomid_n = pick_three_atoms(atomid,moving_mol)

    #atom1 = moving_mol.atom_list[atomid[0]]
    #atom2 = moving_mol.atom_list[atomid[1]]
    #atom3 = moving_mol.atom_list[atomid[2]]
    #atom3 = moving_mol.atom_list[atomid[-1]]
    atom1 = moving_mol.atom_list[atomid_n[0]]
    atom2 = moving_mol.atom_list[atomid_n[1]]
    atom3 = moving_mol.atom_list[atomid_n[2]]
    three_atoms = [atom1,atom2,atom3]
    for a in three_atoms:
        print_atom(a)

    #atom1_fix = fixed_mol.atom_list[atomid[0]] 
    #atom2_fix = fixed_mol.atom_list[atomid[1]] 
    ##atom3_fix = fixed_mol.atom_list[atomid[2]] 
    #atom3_fix = fixed_mol.atom_list[atomid[-1]] 
    atom1_fix = fixed_mol.atom_list[atomid_n[0]]
    atom2_fix = fixed_mol.atom_list[atomid_n[1]]
    atom3_fix = fixed_mol.atom_list[atomid_n[2]]

    three_atoms_fix = [atom1_fix,atom2_fix,atom3_fix]

    for a in three_atoms_fix:
        print_atom(a)

    if same_atoms(three_atoms,three_atoms_fix): # do nothing. 
       print("atoms are aligned.")
       return
    # else move moving_mol to fixed_mol. 

    #exit()
    trans_X1 = atom2.X
    trans_Y1 = atom2.Y
    trans_Z1 = atom2.Z 

    translate(moving_mol.atom_list,-trans_X1,-trans_Y1,-trans_Z1)

    for a in three_atoms:
        print_atom(a)

    #step one.# make basis set which contains . 
    # project onto the z,y plan.

    M1, inverseM1 = basis_set_rot_matrix(three_atoms)
    M2, inverseM2 = basis_set_rot_matrix(three_atoms_fix)

    #rotate(three_atoms,M)
    print ("start")
    for a in three_atoms:
        print_atom(a)

    print ("rotate")
    rotate(moving_mol.atom_list,inverseM1)
    for a in three_atoms:
        print_atom(a)
    
    rotate(moving_mol.atom_list,M_mirror)

    print ("rotate back")
    rotate(moving_mol.atom_list,M2)
    #translate(three_atoms,atom2.X,atom2.Y,atom2.Z)
    trans_X2 = atom2_fix.X
    trans_Y2 = atom2_fix.Y
    trans_Z2 = atom2_fix.Z 
    translate(moving_mol.atom_list,trans_X2,trans_Y2,trans_Z2)
    print('compare_atoms:')
    print('moved:')
    for a in three_atoms:
        print_atom(a)
    print('fixed:')
    for a in three_atoms_fix:
        print_atom(a)
    print(" ")
    #exit()
    return

def check_atoms_aligned(mol1,mol2,atomid):
    atomlist1 = []
    atomlist2 = []
    for i in atomid: 
        atomlist1.append(mol1.atom_list[i])
        atomlist2.append(mol2.atom_list[i])
    bool_val = same_atoms(atomlist1,atomlist2) # do nothing. 
    return bool_val

def main():
   
    #in_mol2_file       = sys.argv[1]
    #in_ring_file       = sys.argv[2]
    #in_omega_mol2_file = sys.argv[3]
    #out_omega_mol2_file = sys.argv[4]
    #in_mol2_file       = sys.argv[1]
    in_ring_file       = sys.argv[1]
    in_omega_mol2_file = sys.argv[2]
    out_omega_mol2_file = sys.argv[3]
   
    
    mollist = mol2.read_Mol2_file(in_omega_mol2_file) 
    #mol  = mol2.read_Mol2_file(in_mol2_file)[0] 
    mol  = copy.deepcopy(mollist[0]) 
    ring = mol2.read_Mol2_file(in_ring_file)[0] 

    match = get_ring_atoms(ring,mol)

    #atoms = match.keys()
    atoms = []
    for i in match.keys():
        atoms.append(i)
        print(i,match[i])

    frist = True
    for mol_ele in mollist: 
        align(mol_ele,mol,atoms,-1.0) # align mol_ele, to mol with sel 
        if not check_atoms_aligned(mol_ele,mol,atoms):
           print ("Try agian...")
           align(mol_ele,mol,atoms,+1.0) # align mol_ele, to mol with sel 
           if check_atoms_aligned(mol_ele,mol,atoms):
              print ("Error: neither worked.")
 
        if frist:
           mol2.write_mol2(mol_ele,out_omega_mol2_file)
           frist = False
        else:
           mol2.append_mol2(mol_ele,out_omega_mol2_file)

main()

