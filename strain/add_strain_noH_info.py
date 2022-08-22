#!/usr/bin/env python3.7
import os, sys, time

# ben 8/6/20
# torsion strain module is no longer separate from DOCK, like it was before
from Torsion_Strain import calc_strain_mol2

###
## written by Jiankun Lyu (and Benjamin Tingle!)

def main():

    pose_file = sys.argv[1]

    # added by ben
    base_name = pose_file.split('/')[-1]
    base_path = '/'.join(pose_file.split('/')[0:-1])
    mol_out_name = base_name.split(".mol2")[0] + ".strain.mol2"
    mol_out_path = os.path.join(base_path, mol_out_name)

    # ben 8/6/20
    # changed strain calculation to function invocation
    # csv data stored in memory rather than to temp file

    xml_path = os.path.join(os.path.dirname(__file__), "TL_2.1_VERSION_6.xml")
    
    start = time.time()
    csv_data = calc_strain_mol2(pose_file, xml_path)
    print(time.time() - start)

    tE_list = []
    pE_list = []
    for line in csv_data:
        if line[1] == "NA":
            tE = -1.0
            pE = -1.0
        else:
            tE = float(line[1])
            pE = float(line[5])
        tE_list.append(tE)
        pE_list.append(pE)

    output = open(mol_out_path, 'w')
    open_pose = open(pose_file, 'r')
    count = 0

    for line in open_pose:
        if line.find("USER_CHARGES") > -1:
            output.write(line)
            #output.write("TotalStrain = %.2f\n" % tE_list[count])
            output.write("Strain = %.2f %.2f\n" % (tE_list[count],pE_list[count]))
            count = count + 1
        else:
            output.write(line)

    open_pose.close()
    output.close()

main()
