#!/usr/bin/env python3.7
import os, sys

###
## written by Jiankun Lyu

def main():

    pwd = os.getcwd()+"/"
    if len(sys.argv) != 2:
        print("python add_strain_info.py multi_pose.mol2")
        sys.exit()

    count = 0
    pose_file = sys.argv[1]
    #strain_filename = sys.argv[2]
    mol_out_name = pose_file.split(".mol2")[0]
    # run strain on the multi-mol2 file
    # http://wiki.docking.org/index.php/Strain_Filtering
    # change dir to where you install strain filter packages
    strainfilter_path = "/mnt/nfs/ex5/work/jklyu/strainfilter_noH_s"
    pwd = os.getcwd()
    os.chdir(strainfilter_path)
    os.system("python Torsion_Strain.py "+pwd+"/"+pose_file)
    os.chdir(pwd)

    straincal_file = open(mol_out_name+"_Torsion_Strain.csv",'r')
    #total_E_list = []
    #max_E_per_list = []
    #flag_emtpy_file = False
    tE_list = []
    pE_list = []
    for line in straincal_file:
        splitline = line.split(',')
        if len(splitline) <= 2:
            tE = -1.0
            pE = -1.0
        else:
            tE = float(splitline[1])
            pE = float(splitline[5])
        tE_list.append(tE)
        pE_list.append(pE)

    #strain_file = open(strain_filename,'r')
    #for line in strain_file:
    #    splitline = line.strip().split()
    #    tE = float(splitline[1])
    #    tE_list.append(tE)
    #    pE = float(splitline[2])
    #    pE_list.append(pE)
    straincal_file.close()

    output = open(mol_out_name+".strain.mol2", 'w')
    open_pose = open(pose_file, 'r')
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
