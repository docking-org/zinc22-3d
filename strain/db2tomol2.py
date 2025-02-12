#from RigFragner import db2
#from RigFragner import mol2
from rigfragner_py37 import db2
from rigfragner_py37 import mol2
from rdkit import Chem
import sys, os
import gzip

def main():
    if len(sys.argv) != 2: # if no input
        print("Syntax: python "+sys.argv[0]+" db2file.db2.gz")
        return

    db2_filename  = sys.argv[1]
    #outfileprefix  = sys.argv[2]

    print("db2_filename: " + db2_filename)
    #print "outfileprefix: " + outfileprefix
    db2gz_name = os.path.basename(db2_filename)
    db2gz_path = os.path.dirname(db2_filename)
    protomer_id = db2gz_name.split('.')[0]
    #flag_hac_details = os.path.isfile("hac_details.txt")
    #if flag_hac_details == True:
    #    outputfile = open("hac_details.txt",'a')
    #else:
    #    outputfile = open("hac_details.txt",'w')
    # read in db2.gz file one by one
    file_db2 = gzip.open(db2_filename,'rt')
    EOF = False
    read_okay = False
    count2 = 0
    #exit()
    while True:
        #lines,name,read_okay,EOF = read_one_hierary(file2)
        # read in one hierarchy
        db2mol,lines,name,read_okay,EOF = db2.read_one_hierarchy_from_db2(file_db2)
	#print lines
        #print read_okay
        if (EOF == True):
            print("file finished")
            break
        if (read_okay == False):
            print("%s.%s is broken" % (name,count2))
            continue
        else:
            #pos = file_db2.tell()
	    #print "%s pos: %d" % (name, pos)
            #outputfile.write("%s %d\n" % (name, pos))
            # find the rigid fragment part and return it in a mol2mol format
            #mol2rig = db2.detect_rig_frag_for_one_hierarchy(db2mol)
            mol2mols = db2.convert_db2_to_mol2_from_one_hierarchy(db2mol,name)
            #mol2mol = convert_1st_conformer_to_mol2_from_one_hierarchy(db2mol)
            # write it out
            #outfilename1 = outfileprefix + ".rig.mol2"
            #prefix = name +"."+protomer_id+"."+str(count2)
            prefix = protomer_id+"."+str(count2)
            #outfilename1 = prefix+".rig.mol2"
            outfilename1 = prefix+".mol2"
            flag = os.path.isfile(outfilename1)
            if flag == True:
                os.remove(outfilename1)
            #print outfilename1
            #outfilename2 = outfileprefix + ".conf.mol2"
            #print outfilename2
            for mol2mol in mol2mols:
                flag_exist = os.path.isfile(outfilename1)
                if (flag_exist == False):
                    mol2.write_mol2(mol2mol,outfilename1)
                    #mol2.write_mol2(mol2mol,outfilename2)
                else:
                    #mol2.append_mol2(mol2mol,outfilename2)
                    mol2.append_mol2(mol2mol,outfilename1)
            # convert .rig.mol2 to .rig.smi
            #os.system("/nfs/soft/openbabel/current/bin/obabel -imol2 "+prefix+".rig.mol2 -osmi -O "+prefix+".rig.smi")
            # rdkit read in then calculate hac
            #smi_file = open(prefix+".rig.smi",'r')
            #output_hac = open(outfileprefix+".hac",'w')
            #for line in smi_file:
            #    splitline = line.strip().split('\t')
            #    if (splitline[0] == "smiles"):
            #        continue
            #    #print(splitline[0])
            #    try:
            #        m = Chem.MolFromSmiles(splitline[0])
            #        outputfile.write(splitline[0]+" "+prefix+" "+str(m.GetNumHeavyAtoms())+"\n")
            #    except:
            #        print(prefix)
            #os.remove(prefix+".rig.smi")
            #os.remove(prefix+".rig.mol2")
            count2 = count2 + 1
    file_db2.close()
    return;
#################################################################################################################
#################################################################################################################
main()

