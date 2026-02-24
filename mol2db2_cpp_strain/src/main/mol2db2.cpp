#include "mol2.h"
#include "clash.h"
#include "hydrogens.h"
#include "hierarchy.h"
#include <iostream>
#include <string>
#include <cstring>

using namespace mol2db2;

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " [options]\n";
    std::cout << "\nOptions:\n";
    std::cout << "  -v, --verbose           Verbose output\n";
    std::cout << "  --covalent              Process covalent ligands (remove SiH3 dummy)\n";
    std::cout << "  -a, --time              Show timing information\n";
    std::cout << "  -m, --mol2 FILE         Input mol2 file (default: db.mol2.gz)\n";
    std::cout << "  -n, --namemol2 FILE     Name file (default: name.txt)\n";
    std::cout << "  -s, --solv FILE         Solvation file (default: db.solv)\n";
    std::cout << "  -o, --db FILE           Output db2.gz file (default: db.db2.gz)\n";
    std::cout << "  -t, --atomtype FILE     Atom type conversion file\n";
    std::cout << "  -c, --colortable FILE   Color table file\n";
    std::cout << "  -d, --clash FILE        Clash parameter file\n";
    std::cout << "  -y, --hydrogens FILE    Hydrogen parameter file\n";
    std::cout << "  -r, --noreseth          Don't reset planar hydrogens\n";
    std::cout << "  -z, --norotateh         Don't rotate terminal hydrogens\n";
    std::cout << "  -x, --disttol FLOAT     Distance tolerance (default: 0.001)\n";
    std::cout << "  --limitset NUM          Maximum number of sets (default: 2500000)\n";
    std::cout << "  --limitconf NUM         Maximum number of conformations (default: 200000)\n";
    std::cout << "  --limitcoord NUM        Maximum number of coordinates (default: 1000000)\n";
    std::cout << "  --help                  Show this help message\n";
}

int main(int argc, char* argv[]) {
    /* Default options */
    bool verbose = false;
    bool covalent = false;
    bool timeit = false;
    std::string mol2file = "db.mol2.gz";
    std::string namefile = "name.txt";
    std::string solvfile = "db.solv";
    std::string db2gzfile = "db.db2.gz";
    std::string* atomtypefile = nullptr;
    std::string* colortablefile = nullptr;
    std::string* clashfile = nullptr;
    std::string* hydrogenfile = nullptr;
    bool reseth = true;
    bool rotateh = true;
    double tolerance = 0.001;
    long limitset = 2500000;
    long limitconf = 200000;
    long limitcoord = 1000000;
    
    /* Parse command line arguments */
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "--covalent") {
            covalent = true;
        } else if (arg == "-a" || arg == "--time") {
            timeit = true;
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-m" || arg == "--mol2") {
            if (i + 1 < argc) {
                mol2file = argv[++i];
            }
        } else if (arg == "-n" || arg == "--namemol2") {
            if (i + 1 < argc) {
                namefile = argv[++i];
            }
        } else if (arg == "-s" || arg == "--solv") {
            if (i + 1 < argc) {
                solvfile = argv[++i];
            }
        } else if (arg == "-o" || arg == "--db") {
            if (i + 1 < argc) {
                db2gzfile = argv[++i];
            }
        } else if (arg == "-t" || arg == "--atomtype") {
            if (i + 1 < argc) {
                atomtypefile = new std::string(argv[++i]);
            }
        } else if (arg == "-c" || arg == "--colortable") {
            if (i + 1 < argc) {
                colortablefile = new std::string(argv[++i]);
            }
        } else if (arg == "-d" || arg == "--clash") {
            if (i + 1 < argc) {
                clashfile = new std::string(argv[++i]);
            }
        } else if (arg == "-y" || arg == "--hydrogens") {
            if (i + 1 < argc) {
                hydrogenfile = new std::string(argv[++i]);
            }
        } else if (arg == "-r" || arg == "--noreseth") {
            reseth = false;
        } else if (arg == "-z" || arg == "--norotateh") {
            rotateh = false;
        } else if (arg == "-x" || arg == "--disttol") {
            if (i + 1 < argc) {
                tolerance = std::stod(argv[++i]);
            }
        } else if (arg == "--limitset") {
            if (i + 1 < argc) {
                limitset = std::stol(argv[++i]);
            }
        } else if (arg == "--limitconf") {
            if (i + 1 < argc) {
                limitconf = std::stol(argv[++i]);
            }
        } else if (arg == "--limitcoord") {
            if (i + 1 < argc) {
                limitcoord = std::stol(argv[++i]);
            }
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    if (verbose) {
        std::cout << "mol2db2 - Converting MOL2 to DB2 format\n";
        std::cout << "Input MOL2: " << mol2file << "\n";
        std::cout << "Input SOLV: " << solvfile << "\n";
        std::cout << "Output DB2: " << db2gzfile << "\n";
    }
    
    try {
        /* Read mol2 file */
        if (verbose) {
            std::cout << "Reading MOL2 file..." << std::endl;
        }
        Mol2 mol2data(&mol2file, &namefile, nullptr);
        
        /* Convert dock types and add colors */
        if (verbose) {
            std::cout << "Converting atom types..." << std::endl;
        }
        mol2data.convertDockTypes(atomtypefile);
        mol2data.addColors(colortablefile);
        
        /* Read solvation data */
        if (verbose) {
            std::cout << "Reading solvation data..." << std::endl;
        }
        mol2data.solvData.readSolvFile(solvfile);
        
        /* Handle covalent ligands */
        if (covalent) {
            if (verbose) {
                std::cout << "Processing covalent ligand..." << std::endl;
            }
            auto result = mol2data.removeCovalentDummyAtom();
            std::string covAtomType = std::get<0>(result);
            std::vector<int> indicesList = std::get<1>(result);
            
            mol2data.recolorCovalentAttachment(covAtomType);
            
            /* Remove corresponding solv entries */
            for (int idx : indicesList) {
                if (idx < (int)mol2data.solvData.charge.size()) {
                    mol2data.solvData.charge.erase(
                        mol2data.solvData.charge.begin() + idx);
                    mol2data.solvData.polarSolv.erase(
                        mol2data.solvData.polarSolv.begin() + idx);
                    mol2data.solvData.surface.erase(
                        mol2data.solvData.surface.begin() + idx);
                    mol2data.solvData.apolarSolv.erase(
                        mol2data.solvData.apolarSolv.begin() + idx);
                    mol2data.solvData.solv.erase(
                        mol2data.solvData.solv.begin() + idx);
                }
            }
        }
        
        /* Create clash decider */
        Clash clashDecider(clashfile);
        
        /* Handle hydrogens */
        Hydrogens hydrogenRotater(hydrogenfile);
        
        if (rotateh || reseth) {
            if (verbose) {
                std::cout << "Processing terminal hydrogens..." << std::endl;
            }
            
            hydrogenRotater.findTerminalHydrogens(&mol2data);
            
            if (verbose) {
                std::cout << "Found " << mol2data.hydrogensToRotate 
                         << " rotatable hydrogens" << std::endl;
            }
            
            if (reseth && mol2data.hydrogensToRotate > 0) {
                if (verbose) {
                    std::cout << "Resetting planar hydrogens..." << std::endl;
                }
                hydrogenRotater.resetHydrogens(&mol2data);
            }
            
            if (rotateh && mol2data.hydrogensToRotate > 0) {
                if (verbose) {
                    std::cout << "Rotating hydrogens..." << std::endl;
                }
                //mol2data = hydrogenRotater.rotateHydrogens(mol2data);
                hydrogenRotater.rotateHydrogens_teb(mol2data);
            }
        }
        
        if (verbose) {
            std::cout << "Total conformations: " << mol2data.atomXyz.size() 
                     << std::endl;
        }
        
        /* Build hierarchy */
        if (verbose) {
            std::cout << "Building hierarchy..." << std::endl;
        }
        
        Hierarchy hierarchy(mol2data, clashDecider, tolerance, verbose,
                            timeit, limitset, limitconf, limitcoord);
        
        /* Write output */
        if (verbose) {
            std::cout << "Writing DB2 file..." << std::endl;
        }
        
        hierarchy.write(db2gzfile, verbose, timeit, limitset, "w");
        
        if (verbose) {
            std::cout << "Conversion complete!" << std::endl;
        }
        
        /* Clean up */
        delete atomtypefile;
        delete colortablefile;
        delete clashfile;
        delete hydrogenfile;
        
        return 0;
        
    } catch (const TooBigError& e) {
        std::cerr << "Error: Hierarchy too large!\n";
        std::cerr << "  Conformations: " << e.getConfs() << "\n";
        std::cerr << "  Sets: " << e.getSets() << "\n";
        std::cerr << "  Coordinates: " << e.getCoords() << "\n";
        std::cerr << "Try adjusting --limitset, --limitconf, or --limitcoord\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
