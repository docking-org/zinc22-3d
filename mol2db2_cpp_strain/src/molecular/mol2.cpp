#include "mol2.h"
#include "geometry.h"
#include "floydwarshall.h"
#include "shortestpaths.h"
#include "divisive_clustering.h"
#include "unionfind.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <zlib.h>

namespace mol2db2 {

// ============================================================================
// AtomConverter Implementation
// ============================================================================

std::unordered_map<std::string, int> AtomConverter::getDefaultConvertTypes() {
    std::unordered_map<std::string, int> types;
    types["C.3"] = 5;
    types["C.2"] = 1;
    types["C.ar"] = 1;
    types["C.1"] = 1;
    types["N.3"] = 10;
    types["N.2"] = 8;
    types["N.1"] = 8;
    types["O.3"] = 12;
    types["O.2"] = 11;
    types["S.3"] = 14;
    types["N.ar"] = 8;
    types["P.3"] = 13;
    types["H"] = 6;
    types["H-C"] = 7;
    types["Br"] = 17;
    types["Cl"] = 16;
    types["F"] = 15;
    types["I"] = 18;
    types["S.2"] = 14;
    types["N.pl3"] = 8;
    types["LP"] = 25;
    types["Na"] = 19;
    types["K"] = 19;
    types["Ca"] = 21;
    types["Li"] = 20;
    types["Al"] = 20;
    types["Du"] = 25;
    types["Du.C"] = 25;
    types["Si"] = 24;
    types["N.am"] = 8;
    types["S.o"] = 14;
    types["S.O"] = 14;
    types["S.o2"] = 14;
    types["S.O2"] = 14;
    types["N.4"] = 9;
    types["O.co2"] = 11;
    types["C.cat"] = 1;
    types["H.spc"] = 6;
    types["O.spc"] = 11;
    types["H.t3p"] = 6;
    types["O.t3p"] = 11;
    types["ANY"] = 25;
    types["HEV"] = 25;
    types["HET"] = 25;
    types["HAL"] = 25;
    types["Mg"] = 20;
    types["Cr.oh"] = 25;
    types["Cr.th"] = 25;
    types["Se"] = 25;
    types["Fe"] = 25;
    types["Cu"] = 25;
    types["Zn"] = 26;
    types["Sn"] = 25;
    types["Mo"] = 25;
    types["Mn"] = 25;
    types["Co.oh"] = 25;
    return types;
}

AtomConverter::AtomConverter(const std::string* parameterFileName) {
    if (parameterFileName != nullptr) {
        std::ifstream parameterFile(*parameterFileName);
        if (parameterFile.is_open()) {
            std::string line;
            while (std::getline(parameterFile, line)) {
                std::istringstream iss(line);
                std::string key;
                int value;
                if (iss >> key >> value) {
                    convertTypes[key] = value;
                }
            }
            parameterFile.close();
        } else {
            convertTypes = getDefaultConvertTypes();
        }
    } else {
        convertTypes = getDefaultConvertTypes();
    }
    
    /* Build special keys map */
    for (auto it = convertTypes.begin(); it != convertTypes.end(); ++it) {
        const std::string& key = it->first;
        size_t pos = key.find('-');
        if (pos != std::string::npos) {
            std::string firstPart = key.substr(0, pos);
            specialKeys[key] = firstPart;
        }
    }
}

void AtomConverter::printParameters() const {
    std::vector<std::pair<std::string, int>> items(convertTypes.begin(), 
                                                     convertTypes.end());
    std::sort(items.begin(), items.end());
    
    for (size_t i = 0; i < items.size(); i++) {
        std::cout << items[i].first << " " << items[i].second << std::endl;
    }
}

int AtomConverter::convertMol2atomNum(const Mol2& mol2data, int atomNum) const {
    int actualNum = atomNum - 1;
    const std::string& actualName = mol2data.atomType[actualNum];
    
    /* Check if this type needs special bond checking */
    bool needsSpecialCheck = false;
    for (auto it = specialKeys.begin(); it != specialKeys.end(); ++it) {
        if (it->second == actualName) {
            needsSpecialCheck = true;
            break;
        }
    }
    
    if (!needsSpecialCheck) {
        auto it = convertTypes.find(actualName);
        if (it != convertTypes.end()) {
            return it->second;
        }
        return 25; /* Default unknown type */
    } else {
        /* Check each special key */
        for (auto it = specialKeys.begin(); it != specialKeys.end(); ++it) {
            const std::string& key = it->first;
            size_t pos = key.find('-');
            if (pos != std::string::npos) {
                std::string bondedAtom = key.substr(pos + 1);
                if (mol2data.bondedTo(atomNum, bondedAtom, 1, nullptr, false, nullptr)) {
                    return convertTypes.at(key);
                }
            }
        }
        /* No special match found */
        auto it = convertTypes.find(actualName);
        if (it != convertTypes.end()) {
            return it->second;
        }
        return 25;
    }
}

// ============================================================================
// ColorTable Implementation
// ============================================================================

std::string ColorTable::getDefaultColorDefault() {
    return "neutral";
}

std::unordered_map<std::string, int> ColorTable::getDefaultColorInts() {
    std::unordered_map<std::string, int> colors;
    colors["positive"] = 1;
    colors["negative"] = 2;
    colors["acceptor"] = 3;
    colors["donor"] = 4;
    colors["ester_o"] = 5;
    colors["amide_o"] = 6;
    colors["neutral"] = 7;
    return colors;
}

std::vector<std::variant<ColorTable::Rule2, ColorTable::Rule4>> 
ColorTable::getDefaultRulesTable() {
    std::vector<std::variant<Rule2, Rule4>> rules;
    
    rules.push_back(Rule2{"N.4", "positive"});
    rules.push_back(Rule2{"O.co2", "negative"});
    rules.push_back(Rule2{"O.2", "acceptor"});
    rules.push_back(Rule2{"O.3", "acceptor"});
    rules.push_back(Rule2{"S.2", "acceptor"});
    rules.push_back(Rule2{"N.ar", "acceptor"});
    rules.push_back(Rule4{"P.3", 1, "O.co2", "negative"});
    rules.push_back(Rule4{"S.o2", 1, "O.co2", "negative"});
    rules.push_back(Rule4{"N.2", 1, "H", "donor"});
    rules.push_back(Rule4{"N.am", 1, "H", "donor"});
    rules.push_back(Rule4{"N.pl3", 1, "H", "donor"});
    rules.push_back(Rule4{"O.3", 1, "H", "donor"});
    rules.push_back(Rule4{"N.ar", -1, "H", "acceptor"});
    rules.push_back(Rule4{"N.ar", -1, "C.3", "acceptor"});
    rules.push_back(Rule4{"N.ar", 1, "H", "donor"});
    rules.push_back(Rule4{"O.3", 1, "H", "donor"});
    rules.push_back(Rule4{"O.2", 2, "O.3", "ester_o"});
    rules.push_back(Rule4{"O.2", 2, "N.pl3", "amide_o"});
    rules.push_back(Rule4{"O.2", 2, "N.am", "amide_o"});
    rules.push_back(Rule4{"O.2", 2, "N.3", "amide_o"});
    
    return rules;
}

ColorTable::ColorTable(const std::string* parameterFileName) {
    if (parameterFileName != nullptr) {
        std::ifstream parameterFile(*parameterFileName);
        if (parameterFile.is_open()) {
            int phase = 0;
            std::string line;
            
            while (std::getline(parameterFile, line)) {
                std::istringstream iss(line);
                std::vector<std::string> tokens;
                std::string token;
                while (iss >> token) {
                    tokens.push_back(token);
                }
                
                if (tokens.empty()) continue;
                
                if (phase == 0) {
                    defaultColor = tokens[0];
                    phase = 1;
                } else if (phase == 1) {
                    if (tokens.size() == 1 && tokens[0] == "rules") {
                        phase = 2;
                    } else if (tokens.size() >= 2) {
                        colorInts[tokens[0]] = std::stoi(tokens[1]);
                    }
                } else if (phase == 2) {
                    if (tokens.size() == 2) {
                        rulesTable.push_back(Rule2{tokens[0], tokens[1]});
                    } else if (tokens.size() == 4) {
                        rulesTable.push_back(Rule4{tokens[0], std::stoi(tokens[1]), 
                                                   tokens[2], tokens[3]});
                    }
                }
            }
            
            parameterFile.close();
        } else {
            defaultColor = getDefaultColorDefault();
            colorInts = getDefaultColorInts();
            rulesTable = getDefaultRulesTable();
        }
    } else {
        defaultColor = getDefaultColorDefault();
        colorInts = getDefaultColorInts();
        rulesTable = getDefaultRulesTable();
    }
}

void ColorTable::printParameters() const {
    std::cout << defaultColor << std::endl;
    
    std::vector<std::pair<std::string, int>> items(colorInts.begin(), 
                                                     colorInts.end());
    std::sort(items.begin(), items.end(), 
              [](const auto& a, const auto& b) { return a.second < b.second; });
    
    for (size_t i = 0; i < items.size(); i++) {
        std::cout << items[i].first << " " << items[i].second << std::endl;
    }
    
    std::cout << "rules" << std::endl;
    for (size_t i = 0; i < rulesTable.size(); i++) {
        if (std::holds_alternative<Rule2>(rulesTable[i])) {
            const Rule2& rule = std::get<Rule2>(rulesTable[i]);
            std::cout << std::get<0>(rule) << " " << std::get<1>(rule) << std::endl;
        } else {
            const Rule4& rule = std::get<Rule4>(rulesTable[i]);
            std::cout << std::get<0>(rule) << " " << std::get<1>(rule) << " "
                     << std::get<2>(rule) << " " << std::get<3>(rule) << std::endl;
        }
    }
}

int ColorTable::convertMol2color(const Mol2& mol2data, int atomNum) const {
    int actualNum = atomNum - 1;
    const std::string& actualName = mol2data.atomType[actualNum];
    std::string lastColorFound = defaultColor;
    
    /* Apply rules in order */
    for (size_t i = 0; i < rulesTable.size(); i++) {
        if (std::holds_alternative<Rule2>(rulesTable[i])) {
            const Rule2& rule = std::get<Rule2>(rulesTable[i]);
            if (actualName.find(std::get<0>(rule)) == 0) {
                lastColorFound = std::get<1>(rule);
            }
        } else {
            const Rule4& rule = std::get<Rule4>(rulesTable[i]);
            std::string atomType = std::get<0>(rule);
            int bondsAway = std::get<1>(rule);
            std::string otherType = std::get<2>(rule);
            std::string color = std::get<3>(rule);
            
            if (actualName.find(atomType) == 0) {
                if (bondsAway == -1) {
                    /* Not bonded to */
                    if (!mol2data.bondedTo(atomNum, otherType, 1, nullptr, 
                                          false, nullptr)) {
                        lastColorFound = color;
                    }
                } else {
                    /* Bonded to at specific distance */
                    if (mol2data.bondedTo(atomNum, otherType, bondsAway, 
                                         nullptr, false, nullptr)) {
                        lastColorFound = color;
                    }
                }
            }
        }
    }
    
    auto it = colorInts.find(lastColorFound);
    if (it != colorInts.end()) {
        return it->second;
    }
    return 7; /* Default to neutral */
}

// Continued in Part 2...
// mol2.cpp Part 2 - SolvData and Mol2 basic operations

// ============================================================================
// SolvData Implementation
// ============================================================================

SolvData::SolvData() 
    : totalAtoms(0), totalCharge(0.0), totalPolarSolv(0.0),
      totalSurface(0.0), totalApolarSolv(0.0), totalSolv(0.0) {
}

void SolvData::readSolvFile(const std::string& filename) {
    std::ifstream solvfile(filename);
    if (!solvfile.is_open()) {
        std::cerr << "Warning: Could not open solv file: " << filename << std::endl;
        return;
    }
    
    std::string line;
    bool first = true;
    
    while (std::getline(solvfile, line)) {
        std::istringstream iss(line);
        
        if (first) {
            /* First line has summary data */
            if (!(iss >> name >> totalAtoms >> totalCharge >> totalPolarSolv 
                  >> totalSurface >> totalApolarSolv >> totalSolv)) {
                std::cerr << "Error reading solv file first line" << std::endl;
                break;
            }
            first = false;
        } else {
            /* Per-atom data */
            double chg, pol, surf, apol, solv_val;
            if (iss >> chg >> pol >> surf >> apol >> solv_val) {
                charge.push_back(chg);
                polarSolv.push_back(pol);
                surface.push_back(surf);
                apolarSolv.push_back(apol);
                solv.push_back(solv_val);
            } else {
                /* Might be multiple molecules in file - stop here */
                break;
            }
        }
    }
    
    solvfile.close();
}

// ============================================================================
// Mol2 Implementation
// ============================================================================

Mol2::Mol2() {
    blankNew();
}

// TEB added. 
Mol2::~Mol2(){

    atomNum.clear();
    atomName.clear();
    atomXyz.clear();
    inputEnergy.clear();
    inputTotalStrain.clear();
    inputMaxStrain.clear();
    inputHydrogens.clear();
    atomType.clear();
    atomCharge.clear();
    atomBonds.clear();
    
    bondNum.clear();
    bondStart.clear();
    bondEnd.clear();
    bondType.clear();
    
    
    dockNum.clear();
    colorNum.clear();
    
    hydrogenRotAngles.clear();
    dihedrals.clear();
    rotAngles.clear();
    
    
    bondDists.clear();
    bondDistsOrderKeys.clear();
    
    colorConverter.reset();
}
// TEB added. 
// this is the copy constructor
//Mol2::Mol2(const Mol2& mol){
Mol2::Mol2(const Mol2& oldM){
   name       =  oldM.name;
   protName   =  oldM.protName;
   smiles     =  oldM.smiles;
   longname   =  oldM.longname;
                      
   atomNum    =  oldM.atomNum;
   atomName   =  oldM.atomName;
   atomType   =  oldM.atomType;
   atomCharge =  oldM.atomCharge;
   atomBonds  =  oldM.atomBonds;
                      
   bondNum    =  oldM.bondNum;
   bondStart  =  oldM.bondStart;
   bondEnd    =  oldM.bondEnd;
   bondType   =  oldM.bondType;

    /* These are copied for manipulation */
    atomXyz        = oldM.atomXyz;
    inputEnergy      = oldM.inputEnergy;
    inputTotalStrain = oldM.inputTotalStrain;
    inputMaxStrain   = oldM.inputMaxStrain;
    inputHydrogens   = oldM.inputHydrogens;
    xyzCount       = oldM.xyzCount;
    origXyzCount   = oldM.origXyzCount;
                          
    dockNum        = oldM.dockNum;
    colorNum       = oldM.colorNum;
    if (colorConverter) {
        colorConverter.reset(new ColorTable(*(oldM.colorConverter)));
    }

    solvData          = oldM.solvData;
                             
    hydrogenRotAngles = oldM.hydrogenRotAngles;
    dihedrals         = oldM.dihedrals;
    rotAngles         = oldM.rotAngles;
    hydrogensToRotate = oldM.hydrogensToRotate;

}

void Mol2::blankNew() {
    name = "fake";
    protName = "fake";
    smiles = "fake";
    longname = "fake";

    atomNum.clear();
    atomName.clear();
    atomXyz.clear();
    inputEnergy.clear();
    inputTotalStrain.clear();
    inputMaxStrain.clear();
    inputHydrogens.clear();
    atomType.clear();
    atomCharge.clear();
    atomBonds.clear();
    
    bondNum.clear();
    bondStart.clear();
    bondEnd.clear();
    bondType.clear();
    
    xyzCount = -1;
    origXyzCount = -1;
    
    dockNum.clear();
    colorNum.clear();
    
    hydrogensToRotate = 0;
    hydrogenRotAngles.clear();
    dihedrals.clear();
    rotAngles.clear();
    
    bondDists.clear();
    bondDistsOrderKeys.clear();
    
    phase = 0;
}

Mol2::Mol2(const std::string* mol2FileName, 
           const std::string* nameFileName,
           const std::vector<std::string>* mol2text) {
    blankNew();
    
    /* Read from text if provided */
    if (mol2text != nullptr) {
        for (size_t i = 0; i < mol2text->size(); i++) {
            processLine((*mol2text)[i]);
        }
        xyzCount++;
    }
    
    /* Read name file */
    if (nameFileName != nullptr) {
        std::ifstream namefile(*nameFileName);
        if (namefile.is_open()) {
            std::string firstLine;
            std::getline(namefile, firstLine);
            
            std::istringstream iss(firstLine);
            std::vector<std::string> tokens;
            std::string token;
            while (iss >> token) {
                tokens.push_back(token);
            }
            
            if (!tokens.empty()) {
                if (tokens[0] == "name.txt" && tokens.size() >= 6) {
                    name = tokens[2];
                    protName = "none";
                    smiles = tokens[3];
                    longname = tokens[5];
                } else if (tokens[0] == "name.cxcalc.txt" && tokens.size() >= 8) {
                    name = tokens[2];
                    protName = tokens[3];
                    smiles = tokens[4];
                    longname = tokens[7];
                } else if (tokens.size() == 7) {
                    /* Decoys format */
                    name = tokens[1];
                    protName = tokens[2];
                    smiles = tokens[3];
                    longname = tokens[6];
                } else if (tokens.size() == 5) {
                    /* Ligands format */
                    name = tokens[1];
                    protName = "none";
                    smiles = tokens[2];
                    longname = tokens[4];
                } else if (tokens.size() == 3) {
                    /* dbgen format */
                    name = tokens[0];
                    protName = "none";
                    smiles = tokens[1];
                    longname = tokens[2];
                }
            }
            
            namefile.close();
        }
    }
    
    /* Read mol2 file */
    if (mol2FileName != nullptr) {
        /* Check if reading from stdin */
        if (*mol2FileName == "-") {
            std::string line;
            while (std::getline(std::cin, line)) {
                processLine(line);
            }
        } else {
        /* Check if gzipped */
        bool isGzipped = (mol2FileName->size() > 3 &&
                         mol2FileName->substr(mol2FileName->size() - 3) == ".gz");

        if (isGzipped) {
            gzFile gzfile = gzopen(mol2FileName->c_str(), "rb");
            if (gzfile != NULL) {
                char buffer[4096];
                while (gzgets(gzfile, buffer, sizeof(buffer)) != NULL) {
                    processLine(std::string(buffer));
                }
                gzclose(gzfile);
            }
        } else {
            std::ifstream mol2file(*mol2FileName);
            if (mol2file.is_open()) {
                std::string line;
                while (std::getline(mol2file, line)) {
                    processLine(line);
                }
                mol2file.close();
            }
        }
        } // end else (not stdin)
        
        xyzCount++;
        origXyzCount = xyzCount;
        
        /* Ensure inputEnergy, strain, and inputHydrogens have correct size */
        while ((int)inputEnergy.size() < xyzCount) {
            inputEnergy.push_back(9999.99);
        }
        while ((int)inputTotalStrain.size() < xyzCount) {
            inputTotalStrain.push_back(9999.99);
        }
        while ((int)inputMaxStrain.size() < xyzCount) {
            inputMaxStrain.push_back(9999.99);
        }
        while ((int)inputHydrogens.size() < xyzCount) {
            inputHydrogens.push_back(0);
        }
    }
}

void Mol2::processLine(const std::string& line) {
    if (line.substr(0, 17) == "@<TRIPOS>MOLECULE") {
        phase = 1;
    } else if (line.substr(0, 13) == "@<TRIPOS>ATOM") {
        phase = 2;
        xyzCount++;
        atomXyz.push_back(std::vector<std::array<double, 3>>());
    } else if (line.substr(0, 13) == "@<TRIPOS>BOND") {
        phase = 3;
    } else if (line.substr(0, 9) == "@<TRIPOS>") {
        phase = 0;
    } else if (line.length() <= 1 || line[0] == '#') {
        /* Comment or empty line */
    } else {
        if (phase == 1) {
            /* Header phase */
            if (name == "fake") {
                std::istringstream iss(line);
                iss >> name;
            }
            
            /* Check for energy line */
            if (line.find("mmff94s") == 0) {
                std::istringstream iss(line);
                std::string dummy, eq;
                double energy;
                if (iss >> dummy >> eq >> energy) {
                    inputEnergy.push_back(energy);
                    inputHydrogens.push_back(0);
                }
                phase = 0;
            }
            /* Check for strain line (format: "Strain = <totalStrain> <maxStrain>") */
            else if (line.find("Strain") == 0) {
                std::istringstream iss(line);
                std::string dummy, eq;
                double totalStrain, maxStrain;
                if (iss >> dummy >> eq >> totalStrain >> maxStrain) {
                    inputTotalStrain.push_back(totalStrain);
                    inputMaxStrain.push_back(maxStrain);
                    inputHydrogens.push_back(0);
                }
                phase = 0;
            }
        } else if (phase == 2) {
            /* Atom phase */
            std::istringstream iss(line);
            int num;
            std::string atomname, atomtype;
            double x, y, z, charge;
            int resnum;
            std::string resname;
            
            if (xyzCount == 0) {
                /* First conformation - read all atom data */
                if (iss >> num >> atomname >> x >> y >> z >> atomtype >> 
                    resnum >> resname >> charge) {
                    atomNum.push_back(num);
                    atomName.push_back(atomname);
                    atomType.push_back(atomtype);
                    atomCharge.push_back(charge);
                    atomBonds.push_back(std::vector<std::pair<int, std::string>>());
                    
                    atomXyz[xyzCount].push_back(std::array<double, 3>{x, y, z});
                }
            } else {
                /* Later conformations - only read coordinates */
                if (iss >> num >> atomname >> x >> y >> z) {
                    atomXyz[xyzCount].push_back(std::array<double, 3>{x, y, z});
                }
            }
        } else if (phase == 3 && xyzCount == 0) {
            /* Bond phase - only read for first molecule */
            std::istringstream iss(line);
            int num, start, end;
            std::string type;
            
            if (iss >> num >> start >> end >> type) {
                bondNum.push_back(num);
                bondStart.push_back(start);
                bondEnd.push_back(end);
                bondType.push_back(type);
                
                atomBonds[start - 1].push_back(std::make_pair(end - 1, type));
                atomBonds[end - 1].push_back(std::make_pair(start - 1, type));
            }
        }
    }
}

Mol2 Mol2::copy() const {
//Mol2 copy() const {
    Mol2 newM;
    newM.name = name;
    newM.protName = protName;
    newM.smiles = smiles;
    newM.longname = longname;
    
    newM.atomNum = atomNum;
    newM.atomName = atomName;
    newM.atomType = atomType;
    newM.atomCharge = atomCharge;
    newM.atomBonds = atomBonds;
    
    newM.bondNum = bondNum;
    newM.bondStart = bondStart;
    newM.bondEnd = bondEnd;
    newM.bondType = bondType;
    
    /* These are copied for manipulation */
    newM.atomXyz = atomXyz;
    newM.inputEnergy      = inputEnergy;
    newM.inputTotalStrain = inputTotalStrain;
    newM.inputMaxStrain   = inputMaxStrain;
    newM.inputHydrogens   = inputHydrogens;
    newM.xyzCount = xyzCount;
    newM.origXyzCount = origXyzCount;
    
    newM.dockNum = dockNum;
    newM.colorNum = colorNum;
    if (colorConverter) {
        newM.colorConverter.reset(new ColorTable(*colorConverter));
    }
    
    newM.solvData = solvData;
    
    newM.hydrogenRotAngles = hydrogenRotAngles;
    newM.dihedrals = dihedrals;
    newM.rotAngles = rotAngles;
    newM.hydrogensToRotate = hydrogensToRotate;
    
    /* Don't copy cached bond distances */
    
    return newM;
}

// Continued in Part 3...
// mol2.cpp Part 3 - Bond queries and operations

void Mol2::keepConfsOnly(int first, int last) {
    if (last > xyzCount) {
        last = xyzCount;
    }
    
    xyzCount = last - first;
    
    std::vector<std::vector<std::array<double, 3>>> newXyz;
    std::vector<double> newEnergy;
    std::vector<double> newTotalStrain;
    std::vector<double> newMaxStrain;
    std::vector<int> newHydrogens;

    for (int i = first; i < last; i++) {
        newXyz.push_back(atomXyz[i]);
        newEnergy.push_back(i < (int)inputEnergy.size() ? inputEnergy[i] : 9999.99);
        newTotalStrain.push_back(i < (int)inputTotalStrain.size() ? inputTotalStrain[i] : 9999.99);
        newMaxStrain.push_back(i < (int)inputMaxStrain.size() ? inputMaxStrain[i] : 9999.99);
        newHydrogens.push_back(i < (int)inputHydrogens.size() ? inputHydrogens[i] : 0);
    }

    atomXyz = newXyz;
    inputEnergy      = newEnergy;
    inputTotalStrain = newTotalStrain;
    inputMaxStrain   = newMaxStrain;
    inputHydrogens   = newHydrogens;
}

int Mol2::bondsBetween(int atomNum, int atomOther) const {
    return bondsBetweenActual(atomNum - 1, atomOther - 1);
}

int Mol2::bondsBetweenActual(int actualNum, int actualOther) const {
    if (bondDists.empty()) {
        const_cast<Mol2*>(this)->calcBondDists();
    }
    
    int row = bondDistsOrderKeys.at(actualNum);
    int col = bondDistsOrderKeys.at(actualOther);
    return bondDists[row][col];
}

void Mol2::calcBondDists() {
    std::unordered_map<int, std::vector<std::pair<int, int>>> neighbors;
    
    for (size_t i = 0; i < atomNum.size(); i++) {
        neighbors[i] = std::vector<std::pair<int, int>>();
        for (size_t j = 0; j < atomBonds[i].size(); j++) {
            neighbors[i].push_back(std::make_pair(atomBonds[i][j].first, 1));
        }
    }
    
    floydWarshall(neighbors, &bondDists, &bondDistsOrderKeys);
}

std::vector<int> Mol2::distFromAtoms(const std::vector<int>& atoms) const {
    std::unordered_map<int, std::vector<std::pair<int, int>>> neighbors;
    std::vector<int> nodes;
    
    for (size_t i = 0; i < atomNum.size(); i++) {
        nodes.push_back(i);
        neighbors[i] = std::vector<std::pair<int, int>>();
        for (size_t j = 0; j < atomBonds[i].size(); j++) {
            neighbors[i].push_back(std::make_pair(atomBonds[i][j].first, 1));
        }
    }
    
    std::unordered_map<int, int> dists = shortestPaths(nodes, neighbors, 0, atoms);
    
    std::vector<int> result(atomNum.size());
    for (size_t i = 0; i < atomNum.size(); i++) {
        auto it = dists.find(i);
        if (it != dists.end()) {
            result[i] = it->second;
        } else {
            result[i] = 999999999;
        }
    }
    
    return result;
}

bool Mol2::bondedTo(int atomNum, const std::string& firstName, 
                    int bondsAway, const std::string* lastBond,
                    bool returnAtom, int* returnedAtomNum) const {
    return bondedToActual(atomNum - 1, firstName, bondsAway, lastBond, 
                         returnAtom, returnedAtomNum);
}

std::map<int, std::vector<int>> Mol2::bondedToAll(
    int atomNum, const std::string& firstName,
    int bondsAway, const std::string* lastBond) const {
    return bondedToActualAll(atomNum - 1, firstName, bondsAway, lastBond);
}

bool Mol2::bondedToActual(int actualNum, const std::string& firstName,
                           int bondsAway, const std::string* lastBond,
                           bool returnAtom, int* returnedAtomNum) const {
    std::map<int, std::vector<int>> bondedAwayNums = 
        bondedToActualAll(actualNum, firstName, bondsAway, lastBond);
    
    auto it = bondedAwayNums.find(bondsAway);
    if (it == bondedAwayNums.end()) {
        if (returnAtom && returnedAtomNum) {
            *returnedAtomNum = -1;
        }
        return false;
    }
    
    for (size_t i = 0; i < it->second.size(); i++) {
        int anAtomNum = it->second[i];
        if (firstName.empty() || atomType[anAtomNum].find(firstName) == 0) {
            if (returnAtom && returnedAtomNum) {
                *returnedAtomNum = atomNum[anAtomNum];
            }
            return true;
        }
    }
    
    if (returnAtom && returnedAtomNum) {
        *returnedAtomNum = -1;
    }
    return false;
}

std::map<int, std::vector<int>> Mol2::bondedToActualAll(
    int actualNum, const std::string& firstName,
    int bondsAway, const std::string* lastBond) const {
    
    std::map<int, std::vector<int>> bondedAwayNums;
    bondedAwayNums[0].push_back(actualNum);
    
    int checked = 0;
    while (checked < bondsAway) {
        for (size_t i = 0; i < bondedAwayNums[checked].size(); i++) {
            int startNum = bondedAwayNums[checked][i];
            
            for (size_t j = 0; j < atomBonds[startNum].size(); j++) {
                int otherAtom = atomBonds[startNum][j].first;
                std::string bondType = atomBonds[startNum][j].second;
                
                /* Check if already in a previous list */
                bool okayToAdd = true;
                for (auto it = bondedAwayNums.begin(); 
                     it != bondedAwayNums.end(); ++it) {
                    for (size_t k = 0; k < it->second.size(); k++) {
                        if (it->second[k] == otherAtom) {
                            okayToAdd = false;
                            break;
                        }
                    }
                    if (!okayToAdd) break;
                }
                
                /* Check bond type constraint */
                if (okayToAdd && lastBond != nullptr && 
                    checked + 1 == bondsAway) {
                    if (*lastBond != "*") {
                        if (bondType.find(*lastBond) == std::string::npos) {
                            okayToAdd = false;
                        }
                    }
                }
                
                if (okayToAdd) {
                    bondedAwayNums[checked + 1].push_back(otherAtom);
                }
            }
        }
        checked++;
    }
    
    return bondedAwayNums;
}

bool Mol2::isAtomBondedOtherThan(int atomNum, const std::set<int>& count,
                                 const std::set<std::string>& otherThan) const {
    std::map<int, std::vector<int>> bondedAwayNums = 
        bondedToAll(atomNum, "", 1, nullptr);
    
    int otherThanCount = 0;
    auto it = bondedAwayNums.find(1);
    if (it != bondedAwayNums.end()) {
        for (size_t i = 0; i < it->second.size(); i++) {
            int atomNumActual = it->second[i];
            if (otherThan.find(atomType[atomNumActual]) == otherThan.end()) {
                otherThanCount++;
            }
        }
    }
    
    return (count.find(otherThanCount) == count.end());
}

void Mol2::convertDockTypes(const std::string* parameterFileName) {
    AtomConverter converter(parameterFileName);
    dockNum.clear();
    
    for (size_t i = 0; i < atomNum.size(); i++) {
        dockNum.push_back(converter.convertMol2atomNum(*this, atomNum[i]));
    }
}

void Mol2::addColors(const std::string* parameterFileName) {
    colorConverter.reset(new ColorTable(parameterFileName));
    colorNum.clear();
    
    for (size_t i = 0; i < atomNum.size(); i++) {
        colorNum.push_back(colorConverter->convertMol2color(*this, atomNum[i]));
    }
}

std::array<double, 3> Mol2::getXyz(int xyzCount, int atomNum) const {
    for (size_t i = 0; i < this->atomNum.size(); i++) {
        if (this->atomNum[i] == atomNum) {
            return atomXyz[xyzCount][i];
        }
    }
    return std::array<double, 3>{0.0, 0.0, 0.0};
}

std::vector<std::array<double, 3>> Mol2::getXyzManyConfs(
    const std::vector<int>& xyzCounts, int atomIndex) const {
    
    std::vector<std::array<double, 3>> xyzConfs;
    for (size_t i = 0; i < xyzCounts.size(); i++) {
        xyzConfs.push_back(atomXyz[xyzCounts[i]][atomIndex]);
    }
    return xyzConfs;
}

void Mol2::addSolvDataPartialCharges(const std::vector<double>& partialCharges) {
    for (size_t i = 0; i < atomCharge.size() && i < partialCharges.size(); i++) {
        atomCharge[i] = partialCharges[i];
    }
}

// Continued in Part 4...
// mol2.cpp Part 4 - RMSD, clustering, and file I/O

double Mol2::getRMSD(int xyzOne, int xyzTwo) const {
    double sumSquared = 0.0;
    
    for (size_t i = 0; i < atomXyz[xyzOne].size(); i++) {
        sumSquared += distL2Squared3(atomXyz[xyzOne][i], atomXyz[xyzTwo][i]);
    }
    
    double rmsd = sqrt(sumSquared / atomXyz[xyzOne].size());
    return rmsd;
}

std::map<int, std::map<int, double>> Mol2::getRMSDtable(bool forceRedo) {
    if (rmsdTable.empty() || forceRedo) {
        rmsdTable.clear();
        rmsdList.clear();
        
        for (int i = 0; i < xyzCount; i++) {
            rmsdTable[i] = std::map<int, double>();
        }
        
        for (int i = 0; i < xyzCount; i++) {
            for (int j = i + 1; j < xyzCount; j++) {
                double rmsd = getRMSD(i, j);
                rmsdTable[i][j] = rmsd;
                rmsdTable[j][i] = rmsd;
                rmsdList.push_back(std::make_tuple(rmsd, i, j));
            }
        }
        
        std::sort(rmsdList.begin(), rmsdList.end(),
                 [](const auto& a, const auto& b) { 
                     return std::get<0>(a) < std::get<0>(b); 
                 });
    }
    
    return rmsdTable;
}

std::vector<std::tuple<double, int, int>> Mol2::getRMSDlist() {
    getRMSDtable();
    return rmsdList;
}

std::vector<std::vector<int>> Mol2::getRMSDclusters(
    const double* rmsdCutoff, int numClusters) {
    
    getRMSDtable();
    
    UnionFind clusters;
    for (int i = 0; i < xyzCount; i++) {
        clusters.find(i);
    }
    
    double cutoff = (rmsdCutoff != nullptr) ? *rmsdCutoff : 
                    (rmsdList.empty() ? 0.0 : std::get<0>(rmsdList.back()) + 1.0);
    
    for (size_t i = 0; i < rmsdList.size(); i++) {
        if (std::get<0>(rmsdList[i]) > cutoff) {
            break;
        }
        clusters.unionSets(std::get<1>(rmsdList[i]), std::get<2>(rmsdList[i]));
    }
    
    return clusters.toLists();
}

std::vector<std::vector<int>> Mol2::getRMSDclustersAll(
    const double* rmsdCutoff, int numClusters) {
    
    getRMSDtable();
    
    UnionFind clusters;
    for (int i = 0; i < xyzCount; i++) {
        clusters.find(i);
    }
    
    double cutoff = (rmsdCutoff != nullptr) ? *rmsdCutoff : 
                    (rmsdList.empty() ? 0.0 : std::get<0>(rmsdList.back()) + 1.0);
    
    for (size_t i = 0; i < rmsdList.size(); i++) {
        double thisRMSD = std::get<0>(rmsdList[i]);
        if (thisRMSD > cutoff) {
            break;
        }
        
        int conf1 = std::get<1>(rmsdList[i]);
        int conf2 = std::get<2>(rmsdList[i]);
        
        if (clusters.different(conf1, conf2)) {
            /* All linkage - check all pairs */
            bool combine = true;
            std::vector<int> cluster1 = clusters.getList(conf1);
            std::vector<int> cluster2 = clusters.getList(conf2);
            
            for (size_t j = 0; j < cluster1.size() && combine; j++) {
                for (size_t k = 0; k < cluster2.size() && combine; k++) {
                    if (rmsdTable[cluster1[j]][cluster2[k]] > thisRMSD) {
                        combine = false;
                    }
                }
            }
            
            if (combine) {
                clusters.unionSets(conf1, conf2);
            }
        }
    }
    
    return clusters.toLists();
}

std::vector<std::vector<int>> Mol2::divisiveClustering() {
    int numClusters = std::min(20, std::max(1, origXyzCount / 3));
    return mol2db2::divisiveClustering(atomXyz, numClusters);
}

void Mol2::writeMol2File(std::ostream& outFile, 
                         const std::vector<int>* whichXyz) const {
    std::vector<int> xyzToWrite;
    if (whichXyz != nullptr) {
        xyzToWrite = *whichXyz;
    } else {
        for (int i = 0; i < xyzCount; i++) {
            xyzToWrite.push_back(i);
        }
    }
    
    for (size_t idx = 0; idx < xyzToWrite.size(); idx++) {
        int oneXyz = xyzToWrite[idx];
        
        outFile << "@<TRIPOS>MOLECULE\n";
        if (protName != "fake") {
            outFile << name << " " << protName << "\n";
        } else {
            outFile << name << "\n";
        }
        
        outFile << atomNum.size() << " " << bondNum.size() 
                << " 0 0 0\n";
        outFile << "SMALL\nUSER_CHARGES\n\n";
        
        if (oneXyz < (int)inputEnergy.size()) {
            outFile << "mmff94s_NoEstat = " << inputEnergy[oneXyz] << "\n";
        }
        
        outFile << "@<TRIPOS>ATOM\n";
        for (size_t i = 0; i < atomNum.size(); i++) {
            outFile << atomNum[i] << " " << atomName[i] << "    "
                   << atomXyz[oneXyz][i][0] << "  "
                   << atomXyz[oneXyz][i][1] << "  "
                   << atomXyz[oneXyz][i][2] << " "
                   << atomType[i] << "     1 <0>       "
                   << atomCharge[i] << "\n";
        }
        
        outFile << "@<TRIPOS>BOND\n";
        for (size_t i = 0; i < bondNum.size(); i++) {
            outFile << bondNum[i] << " " << bondStart[i] << " "
                   << bondEnd[i] << " " << bondType[i] << "\n";
        }
    }
}

void Mol2::writeMol2(const std::string& outName, 
                     const std::vector<int>* whichXyz) const {
    std::ofstream outFile(outName);
    if (outFile.is_open()) {
        writeMol2File(outFile, whichXyz);
        outFile.close();
    } else {
        std::cerr << "Error opening output file: " << outName << std::endl;
    }
}

void Mol2::deleteBond(int bondInd) {
    int nbonds = bondNum.size();
    
    /* Update bond numbers */
    for (int i = bondInd + 1; i < nbonds; i++) {
        bondNum[i]--;
    }
    
    /* Delete the bond */
    bondNum.erase(bondNum.begin() + bondInd);
    bondStart.erase(bondStart.begin() + bondInd);
    bondEnd.erase(bondEnd.begin() + bondInd);
    bondType.erase(bondType.begin() + bondInd);
}

void Mol2::deleteAtom(int atomInd) {
    int natoms = atomNum.size();
    int atomn = atomNum[atomInd];
    
    /* Update atom numbers */
    for (int i = atomInd + 1; i < natoms; i++) {
        atomNum[i]--;
    }
    
    /* Remove atom from tables */
    atomNum.erase(atomNum.begin() + atomInd);
    atomType.erase(atomType.begin() + atomInd);
    atomBonds.erase(atomBonds.begin() + atomInd);
    atomName.erase(atomName.begin() + atomInd);
    atomCharge.erase(atomCharge.begin() + atomInd);
    
    if (!dockNum.empty()) {
        dockNum.erase(dockNum.begin() + atomInd);
    }
    if (!colorNum.empty()) {
        colorNum.erase(colorNum.begin() + atomInd);
    }
    
    for (size_t i = 0; i < atomXyz.size(); i++) {
        atomXyz[i].erase(atomXyz[i].begin() + atomInd);
    }
    
    /* Delete bonds containing this atom */
    int nbonds = bondNum.size();
    for (int i = 0; i < nbonds; ) {
        if (bondStart[i] == atomn || bondEnd[i] == atomn) {
            deleteBond(i);
            nbonds--;
        } else {
            i++;
        }
    }
    
    /* Update bond atom numbering */
    for (size_t i = 0; i < bondStart.size(); i++) {
        if (bondStart[i] >= atomn) bondStart[i]--;
        if (bondEnd[i] >= atomn) bondEnd[i]--;
    }
    
    /* Update atomBonds numbering */
    for (size_t i = 0; i < atomBonds.size(); i++) {
        for (size_t j = 0; j < atomBonds[i].size(); j++) {
            if (atomBonds[i][j].first >= atomInd) {
                atomBonds[i][j].first--;
            }
        }
    }
}

std::tuple<std::string, std::vector<int>> Mol2::removeCovalentDummyAtom() {
    std::vector<int> indicesList;
    int si_index = -1;
    std::string covAtomType;
    
    /* Find Si atom */
    for (size_t i = 0; i < atomType.size(); i++) {
        if (atomType[i] == "Si") {
            si_index = i;
            atomType[i] = "del";
            break;
        }
    }
    
    if (si_index == -1) {
        return std::make_tuple("", indicesList);
    }
    
    /* Find hydrogens and covalent attachment point */
    for (size_t i = 0; i < atomBonds[si_index].size(); i++) {
        int bondedAtom = atomBonds[si_index][i].first;
        if (atomType[bondedAtom] == "H") {
            atomType[bondedAtom] = "del";
        } else {
            covAtomType = atomType[bondedAtom];
            atomType[bondedAtom] = "cov";
            
            /* Remove bond to Si */
            for (size_t j = 0; j < atomBonds[bondedAtom].size(); ) {
                if (atomBonds[bondedAtom][j].first == si_index) {
                    atomBonds[bondedAtom].erase(atomBonds[bondedAtom].begin() + j);
                } else {
                    j++;
                }
            }
        }
    }
    
    /* Delete marked atoms */
    int natoms = atomNum.size();
    for (int i = 0; i < natoms; ) {
        if (atomType[i] == "del") {
            indicesList.push_back(i);
            deleteAtom(i);
            natoms--;
        } else {
            i++;
        }
    }
    
    return std::make_tuple(covAtomType, indicesList);
}

void Mol2::recolorCovalentAttachment(const std::string& covAtomType) {
    int cov_index = -1;
    
    /* Find covalent attachment atom */
    for (size_t i = 0; i < atomType.size(); i++) {
        if (atomType[i] == "cov") {
            cov_index = i;
            break;
        }
    }
    
    if (cov_index == -1) {
        return;
    }
    
    /* Color it and neighbors */
    if (cov_index < (int)colorNum.size()) {
        colorNum[cov_index] = 8;
        
        int color = 9;
        for (size_t i = 0; i < atomBonds[cov_index].size(); i++) {
            int neighbor = atomBonds[cov_index][i].first;
            if (neighbor < (int)colorNum.size()) {
                colorNum[neighbor] = color;
                color++;
            }
        }
    }
    
    /* Reset atom type */
    atomType[cov_index] = covAtomType;
}

// Utility function for reading DOCK mol2 files
std::tuple<std::vector<Mol2>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>>
readDockMol2file(const std::string& mol2FileName, bool recdes,
                 bool ligdes, bool charge, bool elec) {
    
    std::vector<Mol2> mol2data;
    std::vector<double> mol2rd, mol2ld, mol2charge, mol2elec;
    std::vector<std::string> mol2lines;
    std::vector<std::string> mol2names;
    
    std::ifstream mol2file(mol2FileName);
    if (!mol2file.is_open()) {
        return std::make_tuple(mol2data, mol2rd, mol2ld, mol2charge, mol2elec);
    }
    
    std::string line;
    while (std::getline(mol2file, line)) {
        if (line.substr(0, 17) == "@<TRIPOS>MOLECULE") {
// mol2.cpp Part 5 - Complete readDockMol2file function
// (This continues from Part 4)
            if (!mol2lines.empty()) {
                Mol2 newMol2(nullptr, nullptr, &mol2lines);
                mol2data.push_back(newMol2);
            }
            mol2lines.clear();
        }
        
        if (!line.empty() && line[0] != '#') {
            mol2lines.push_back(line);
        }
        
        /* Parse special comment lines for DOCK output */
        if (line.find("##########                 Name:") == 0) {
            std::istringstream iss(line);
            std::string dummy1, dummy2, name;
            iss >> dummy1 >> dummy2 >> name;
            mol2names.push_back(name);
        }
        
        if (ligdes && line.find("##########  Ligand Polar Desolv:") == 0) {
            std::istringstream iss(line);
            std::string dummy;
            double val;
            iss >> dummy >> dummy >> dummy >> dummy >> val;
            mol2ld.push_back(val);
        }
        
        if (recdes && line.find("########## Receptor Desolvation:") == 0) {
            std::istringstream iss(line);
            std::string dummy;
            double val;
            iss >> dummy >> dummy >> dummy >> val;
            mol2rd.push_back(val);
        }
        
        if (charge && line.find("##########        Ligand Charge:") == 0) {
            std::istringstream iss(line);
            std::string dummy;
            double val;
            iss >> dummy >> dummy >> dummy >> val;
            mol2charge.push_back(val);
        }
        
        if (elec && line.find("##########        Electrostatic:") == 0) {
            std::istringstream iss(line);
            std::string dummy;
            double val;
            iss >> dummy >> dummy >> val;
            mol2elec.push_back(val);
        }
    }
    
    /* Don't forget last molecule */
    if (!mol2lines.empty()) {
        Mol2 newMol2(nullptr, nullptr, &mol2lines);
        mol2data.push_back(newMol2);
    }
    
    mol2file.close();
    
    /* Set names */
    for (size_t i = 0; i < mol2data.size() && i < mol2names.size(); i++) {
        mol2data[i].name = mol2names[i];
    }
    
    return std::make_tuple(mol2data, mol2rd, mol2ld, mol2charge, mol2elec);
}

} // namespace mol2db2
