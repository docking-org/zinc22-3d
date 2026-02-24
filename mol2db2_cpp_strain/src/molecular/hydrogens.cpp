#include "hydrogens.h"
#include "mol2.h"
#include "geometry.h"
#include "combinatorics.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

namespace mol2db2 {

std::vector<Hydrogens::Rule> Hydrogens::getDefaultRules() {
    std::vector<Rule> rules;
    rules.push_back(std::make_tuple(1, "C.ar", "1", "S", "1", "H", "180"));
    rules.push_back(std::make_tuple(1, "C.ar", "1", "O", "1", "H", "180"));
    rules.push_back(std::make_tuple(1, "C.1", "1", "S", "1", "H", "-"));
    rules.push_back(std::make_tuple(1, "C.1", "1", "O", "1", "H", "-"));
    rules.push_back(std::make_tuple(1, "C", "1", "S", "1", "H", "120,240"));
    rules.push_back(std::make_tuple(1, "C", "1", "O", "1", "H", "120,240"));
    rules.push_back(std::make_tuple(2, "2", "2", "N", "1", "H", "-"));
    rules.push_back(std::make_tuple(2, "1", "2", "N", "1", "H", "180"));
    return rules;
}

Hydrogens::Hydrogens(const std::string* parameterFileName) {
    if (parameterFileName != nullptr) {
        std::ifstream parameterFile(*parameterFileName);
        if (parameterFile.is_open()) {
            std::string line;
            while (std::getline(parameterFile, line)) {
                std::istringstream iss(line);
                std::vector<std::string> tokens;
                std::string token;
                while (iss >> token) {
                    tokens.push_back(token);
                }
                
                if (tokens.size() == 7) {
                    rules.push_back(std::make_tuple(
                        std::stoi(tokens[0]), tokens[1], tokens[2],
                        tokens[3], tokens[4], tokens[5], tokens[6]));
                }
            }
            parameterFile.close();
        } else {
            rules = getDefaultRules();
        }
    } else {
        rules = getDefaultRules();
    }
}

void Hydrogens::printParameters() const {
    for (size_t i = 0; i < rules.size(); i++) {
        std::cout << std::get<0>(rules[i]) << " "
                 << std::get<1>(rules[i]) << " "
                 << std::get<2>(rules[i]) << " "
                 << std::get<3>(rules[i]) << " "
                 << std::get<4>(rules[i]) << " "
                 << std::get<5>(rules[i]) << " "
                 << std::get<6>(rules[i]) << std::endl;
    }
}

void Hydrogens::findTerminalHydrogens(Mol2* mol2data) {
    mol2data->hydrogenRotAngles.clear();
    mol2data->hydrogensToRotate = 0;
    mol2data->dihedrals.clear();
    
    for (size_t count = 0; count < mol2data->atomNum.size(); count++) {
        int atomNum = mol2data->atomNum[count];
        const std::string& atomType = mol2data->atomType[count];
        std::string result = "-";
        
        for (size_t ruleIdx = 0; ruleIdx < rules.size(); ruleIdx++) {
            const Rule& rule = rules[ruleIdx];
            int ruleType = std::get<0>(rule);
            
            if (ruleType == 1) {
                /* Type 1 rule: (1, "C.1", "1", "O", "1", "H", "-") */
                std::string atomMatch = std::get<5>(rule);
                std::string atom2 = std::get<3>(rule);
                std::string bond1 = std::get<4>(rule);
                std::string atom1 = std::get<1>(rule);
                std::string bond0 = std::get<2>(rule);
                std::string degrees = std::get<6>(rule);
                
                if (atomType.find(atomMatch) == 0) {
                    if (mol2data->bondedTo(atomNum, atom2, 1, &bond1, 
                                          false, nullptr)) {
                        if (mol2data->bondedTo(atomNum, atom1, 2, &bond0, 
                                              false, nullptr)) {
                            result = degrees;
                            break;
                        }
                    }
                }
            } else if (ruleType == 2) {
                /* Type 2 rule: (2, "2", "2", "N", "1", "H", "-") */
                std::string atomMatch = std::get<5>(rule);
                std::string atom2 = std::get<3>(rule);
                std::string bond1 = std::get<4>(rule);
                std::string bond0a = std::get<1>(rule);
                std::string bond0b = std::get<2>(rule);
                std::string degrees = std::get<6>(rule);
                
                if (atomType.find(atomMatch) == 0) {
                    if (mol2data->bondedTo(atomNum, atom2, 1, &bond1, 
                                          false, nullptr)) {
                        if (mol2data->bondedTo(atomNum, "", 2, &bond0a, 
                                              false, nullptr)) {
                            if (mol2data->bondedTo(atomNum, "", 2, &bond0b, 
                                                  false, nullptr)) {
                                result = degrees;
                                break;
                            }
                        }
                    }
                }
            }
        }
        
        mol2data->hydrogenRotAngles.push_back(result);
        if (result != "-") {
            mol2data->hydrogensToRotate++;
        }
    }
}

void Hydrogens::findDihedrals(Mol2* mol2data) {
    if (!mol2data->dihedrals.empty()) {
        return;  /* Already found */
    }
    
    mol2data->dihedrals.clear();
    mol2data->rotAngles.clear();
    
    for (size_t count = 0; count < mol2data->hydrogenRotAngles.size(); count++) {
        const std::string& angles = mol2data->hydrogenRotAngles[count];
        
        if (angles != "-") {
            std::vector<int> dihedral(4, -1);
            int hydrogenNum = mol2data->atomNum[count];
            
            dihedral[3] = hydrogenNum;
            
            /* Find atoms 2, 1, 0 bonds away */
            int atom2 = -1;
            mol2data->bondedTo(hydrogenNum, "", 1, nullptr, true, &atom2);
            dihedral[2] = atom2;
            
            int atom1 = -1;
            mol2data->bondedTo(hydrogenNum, "", 2, nullptr, true, &atom1);
            dihedral[1] = atom1;
            
            int atom0 = -1;
            mol2data->bondedTo(hydrogenNum, "", 3, nullptr, true, &atom0);
            dihedral[0] = atom0;
            
            mol2data->dihedrals[hydrogenNum] = dihedral;
            mol2data->rotAngles[hydrogenNum] = std::vector<double>();
            
            /* Parse angles */
            std::istringstream iss(angles);
            std::string angleStr;
            while (std::getline(iss, angleStr, ',')) {
                mol2data->rotAngles[hydrogenNum].push_back(std::stod(angleStr));
            }
        }
    }
}

std::tuple<std::vector<std::array<double, 3>>, double> 
Hydrogens::getCurDihedral(int atomNum, int xyzCount, const Mol2& mol2data) {
    const std::vector<int>& curDihedralNums = mol2data.dihedrals.at(atomNum);
    std::vector<std::array<double, 3>> curXyz;
    
    for (size_t i = 0; i < curDihedralNums.size(); i++) {
        curXyz.push_back(mol2data.getXyz(xyzCount, curDihedralNums[i]));
    }
    
    double dihedral1 = getDihedralUnited(curXyz);
    return std::make_tuple(curXyz, dihedral1);
}

void Hydrogens::resetHydrogens(Mol2* mol2data) {
    findDihedrals(mol2data);
    
    for (int xyzCount = 0; xyzCount < mol2data->xyzCount; xyzCount++) {
        for (size_t atomIndex = 0; atomIndex < mol2data->atomNum.size(); atomIndex++) {
            int atomNum = mol2data->atomNum[atomIndex];
            
            if (mol2data->dihedrals.find(atomNum) != mol2data->dihedrals.end()) {
                const std::vector<double>& angles = mol2data->rotAngles[atomNum];
                
                /* Only reset planar (180 degree) hydrogens */
                bool has180 = false;
                for (size_t i = 0; i < angles.size(); i++) {
                    if (fabs(angles[i] - 180.0) < 0.01) {
                        has180 = true;
                        break;
                    }
                }
                
                if (has180) {
                    auto result = getCurDihedral(atomNum, xyzCount, *mol2data);
                    std::vector<std::array<double, 3>> curXyz = std::get<0>(result);
                    double dihedral1 = std::get<1>(result);
                    
                    /* Rotate to 0 */
                    double newTheta = 0.0 - dihedral1;
                    std::array<double, 3> newHxyz = rotateAboutLine(
                        curXyz[1], curXyz[2], curXyz[3], newTheta);
                    
                    mol2data->atomXyz[xyzCount][atomIndex] = newHxyz;
                    mol2data->inputHydrogens[xyzCount] = 1;  /* Reset */
                }
            }
        }
    }
}

Mol2 Hydrogens::rotateHydrogens(const Mol2& mol2data) {
    Mol2* mutableMol2 = const_cast<Mol2*>(&mol2data);
    findDihedrals(mutableMol2);
    
    Mol2 rotatedMol2data = mol2data.copy();
    //mol2data.colorConverter.reset();
    rotatedMol2data.atomXyz.clear();
    rotatedMol2data.inputEnergy.clear();
    rotatedMol2data.inputTotalStrain.clear();
    rotatedMol2data.inputMaxStrain.clear();
    rotatedMol2data.inputHydrogens.clear();
    rotatedMol2data.origXyzCount = mol2data.xyzCount;
    
    for (int xyzCount = 0; xyzCount < mol2data.xyzCount; xyzCount++) {
        std::vector<std::vector<std::array<double, 3>>> indexToCoords;
        std::vector<std::array<double, 3>> originals;
        
        for (size_t atomIndex = 0; atomIndex < mol2data.atomNum.size(); atomIndex++) {
            int atomNum = mol2data.atomNum[atomIndex];
            
            if (mol2data.dihedrals.find(atomNum) != mol2data.dihedrals.end()) {
                std::vector<std::array<double, 3>> thisCoords;
                auto result = getCurDihedral(atomNum, xyzCount, mol2data);
                std::vector<std::array<double, 3>> curXyz = std::get<0>(result);
                double dihedral1 = std::get<1>(result);
                
                thisCoords.push_back(curXyz[3]);  /* Original */
                originals.push_back(curXyz[3]);
                
                const std::vector<double>& angles = mol2data.rotAngles.at(atomNum);
                for (size_t i = 0; i < angles.size(); i++) {
                    double newTheta = angles[i] * M_PI / 180.0;  /* To radians */
                    std::array<double, 3> newHxyz = rotateAboutLine(
                        curXyz[1], curXyz[2], curXyz[3], newTheta);
                    thisCoords.push_back(newHxyz);
                }
                
                indexToCoords.push_back(thisCoords);
            }
        }
        
        /* Generate all combinations */
        std::vector<std::vector<std::array<double, 3>>> combinations = 
            allCombinations(indexToCoords);
        
        for (size_t combIdx = 0; combIdx < combinations.size(); combIdx++) {
            const std::vector<std::array<double, 3>>& aCombo = combinations[combIdx];
            
            /* Copy original coordinates */
            std::vector<std::array<double, 3>> atomXyzCopy = 
                mol2data.atomXyz[xyzCount];
            
            /* Replace hydrogen positions */
            int replaceIndex = 0;
            for (size_t atomIndex = 0; atomIndex < mol2data.atomNum.size(); atomIndex++) {
                int atomNum = mol2data.atomNum[atomIndex];
                if (mol2data.dihedrals.find(atomNum) != mol2data.dihedrals.end()) {
                    atomXyzCopy[atomIndex] = aCombo[replaceIndex];
                    replaceIndex++;
                }
            }
            
            /* Determine if this is original or rotated */
            bool isOriginal = (aCombo == originals);
            
            rotatedMol2data.atomXyz.push_back(atomXyzCopy);
            rotatedMol2data.inputEnergy.push_back(mol2data.inputEnergy[xyzCount]);
            rotatedMol2data.inputTotalStrain.push_back(mol2data.inputTotalStrain[xyzCount]);
            rotatedMol2data.inputMaxStrain.push_back(mol2data.inputMaxStrain[xyzCount]);

            if (isOriginal) {
                rotatedMol2data.inputHydrogens.push_back(
                    mol2data.inputHydrogens[xyzCount]);
            } else {
                rotatedMol2data.inputHydrogens.push_back(2);  /* Rotated */
            }
        }
    }
    
    rotatedMol2data.xyzCount = rotatedMol2data.atomXyz.size();
    return rotatedMol2data;
}

void Hydrogens::rotateHydrogens_teb(Mol2& mol2data) {
    Mol2* mutableMol2 = const_cast<Mol2*>(&mol2data);
    findDihedrals(mutableMol2);

    Mol2 orimol2data = mol2data.copy();
    
    //mol2data.colorConverter.reset();
    mol2data.atomXyz.clear();
    mol2data.inputEnergy.clear();
    mol2data.inputTotalStrain.clear();
    mol2data.inputMaxStrain.clear();
    mol2data.inputHydrogens.clear();
    //mol2data.origXyzCount = mol2data.xyzCount;
    mol2data.origXyzCount = orimol2data.xyzCount;

    for (int xyzCount = 0; xyzCount < orimol2data.xyzCount; xyzCount++) {
        std::vector<std::vector<std::array<double, 3>>> indexToCoords;
        std::vector<std::array<double, 3>> originals;

        for (size_t atomIndex = 0; atomIndex < orimol2data.atomNum.size(); atomIndex++) {
            int atomNum = mol2data.atomNum[atomIndex];

            if (orimol2data.dihedrals.find(atomNum) != orimol2data.dihedrals.end()) {
                std::vector<std::array<double, 3>> thisCoords;
                auto result = getCurDihedral(atomNum, xyzCount, orimol2data);
                std::vector<std::array<double, 3>> curXyz = std::get<0>(result);
                double dihedral1 = std::get<1>(result);

                thisCoords.push_back(curXyz[3]);  /* Original */
                originals.push_back(curXyz[3]);

                const std::vector<double>& angles = orimol2data.rotAngles.at(atomNum);
                for (size_t i = 0; i < angles.size(); i++) {
                    double newTheta = angles[i] * M_PI / 180.0;  /* To radians */
                    std::array<double, 3> newHxyz = rotateAboutLine(
                        curXyz[1], curXyz[2], curXyz[3], newTheta);
                    thisCoords.push_back(newHxyz);
                }

                indexToCoords.push_back(thisCoords);
            }
        }

        /* Generate all combinations */
        std::vector<std::vector<std::array<double, 3>>> combinations =
            allCombinations(indexToCoords);

        for (size_t combIdx = 0; combIdx < combinations.size(); combIdx++) {
            const std::vector<std::array<double, 3>>& aCombo = combinations[combIdx];

            /* Copy original coordinates */
            std::vector<std::array<double, 3>> atomXyzCopy =
                orimol2data.atomXyz[xyzCount];

            /* Replace hydrogen positions */
            int replaceIndex = 0;
            for (size_t atomIndex = 0; atomIndex < orimol2data.atomNum.size(); atomIndex++) {
                int atomNum = orimol2data.atomNum[atomIndex];
                if (orimol2data.dihedrals.find(atomNum) != orimol2data.dihedrals.end()) {
                    atomXyzCopy[atomIndex] = aCombo[replaceIndex];
                    replaceIndex++;
                }
            }

            /* Determine if this is original or rotated */
            bool isOriginal = (aCombo == originals);

            mol2data.atomXyz.push_back(atomXyzCopy);
            mol2data.inputEnergy.push_back(orimol2data.inputEnergy[xyzCount]);
            mol2data.inputTotalStrain.push_back(orimol2data.inputTotalStrain[xyzCount]);
            mol2data.inputMaxStrain.push_back(orimol2data.inputMaxStrain[xyzCount]);

            if (isOriginal) {
                mol2data.inputHydrogens.push_back(
                    orimol2data.inputHydrogens[xyzCount]);
            } else {
                mol2data.inputHydrogens.push_back(2);  /* Rotated */
            }
        }
    }

    mol2data.xyzCount = mol2data.atomXyz.size();
    return;
}


void Hydrogens::findAngles(const Mol2& mol2data) {
    Mol2* mutableMol2 = const_cast<Mol2*>(&mol2data);
    findDihedrals(mutableMol2);
    
    for (int xyzCount = 0; xyzCount < mol2data.xyzCount; xyzCount++) {
        for (size_t atomIndex = 0; atomIndex < mol2data.atomNum.size(); atomIndex++) {
            int atomNum = mol2data.atomNum[atomIndex];
            
            if (mol2data.dihedrals.find(atomNum) != mol2data.dihedrals.end()) {
                auto result = getCurDihedral(atomNum, xyzCount, mol2data);
                double dihedral1 = std::get<1>(result);
                std::cout << dihedral1 << std::endl;
            }
        }
    }
}

} // namespace mol2db2
