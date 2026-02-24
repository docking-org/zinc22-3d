#include "clash.h"
#include "mol2.h"
#include "geometry.h"
#include <fstream>
#include <sstream>
#include <iostream>

namespace mol2db2 {

std::vector<Clash::Rule> Clash::getDefaultRules() {
    std::vector<Rule> rules;
    /* Only H-H with minimum distance 1.70 at 2+ bonds away */
    rules.push_back(std::make_tuple("min", 2, 1, "H", "H", 1.70, 1.70*1.70));
    return rules;
}

Clash::Clash(const std::string* parameterFileName) {
    if (parameterFileName != nullptr) {
        std::ifstream parameterFile(*parameterFileName);
        if (!parameterFile.is_open()) {
            std::cerr << "Warning: Could not open clash parameter file: " 
                     << *parameterFileName << std::endl;
            rules = getDefaultRules();
            return;
        }
        
        std::string line;
        while (std::getline(parameterFile, line)) {
            std::istringstream iss(line);
            std::string minmax, atomType1, atomType2;
            int bond, cmpVal;
            double dist;
            
            if (iss >> minmax >> bond >> cmpVal >> atomType1 >> atomType2 >> dist) {
                rules.push_back(std::make_tuple(minmax, bond, cmpVal, 
                                                atomType1, atomType2, 
                                                dist, dist * dist));
            }
        }
        
        parameterFile.close();
    } else {
        rules = getDefaultRules();
    }
}

void Clash::printParameters() const {
    for (size_t i = 0; i < rules.size(); i++) {
        std::cout << std::get<0>(rules[i]) << " "
                 << std::get<1>(rules[i]) << " "
                 << std::get<2>(rules[i]) << " "
                 << std::get<3>(rules[i]) << " "
                 << std::get<4>(rules[i]) << " "
                 << std::get<5>(rules[i]) << std::endl;
    }
}

bool Clash::decideDistanceRules(const Mol2& mol2data,
                                const std::vector<std::array<double, 3>>& xyzData) const {
    int natoms = xyzData.size();

    for (const Rule& rule : rules) {
        const std::string& minmax  = std::get<0>(rule);
        int   bondReq = std::get<1>(rule);
        int   cmpReq  = std::get<2>(rule);
        const std::string& typeA   = std::get<3>(rule);
        const std::string& typeB   = std::get<4>(rule);
        double distSq = std::get<6>(rule);

        /* operator selected from cmpReq: -1 → lt, 0 → eq, 1 → gt */
        auto bondCmp = [&](int bd) -> bool {
            if (cmpReq == -1) return bd < bondReq;
            if (cmpReq ==  0) return bd == bondReq;
            return bd > bondReq;
        };

        if (typeA == typeB) {
            /* Collect atoms matching typeA */
            std::vector<int> atoms;
            for (int i = 0; i < natoms; i++) {
                if (typeA == "*" || mol2data.atomType[i].find(typeA) == 0)
                    atoms.push_back(i);
            }
            int n = atoms.size();
            for (int ai = 0; ai < n; ai++) {
                for (int aj = ai + 1; aj < n; aj++) {
                    int atomA = atoms[ai], atomB = atoms[aj];
                    if (!bondCmp(mol2data.bondsBetweenActual(atomA, atomB))) continue;
                    double d2 = distL2Squared3(xyzData[atomA], xyzData[atomB]);
                    if ((minmax == "min" && d2 < distSq) ||
                        (minmax == "max" && d2 > distSq))
                        return true;
                }
            }
        } else {
            /* Different type pairs — old pairwise loop */
            for (int atomA = 0; atomA < natoms; atomA++) {
                if (typeA != "*" && mol2data.atomType[atomA].find(typeA) != 0) continue;
                for (int atomB = 0; atomB < natoms; atomB++) {
                    if (atomA == atomB) continue;
                    if (typeB != "*" && mol2data.atomType[atomB].find(typeB) != 0) continue;
                    if (!bondCmp(mol2data.bondsBetweenActual(atomA, atomB))) continue;
                    double d2 = distL2Squared3(xyzData[atomA], xyzData[atomB]);
                    if ((minmax == "min" && d2 < distSq) ||
                        (minmax == "max" && d2 > distSq))
                        return true;
                }
            }
        }
    }
    return false;
}

} // namespace mol2db2
