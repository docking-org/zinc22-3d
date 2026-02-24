#include "hierarchy.h"
#include "mol2.h"
#include "clash.h"
#include "geometry.h"
#include "buckets2.h"
#include "unionfind.h"
#include <zlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <functional>

namespace mol2db2 {

int computeBreaks(const TooBigError& limitError,
                  long limitset, long limitconf, long limitcoord) {
    int breaksS = (limitset > 0) ? 
        (int)std::ceil(limitError.getSets() / (double)limitset) : 1;
    int breaksC = (limitconf > 0) ? 
        (int)std::ceil(limitError.getConfs() / (double)limitconf) : 1;
    int breaksX = (limitcoord > 0) ? 
        (int)std::ceil(limitError.getCoords() / (double)limitcoord) : 1;
    
    return std::max({breaksS, breaksC, breaksX});
}

Hierarchy::Hierarchy(Mol2& mol2data,
                     const Clash& clashDecider,
                     double tolerance,
                     bool verbose,
                     bool timeit,
                     long limitset,
                     long limitconf,
                     long limitcoord) {

    this->mol2data   = &mol2data;
    this->numConfs   = 0;
    this->rigidStructureCount = 0;

    /* Check size limits */
    int totalCoords = (int)mol2data.atomXyz.size() *
                      (mol2data.atomXyz.empty() ? 0 : (int)mol2data.atomXyz[0].size());

    if (verbose) {
        std::cout << "Total number of sets (complete confs): "
                  << mol2data.atomXyz.size() << std::endl;
    }

    if ((long)mol2data.atomXyz.size() > limitset)
        throw TooBigError(0, mol2data.atomXyz.size(), totalCoords);
    if ((long)totalCoords > limitcoord)
        throw TooBigError(0, mol2data.atomXyz.size(), totalCoords);

    /* Step 1: identify rigid structures (cycles + rigid bonds) */
    if (verbose) std::cout << "Finding rigid structures..." << std::endl;
    int natoms = mol2data.atomXyz.empty() ? 0 : (int)mol2data.atomXyz[0].size();
    getRigidStructures(natoms, mol2data.atomBonds, verbose);

    /* Step 2: cluster positions per atom and build confs */
    if (verbose) std::cout << "Counting positions / building conformations..." << std::endl;
    countPositionsPython(mol2data.atomXyz, tolerance, verbose);

    /* Step 3: find rigid heavy atoms (needed for R-lines) */
    findRigidComponent(mol2data.atomBonds);
    findRigidHeavy(mol2data.atomType);

    if (verbose)
        std::cout << "Total conformations: " << numConfs << std::endl;

    if ((long)numConfs > limitconf)
        throw TooBigError(numConfs, mol2data.atomXyz.size(), totalCoords);

    /* Step 4: clash detection */
    if (verbose) std::cout << "Identifying clashes..." << std::endl;
    identifyClashSetnumsPython(clashDecider, mol2data);
}

void Hierarchy::findRigidComponent(
    const std::vector<std::vector<std::pair<int, std::string>>>& atomBonds) {
    
    UnionFind clusters;
    
    for (size_t atomNum = 0; atomNum < posCount.size(); atomNum++) {
        if (posCount[atomNum] == 1) {
            for (size_t i = 0; i < atomBonds[atomNum].size(); i++) {
                int otherNum = atomBonds[atomNum][i].first;
                if (otherNum < (int)posCount.size() && posCount[otherNum] == 1) {
                    clusters.unionSets(atomNum, otherNum);
                }
            }
        }
    }
    
    std::vector<std::vector<int>> clusterLists = clusters.toLists();
    
    /* Find largest cluster */
    int maxSize = 0;
    int maxIdx = -1;
    for (size_t i = 0; i < clusterLists.size(); i++) {
        if ((int)clusterLists[i].size() > maxSize) {
            maxSize = clusterLists[i].size();
            maxIdx = i;
        }
    }
    
    if (maxIdx >= 0) {
        rigidComponent = clusterLists[maxIdx];
    }
    
    atomsAssigned.clear();
    for (size_t i = 0; i < rigidComponent.size(); i++) {
        atomsAssigned.insert(rigidComponent[i]);
    }
    
    atomsNotAssigned.clear();
    for (size_t i = 0; i < posCount.size(); i++) {
        if (atomsAssigned.find(i) == atomsAssigned.end()) {
            atomsNotAssigned.insert(i);
        }
    }
}

void Hierarchy::findRigidHeavy(const std::vector<std::string>& atomTypes) {
    heavyRigidCount = 0;
    heavyRigidAtomNums.clear();
    
    for (size_t i = 0; i < rigidComponent.size(); i++) {
        int atomNum = rigidComponent[i];
        if (atomNum < (int)atomTypes.size() && 
            atomTypes[atomNum].find("H") != 0) {
            heavyRigidAtomNums.push_back(atomNum);
            heavyRigidCount++;
        }
    }
}

// ============================================================================

/*
 * getRigidStructures — port of Python hierarchy._getRigidStructures
 *
 * Identifies sets of atoms that are "rigidly connected" (i.e. they always
 * move as a unit).  An atom is part of a rigid structure if it is:
 *   a. In a ring system (cycle), or
 *   b. Connected via a non-single bond (double/triple/aromatic/am), or
 *   c. A terminal atom (exactly one bond neighbour — that bond is treated
 *      as rigid regardless of type).
 *
 * Atoms not part of any such grouping get their own unique structure id.
 * Results stored in this->rigidStructures (indexed by atom) and
 * this->rigidStructureCount.
 */
void Hierarchy::getRigidStructures(
        int natoms,
        const std::vector<std::vector<std::pair<int, std::string>>>& atomBonds,
        bool verbose) {

    int cycleCount = 0;
    std::vector<std::vector<int>> cycles(natoms);
    std::vector<int> parent(natoms, 0);
    std::vector<int> visited(natoms, 0); /* 0=unseen, 1=in-progress, 2=done */

    /* DFS cycle finder (iterative to avoid stack overflow on large mols) */
    std::function<void(int, int)> findCycles = [&](int atom, int prev) {
        if (visited[atom] == 2) return;
        if (visited[atom] == 1) {
            /* Back-edge: trace path back to atom to mark the cycle */
            int curr = prev;
            while (curr != atom) {
                cycles[curr].push_back(cycleCount);
                curr = parent[curr];
            }
            cycles[curr].push_back(cycleCount);
            cycleCount++;
            return;
        }
        visited[atom] = 1;
        parent[atom] = prev;
        for (const auto& [other, btype] : atomBonds[atom]) {
            if (other == prev) continue;
            findCycles(other, atom);
        }
        visited[atom] = 2;
    };

    if (natoms > 0) findCycles(0, -1);

    /* Add "virtual cycles" for rigid bonds and terminal atoms */
    for (int atom = 0; atom < natoms; atom++) {
        if ((int)atomBonds[atom].size() == 1) {
            /* Terminal atom: bond to its sole neighbour is rigid */
            cycles[atom].push_back(cycleCount);
            cycles[atomBonds[atom][0].first].push_back(cycleCount);
            cycleCount++;
            continue;
        }
        for (const auto& [other, btype] : atomBonds[atom]) {
            if (other < atom) continue; /* process each bond once */
            /* Check if already in a common cycle */
            bool skip = false;
            for (int c : cycles[atom]) {
                for (int c2 : cycles[other]) {
                    if (c == c2) { skip = true; break; }
                }
                if (skip) break;
            }
            /* Non-single bonds are rigid */
            if (!skip && btype != "1") {
                cycles[atom].push_back(cycleCount);
                cycles[other].push_back(cycleCount);
                cycleCount++;
            }
        }
    }

    /* Find intersections between cycles */
    std::vector<std::set<int>> intersections(cycleCount);
    for (const auto& atomCycles : cycles) {
        int lc = atomCycles.size();
        for (int i = 0; i < lc; i++) {
            for (int j = i + 1; j < lc; j++) {
                intersections[atomCycles[i]].insert(atomCycles[j]);
                intersections[atomCycles[j]].insert(atomCycles[i]);
            }
        }
    }

    /* DFS over cycle-intersection graph to assign rigid structure ids */
    std::vector<int> rigidMap(cycleCount, 0);
    std::vector<bool> cycleVisited(cycleCount, false);
    int rsc = 0;

    std::function<void(int, int)> findRigidStructures = [&](int cycle, int prev) {
        if (cycleVisited[cycle]) return;
        cycleVisited[cycle] = true;
        rigidMap[cycle] = rsc;
        for (int intersect : intersections[cycle]) {
            if (intersect == prev) continue;
            findRigidStructures(intersect, cycle);
        }
    };

    for (int c = 0; c < cycleCount; c++) {
        if (!cycleVisited[c]) {
            findRigidStructures(c, -1);
            rsc++;
        }
    }

    /* Map atoms → rigid structure id */
    rigidStructures.assign(natoms, 0);
    for (int atom = 0; atom < natoms; atom++) {
        if (!cycles[atom].empty()) {
            rigidStructures[atom] = rigidMap[cycles[atom][0]];
        } else {
            /* Isolated atom (two single bonds, not in any cycle/rigid bond):
             * gets its own unique structure */
            rigidStructures[atom] = rsc++;
        }
    }
    rigidStructureCount = rsc;

    if (verbose) {
        std::cout << "Rigid structure count: " << rigidStructureCount << std::endl;
    }
}

/*
 * countPositionsPython — port of Python hierarchy._countPositions (new algorithm)
 *
 * Uses Buckets2 to cluster per-atom conformation positions, then resolves
 * those clusters into output conformations split along rigid-structure
 * boundaries.  Populates the same 1-indexed data structures used by the
 * writer methods so no writer changes are required.
 */
void Hierarchy::countPositionsPython(
        const std::vector<std::vector<std::array<double, 3>>>& xyzData,
        double tolerance, bool verbose) {

    if (xyzData.empty()) return;
    int nmol2s = (int)xyzData.size();
    int natoms = (int)xyzData[0].size();

    /* ---- Phase 1: cluster per-atom positions with Buckets2 ---- */
    ConfClusterMap confClusters;
    posCount.clear();

    Buckets2 bucketer(tolerance);
    int confNumCtr = 0;
    for (int atom = 0; atom < natoms; atom++) {
        std::vector<std::array<double, 3>> atomXyzs(nmol2s);
        for (int c = 0; c < nmol2s; c++) atomXyzs[c] = xyzData[c][atom];

        auto [npos_local, newConfNum] = bucketer.bucket(atomXyzs, atom, confNumCtr, confClusters);
        posCount.push_back(npos_local); /* number of distinct positions for this atom */
        confNumCtr = newConfNum;
    }

    /* ---- Phase 2: resolve clusters into confs split by rigidStructures ---- */

    /* Sort so the cluster containing ALL confs (the rigid cluster) is first.
     * Python sorts by descending key length: key = sorted conf indices. */
    std::vector<std::pair<std::vector<int>, ConfClusterValue*>> sortedClusters;
    sortedClusters.reserve(confClusters.size());
    for (auto& kv : confClusters)
        sortedClusters.push_back({kv.first, &kv.second});
    std::sort(sortedClusters.begin(), sortedClusters.end(),
              [](const auto& a, const auto& b) { return a.first.size() > b.first.size(); });

    /* Reset output data structures */
    outAtoms = 0;
    outAtomOrigAtom.clear();
    outAtomConfNum.clear();
    outAtomXYZMap.clear();
    confNumAtomList.clear();
    confNums.clear();
    setToConfs.clear();
    for (int s = 0; s < nmol2s; s++) setToConfs[s] = {};

    int globalAtomCnt = 0; /* 0-based running coord counter */
    int confNum_act   = 0; /* 0-based, incremented before first use → 1-based when stored */

    for (auto& [tupleInput, clusterPtr] : sortedClusters) {
        const std::vector<int>& atoms   = clusterPtr->atomIds;
        const std::vector<std::array<double,3>>& xyzlist = clusterPtr->xyzPositions;
        int n = (int)atoms.size();

        /* Sort atoms by rigid structure id (mirroring Python rsatoms sort) */
        std::vector<std::pair<int,int>> rsatoms(n); /* (atomId, rigidStructureId) */
        for (int i = 0; i < n; i++)
            rsatoms[i] = {atoms[i], rigidStructures[atoms[i]]};
        std::sort(rsatoms.begin(), rsatoms.end(),
                  [](const auto& a, const auto& b){ return a.second < b.second; });

        int rs_prev = -1;
        for (int i = 0; i < n; i++) {
            int rs = rsatoms[i].second;

            if (rs == rs_prev) {
                /* Same rigid structure: add to current conf */
                confNumAtomList[confNum_act].push_back(globalAtomCnt + 1);
            } else {
                /* New rigid structure: start a new conf */
                confNum_act++;
                confNums.push_back(confNum_act);
                confNumAtomList[confNum_act] = {globalAtomCnt + 1};

                for (int setno : tupleInput)
                    setToConfs[setno].push_back(confNum_act);
            }

            /* Store coord entry (1-indexed) — uses original-order index i
             * to match Python's atoms[i] / xyzlist[i] behaviour */
            int coordIdx = globalAtomCnt + 1; /* 1-based */
            outAtomOrigAtom[coordIdx]  = atoms[i];      /* original order (Python faithful) */
            outAtomConfNum[coordIdx]   = confNum_act;
            outAtomXYZMap[coordIdx]    = xyzlist[i];    /* representative xyz */

            globalAtomCnt++;
            rs_prev = rs;
        }
    }

    outAtoms = globalAtomCnt;
    numConfs = confNum_act;

    if (verbose) {
        std::cout << "Total output coords: " << outAtoms << std::endl;
        std::cout << "Total output confs:  " << numConfs << std::endl;
    }
}

/*
 * identifyClashSetnumsPython — fast clash detection matching Python's
 * decideDistanceRules() approach.  Checks every input conformation.
 */
void Hierarchy::identifyClashSetnumsPython(const Clash& clashDecider,
                                           const Mol2& mol2data) {
    brokenSets.clear();
    int nconfs = (int)mol2data.atomXyz.size();
    for (int aSet = 0; aSet < nconfs; aSet++) {
        if (clashDecider.decideDistanceRules(mol2data, mol2data.atomXyz[aSet]))
            brokenSets.push_back(aSet);
    }
}

// ============================================================================

void Hierarchy::write(const std::string& db2gzFileName,
                     bool verbose, bool timeit, long limitset,
                     const std::string& writeMode) {

    /* When output is "-", write uncompressed plain text to stdout */
    if (db2gzFileName == "-") {
        colorWriter(std::cout, *mol2data);
        allButSetWriter(std::cout, *mol2data, setToConfs.size());
        setWriter(std::cout, *mol2data);
        std::cout << "E\n";
        std::cout.flush();
        return;
    }

    gzFile outFile = gzopen(db2gzFileName.c_str(),
                           (writeMode == "a") ? "ab" : "wb");
    if (outFile == NULL) {
        std::cerr << "Error opening output file: " << db2gzFileName << std::endl;
        return;
    }

    std::ostringstream buffer;

    /* Write color table if changed */
    colorWriter(buffer, *mol2data);

    /* Write main data */
    allButSetWriter(buffer, *mol2data, setToConfs.size());
    setWriter(buffer, *mol2data);

    buffer << "E\n";

    std::string output = buffer.str();
    gzwrite(outFile, output.c_str(), output.length());
    gzclose(outFile);

    if (verbose) {
        std::cout << db2gzFileName << " written" << std::endl;
    }
}

void Hierarchy::colorWriter(std::ostream& outFile, const Mol2& mol2data) {
    if (!mol2data.colorConverter) return;
    auto defaultInts = ColorTable::getDefaultColorInts();
    if (mol2data.colorConverter->colorInts == defaultInts) return;

    /* Non-default color table — write T lines sorted by color integer */
    std::vector<std::pair<int, std::string>> colors;
    for (const auto& [name, key] : mol2data.colorConverter->colorInts)
        colors.push_back({key, name});
    std::sort(colors.begin(), colors.end());
    for (const auto& [key, name] : colors)
        outFile << "T " << std::setw(2) << key << " " << std::setw(8) << name << "\n";
}

void Hierarchy::allButSetWriter(std::ostream& outFile, const Mol2& mol2data,
                                int setsTotal) {

    /* Molecule header */
    outFile << "M " << std::setw(16) << mol2data.name.substr(0, 16)
            << " " << std::setw(9) << mol2data.protName.substr(0, 9)
            << " " << std::setw(3) << mol2data.atomNum.size()
            << " " << std::setw(3) << mol2data.bondNum.size()
            << " " << std::setw(6) << outAtoms
            << " " << std::setw(6) << numConfs
            << " " << std::setw(6) << setsTotal
            << " " << std::setw(6) << heavyRigidCount
            << " " << std::setw(6) << 5
            << " " << std::setw(6) << 0 << "\n";
    
    /* Solvation data */
    outFile << "M " << std::setw(9) << std::setprecision(4) << std::fixed
           << std::showpos << mol2data.solvData.totalCharge
           << " " << std::setw(10) << std::setprecision(3)
           << mol2data.solvData.totalPolarSolv
           << " " << std::setw(10) << mol2data.solvData.totalApolarSolv
           << " " << std::setw(10) << mol2data.solvData.totalSolv
           << " " << std::setw(9) << std::noshowpos
           << mol2data.solvData.totalSurface << "\n";
    
    /* SMILES and name */
    outFile << "M " << mol2data.smiles.substr(0, 76) << "\n";
    outFile << "M " << mol2data.longname.substr(0, 76) << "\n";
    outFile << "M " << std::setw(10) << std::setprecision(4) << std::showpos
           << 999.999 << std::noshowpos << "\n";
    
    /* Atoms */
    for (size_t i = 0; i < mol2data.atomNum.size(); i++) {
        outFile << "A " << std::setw(3) << mol2data.atomNum[i]
               << " " << std::setw(4) << std::left << mol2data.atomName[i]
               << " " << std::setw(5) << mol2data.atomType[i] << std::right
               << " " << std::setw(2) << mol2data.dockNum[i]
               << " " << std::setw(2) << mol2data.colorNum[i]
               << " " << std::setw(9) << std::setprecision(4) << std::showpos 
               << mol2data.solvData.charge[i]
               << " " << std::setw(10) << std::setprecision(3)
               << mol2data.solvData.polarSolv[i]
               << " " << std::setw(10) << mol2data.solvData.apolarSolv[i]
               << " " << std::setw(10) << mol2data.solvData.solv[i]
               << " " << std::setw(9) << std::noshowpos
               << mol2data.solvData.surface[i] << "\n";
    }
    
    /* Bonds */
    for (size_t i = 0; i < mol2data.bondNum.size(); i++) {
        outFile << "B " << std::setw(3) << mol2data.bondNum[i]
               << " " << std::setw(3) << mol2data.bondStart[i]
               << " " << std::setw(3) << mol2data.bondEnd[i]
               << " " << std::setw(2) << std::left 
               << mol2data.bondType[i] << std::right << "\n";
    }
    
    /* Coordinates */
    for (int xyzNum = 1; xyzNum <= outAtoms; xyzNum++) {
        int atomNum = outAtomOrigAtom[xyzNum];
        int confNum = outAtomConfNum[xyzNum];
        const std::array<double, 3>& xyz = outAtomXYZMap[xyzNum];
        outFile << "X " << std::setw(9) << xyzNum
                << " " << std::setw(3) << (atomNum + 1)
                << " " << std::setw(6) << confNum
                << " " << std::setw(9) << std::setprecision(4) << std::showpos
                << xyz[0]
                << " " << std::setw(9) << xyz[1]
                << " " << std::setw(9) << xyz[2] << std::noshowpos << "\n";
    }

    /* Rigid matching spheres */
    rigidNumSeen = 0;
    for (int rigidNum : heavyRigidAtomNums) {
        rigidNumSeen++;
        int atomColor = mol2data.colorNum[rigidNum];
        std::array<double, 3> xyz = mol2data.atomXyz[0][rigidNum];
        
        outFile << "R " << std::setw(6) << rigidNumSeen
               << " " << std::setw(2) << atomColor
               << " " << std::setw(9) << std::setprecision(4) << std::showpos
               << xyz[0]
               << " " << std::setw(9) << xyz[1]
               << " " << std::setw(9) << xyz[2] << std::noshowpos << "\n";
    }
    
    /* Conformations */
    for (int confNum : confNums) {
        const std::vector<int>& atomList = confNumAtomList[confNum];
        int coordStart = *std::min_element(atomList.begin(), atomList.end());
        int coordEnd   = *std::max_element(atomList.begin(), atomList.end());
        outFile << "C " << std::setw(6) << confNum
                << " " << std::setw(9) << coordStart
                << " " << std::setw(9) << coordEnd << "\n";
    }
}

void Hierarchy::setWriter(std::ostream& outFile, const Mol2& mol2data) {
    std::vector<int> curSets;
    for (auto& pair : setToConfs) {
        curSets.push_back(pair.first);
    }
    std::sort(curSets.begin(), curSets.end());
    
    for (size_t i = 0; i < curSets.size(); i++) {
        int curSet = curSets[i];
        int outSetNum = (i + 1);
        
        const std::vector<int>& curConfs = setToConfs[curSet];
        int totalConfs = curConfs.size();
        
        int totalLines = (int)std::ceil(totalConfs / (double)coSePerLine);
        int lastLineLen = totalConfs % (int)coSePerLine;
        if (lastLineLen == 0) lastLineLen = coSePerLine;
        
        int brokenSet = (std::find(brokenSets.begin(), brokenSets.end(), curSet)
                         != brokenSets.end()) ? 1 : 0;

        /* Always write two strain values (totalStrain + maxStrain), matching
         * the Python strain format.  Fall back to 9999.99 if not present. */
        double totalStrain = (curSet < (int)mol2data.inputTotalStrain.size()) ?
            mol2data.inputTotalStrain[curSet] : 9999.999;
        double maxStrain = (curSet < (int)mol2data.inputMaxStrain.size()) ?
            mol2data.inputMaxStrain[curSet] : 9999.999;

        int outHydro = (curSet < (int)mol2data.inputHydrogens.size()) ?
            mol2data.inputHydrogens[curSet] : 3;

        outFile << "S " << std::setw(6) << outSetNum
                << " " << std::setw(6) << totalLines
                << " " << std::setw(3) << totalConfs
                << " " << std::setw(1) << brokenSet
                << " " << std::setw(1) << outHydro
                << " " << std::setw(11) << std::setprecision(3) << std::showpos
                << totalStrain
                << " " << std::setw(11) << maxStrain << std::noshowpos << "\n";
        
        /* Write conformation lists */
        for (int lineNum = 0; lineNum < totalLines; lineNum++) {
            int numInLine = (lineNum == totalLines - 1) ? lastLineLen : coSePerLine;
            
            outFile << "S " << std::setw(6) << outSetNum
                   << " " << std::setw(6) << (lineNum + 1)
                   << " " << std::setw(1) << numInLine;
            
            for (int i = 0; i < numInLine; i++) {
                int confIdx = lineNum * coSePerLine + i;
                outFile << " " << std::setw(6) << curConfs[confIdx];
            }
            outFile << "\n";
        }
    }
}


void Hierarchy::writeMol2(const std::string& mol2FileName,
                          bool verbose, bool timeit, bool separateClusters) {
    /* Simplified - just write all conformations to one file */
    std::ofstream outFile(mol2FileName);
    if (outFile.is_open()) {
        mol2data->writeMol2File(outFile, nullptr);
        outFile.close();
        if (verbose) {
            std::cout << mol2FileName << " written" << std::endl;
        }
    }
}

} // namespace mol2db2
