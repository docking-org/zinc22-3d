#ifndef HIERARCHY_H
#define HIERARCHY_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <array>

namespace mol2db2 {

/* Forward declarations */
class Mol2;
class Clash;

/*
 * Exception thrown when hierarchy becomes too large
 */
class TooBigError {
public:
    TooBigError(int confs, int sets, int coords)
        : confs(confs), sets(sets), coords(coords) {}
    
    int getConfs() const { return confs; }
    int getCoords() const { return coords; }
    int getSets() const { return sets; }
    
private:
    int confs;
    int coords;
    int sets;
};

/*
 * Hierarchy class - builds conformational hierarchy for ligands
 * Simplified version focusing on core functionality
 */
class Hierarchy {
public:
    /* Constants for output formatting */
    static const int grGrPerLine = 17;
    static const int grCoPerLine = 9;
    static const int coCoPerLine = 9;
    static const int coSePerLine = 8;
    
    /*
     * Constructor - builds hierarchy from mol2 data using the
     * _getRigidStructures + _countPositions algorithm (mol2db2_py3_strain).
     */
    Hierarchy(Mol2& mol2data,
              const Clash& clashDecider,
              double tolerance = 0.001,
              bool verbose = false,
              bool timeit = false,
              long limitset = 9999999999,
              long limitconf = 9999999999,
              long limitcoord = 9999999999);
    
    /*
     * Write hierarchy to db2.gz file
     */
    void write(const std::string& db2gzFileName,
               bool verbose = false,
               bool timeit = false,
               long limitset = 9999999,
               const std::string& writeMode = "w");
    
    /*
     * Write as multi-mol2 file (for debugging)
     */
    void writeMol2(const std::string& mol2fileName,
                   bool verbose = false,
                   bool timeit = false,
                   bool separateClusters = true);
    
private:
    /* Position counting */
    std::vector<int> posCount;

    /* Rigid component */
    std::vector<int> rigidComponent;
    std::set<int> atomsAssigned;
    std::set<int> atomsNotAssigned;
    int heavyRigidCount;
    std::vector<int> heavyRigidAtomNums;

    /* Conformations */
    std::vector<int> confNums;

    /* Sets */
    std::map<int, std::vector<int>> setToConfs;

    /* Output atoms */
    int outAtoms;
    std::map<int, int> outAtomOrigAtom;
    std::map<int, int> outAtomConfNum;
    std::map<int, std::vector<int>> confNumAtomList;

    /* Clashing */
    std::vector<int> brokenSets;

    /* Reference to mol2 data */
    Mol2* mol2data;

    /* Per-atom xyz stored directly (1-indexed: key 1..outAtoms) */
    std::map<int, std::array<double, 3>> outAtomXYZMap;

    /* Rigid structure id per atom (atom index → structure id) */
    std::vector<int> rigidStructures;
    int rigidStructureCount;

    /* Total number of output conformations */
    int numConfs;

    void getRigidStructures(int natoms,
                            const std::vector<std::vector<std::pair<int, std::string>>>& atomBonds,
                            bool verbose);
    void countPositionsPython(const std::vector<std::vector<std::array<double, 3>>>& xyzData,
                              double tolerance, bool verbose);
    void findRigidComponent(const std::vector<std::vector<std::pair<int, std::string>>>& atomBonds);
    void findRigidHeavy(const std::vector<std::string>& atomTypes);
    void identifyClashSetnumsPython(const Clash& clashDecider, const Mol2& mol2data);

    /* Writer helpers */
    void colorWriter(std::ostream& outFile, const Mol2& mol2data);
    void allButSetWriter(std::ostream& outFile, const Mol2& mol2data, int setsTotal);
    void setWriter(std::ostream& outFile, const Mol2& mol2data);

    int rigidNumSeen;
};

/*
 * Helper function to compute number of breaks needed
 */
int computeBreaks(const TooBigError& limitError, 
                  long limitset, long limitconf, long limitcoord);

} // namespace mol2db2

#endif // HIERARCHY_H
