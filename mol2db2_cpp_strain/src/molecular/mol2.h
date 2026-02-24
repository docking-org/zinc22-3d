#ifndef MOL2_H
#define MOL2_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <array>
#include <tuple>
#include <memory>
#include <variant>

namespace mol2db2 {

/* Forward declarations */
class UnionFind;

// ============================================================================
// AtomConverter (formerly sybyl2dock.py)
// ============================================================================
class AtomConverter {
public:
    AtomConverter(const std::string* parameterFileName = nullptr);
    void printParameters() const;
    int convertMol2atomNum(const class Mol2& mol2data, int atomNum) const;

private:
    std::unordered_map<std::string, int> convertTypes;
    std::unordered_map<std::string, std::string> specialKeys;
    
    static std::unordered_map<std::string, int> getDefaultConvertTypes();
};

// ============================================================================
// ColorTable (formerly atom_color_table.py)
// ============================================================================
class ColorTable {
public:
    ColorTable(const std::string* parameterFileName = nullptr);
    void printParameters() const;
    int convertMol2color(const class Mol2& mol2data, int atomNum) const;

    std::string defaultColor;
    std::unordered_map<std::string, int> colorInts;

    static std::unordered_map<std::string, int> getDefaultColorInts();

private:
    /* Rule2: (atom_type, color) */
    /* Rule4: (atom_type, bonds_away, other_atom_type, color) */
    using Rule2 = std::tuple<std::string, std::string>;
    using Rule4 = std::tuple<std::string, int, std::string, std::string>;
    
    std::vector<std::variant<Rule2, Rule4>> rulesTable;
    
    static std::string getDefaultColorDefault();
    static std::vector<std::variant<Rule2, Rule4>> getDefaultRulesTable();
};

// ============================================================================
// SolvData (formerly solv.py)
// ============================================================================
struct SolvData {
    std::string name;
    int totalAtoms;
    double totalCharge;
    double totalPolarSolv;
    double totalSurface;
    double totalApolarSolv;
    double totalSolv;
    
    std::vector<double> charge;
    std::vector<double> polarSolv;
    std::vector<double> apolarSolv;
    std::vector<double> solv;
    std::vector<double> surface;
    
    SolvData();
    void readSolvFile(const std::string& filename);
};

// ============================================================================
// Mol2 (integrates mol2.py + solv.py + sybyl2dock.py + atom_color_table.py)
// ============================================================================
class Mol2 {
public:
    /* Constructors */
    Mol2();
    Mol2(const std::string* mol2FileName, 
         const std::string* nameFileName = nullptr,
         const std::vector<std::string>* mol2text = nullptr);
    ~Mol2();

    // this is the copy constructor
    Mol2(const Mol2&); 

    /* Copy operations */
    Mol2 copy() const; // move outside of class?
    
    /* File I/O */
    void processLine(const std::string& line);
    void writeMol2File(std::ostream& outFile, 
                       const std::vector<int>* whichXyz = nullptr) const;
    void writeMol2(const std::string& outName, 
                   const std::vector<int>* whichXyz = nullptr) const;
    
    /* Bond queries */
    int bondsBetween(int atomNum, int atomOther) const;
    int bondsBetweenActual(int actualNum, int actualOther) const;
    bool bondedTo(int atomNum, const std::string& firstName, 
                  int bondsAway = 1,
                  const std::string* lastBond = nullptr, 
                  bool returnAtom = false,
                  int* returnedAtomNum = nullptr) const;
    std::map<int, std::vector<int>> bondedToAll(int atomNum, 
                                                 const std::string& firstName,
                                                 int bondsAway = 1,
                                                 const std::string* lastBond = nullptr) const;
    bool isAtomBondedOtherThan(int atomNum, const std::set<int>& count,
                               const std::set<std::string>& otherThan) const;
    
    /* Distance calculations */
    void calcBondDists();
    std::vector<int> distFromAtoms(const std::vector<int>& atoms) const;
    
    /* Type conversions */
    void convertDockTypes(const std::string* parameterFileName = nullptr);
    void addColors(const std::string* parameterFileName = nullptr);
    
    /* Conformation operations */
    int countConfs() const { return atomXyz.size(); }
    void keepConfsOnly(int first, int last);
    std::array<double, 3> getXyz(int xyzCount, int atomNum) const;
    std::vector<std::array<double, 3>> getXyzManyConfs(
        const std::vector<int>& xyzCounts, int atomIndex) const;
    
    /* RMSD and clustering */
    double getRMSD(int xyzOne, int xyzTwo) const;
    std::map<int, std::map<int, double>> getRMSDtable(bool forceRedo = false);
    std::vector<std::tuple<double, int, int>> getRMSDlist();
    std::vector<std::vector<int>> getRMSDclusters(
        const double* rmsdCutoff = nullptr, int numClusters = 1);
    std::vector<std::vector<int>> getRMSDclustersAll(
        const double* rmsdCutoff = nullptr, int numClusters = 1);
    std::vector<std::vector<int>> divisiveClustering();
    
    /* Solvation data */
    void addSolvDataPartialCharges(const std::vector<double>& partialCharges);
    
    /* Atom/bond deletion (for covalent ligands) */
    void deleteBond(int bondInd);
    void deleteAtom(int atomInd);
    std::tuple<std::string, std::vector<int>> removeCovalentDummyAtom();
    void recolorCovalentAttachment(const std::string& covAtomType);
    
    /* Public data members */
    std::string name;
    std::string protName;
    std::string smiles;
    std::string longname;
    
    std::vector<int> atomNum;
    std::vector<std::string> atomName;
    std::vector<std::string> atomType;
    std::vector<double> atomCharge;
    std::vector<std::vector<std::pair<int, std::string>>> atomBonds;
    
    std::vector<int> bondNum;
    std::vector<int> bondStart;
    std::vector<int> bondEnd;
    std::vector<std::string> bondType;
    
    std::vector<std::vector<std::array<double, 3>>> atomXyz;
    std::vector<double> inputEnergy;
    std::vector<double> inputTotalStrain;
    std::vector<double> inputMaxStrain;
    std::vector<int> inputHydrogens;
    
    int xyzCount;
    int origXyzCount;
    
    /* Dock atom types */
    std::vector<int> dockNum;
    
    /* Color data */
    std::vector<int> colorNum;
    std::unique_ptr<ColorTable> colorConverter;
    
    /* Solvation data */
    SolvData solvData;
    
    /* Hydrogen rotation data */
    std::vector<std::string> hydrogenRotAngles;
    std::map<int, std::vector<int>> dihedrals;
    std::map<int, std::vector<double>> rotAngles;
    int hydrogensToRotate;
    
private:
    void blankNew();
    std::map<int, std::vector<int>> bondedToActualAll(
        int actualNum, const std::string& firstName,
        int bondsAway = 1, const std::string* lastBond = nullptr) const;
    bool bondedToActual(int actualNum, const std::string& firstName,
                        int bondsAway = 1, const std::string* lastBond = nullptr,
                        bool returnAtom = false, 
                        int* returnedAtomNum = nullptr) const;
    
    /* Bond distance caching */
    mutable std::vector<std::vector<int>> bondDists;
    mutable std::unordered_map<int, int> bondDistsOrderKeys;
    
    /* RMSD caching */
    mutable std::map<int, std::map<int, double>> rmsdTable;
    mutable std::vector<std::tuple<double, int, int>> rmsdList;
    
    /* File reading state */
    int phase;
};

/* Utility function (from readDockMol2file in mol2.py) */
std::tuple<std::vector<Mol2>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>>
readDockMol2file(const std::string& mol2FileName, bool recdes = false,
                 bool ligdes = false, bool charge = false, bool elec = false);

} // namespace mol2db2

#endif // MOL2_H
