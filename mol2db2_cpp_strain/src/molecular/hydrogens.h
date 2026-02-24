#ifndef HYDROGENS_H
#define HYDROGENS_H

#include <string>
#include <vector>
#include <tuple>

namespace mol2db2 {

/* Forward declaration */
class Mol2;

/*
 * Handles terminal rotatable hydrogen definitions.
 * Rules have the format:
 * (type, Atom, bond, Atom, bond, Atom, Degrees)
 * 
 * Type 1: (1, "C.ar", "1", "S", "1", "H", "180")
 *   Atom - Atom name like "C.ar" or "C" which matches all starting with C
 *   bond - "1", "2", "3" or "ar" or "*" which means any
 *   Degrees - "120,240" or "180" or "-" (no rotation)
 * 
 * Type 2: (2, bond, bond, Atom, bond, Atom, Degrees)
 *   where the two bonds at front are both applied to the same Atom (the first)
 * 
 * Rules are processed in order so you can exclude certain things with "-"
 * All things not specified are "-" (no rotation)
 */
class Hydrogens {
public:
    Hydrogens(const std::string* parameterFileName = nullptr);
    
    void printParameters() const;
    
    /* Find all terminal hydrogens that can be rotated */
    void findTerminalHydrogens(Mol2* mol2data);
    
    /* Reset planar hydrogens to 0 degrees */
    void resetHydrogens(Mol2* mol2data);
    
    /* Rotate hydrogens to all specified angles, returns new Mol2 */
    Mol2 rotateHydrogens(const Mol2& mol2data);
    void rotateHydrogens_teb(Mol2& mol2data); // teb added. modify the reference mol2data.  
    
    /* Find current angles of terminal hydrogens */
    void findAngles(const Mol2& mol2data);
    
private:
    /* Rule format: (type, atom1/bond1, bond1/bond2, atom2, bond2, atom3, degrees) */
    using Rule = std::tuple<int, std::string, std::string, std::string, 
                           std::string, std::string, std::string>;
    
    std::vector<Rule> rules;
    
    /* Find dihedrals for rotation */
    void findDihedrals(Mol2* mol2data);
    
    /* Get current dihedral angle */
    std::tuple<std::vector<std::array<double, 3>>, double> 
    getCurDihedral(int atomNum, int xyzCount, const Mol2& mol2data);
    
    /* Default rules */
    static std::vector<Rule> getDefaultRules();
};

} // namespace mol2db2

#endif // HYDROGENS_H
