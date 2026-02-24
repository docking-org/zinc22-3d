#ifndef CLASH_H
#define CLASH_H

#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <variant>

namespace mol2db2 {

/* Forward declaration */
class Mol2;

/*
 * Clash detection class.
 * Holds parameters that determine what a clashed conformation contains.
 * 
 * Rules have the format:
 * ("min|max", bond, cmp, "X", "Y", dist, dist_squared)
 * 
 * min|max - whether the constraint is a minimum or maximum
 * bond - how many bonds are between the atoms
 * cmp - how to use the bond distance:
 *       -1 means there must be less than that many bonds
 *        0 means there must be exactly that many bonds
 *        1 means there must be more than that many bonds
 * X - atom type X, "*" means all
 * Y - atom type Y, "*" means all
 * dist - float that is the distance constraint
 * dist_squared - precomputed square of dist for speed
 */
class Clash {
public:
    Clash(const std::string* parameterFileName = nullptr);
    
    void printParameters() const;
    
    /*
     * Decide if a conformation is clashed (new fast algorithm).
     * Only computes distances between same-type atom pairs, matching Python's
     * decideDistanceRules behaviour.
     * Returns true if clashed, false if not clashed.
     */
    bool decideDistanceRules(const Mol2& mol2data,
                             const std::vector<std::array<double, 3>>& xyzData) const;

private:
    /* Rule format: (min/max, bond, cmp, atomType1, atomType2, dist, dist^2) */
    using Rule = std::tuple<std::string, int, int, std::string, std::string, 
                            double, double>;
    
    std::vector<Rule> rules;
    
    /* Default rules: only H-H with minimum distance 1.70 at 2+ bonds away */
    static std::vector<Rule> getDefaultRules();
};

} // namespace mol2db2

#endif // CLASH_H
