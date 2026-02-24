#ifndef UNIONFIND_H
#define UNIONFIND_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>

namespace mol2db2 {

/*
 * Union-Find data structure with path compression and union by rank.
 * Based on CLRS algorithms textbook.
 */
class UnionFind {
public:
    UnionFind();
    
    std::unordered_map<int, int> ranks;
    std::unordered_map<int, int> parents;

    /* Print internal state (for debugging) */
    void printPar() const;
    
    /* Check if name has been seen before */
    bool check(int name) const;
    
    /* Find root of name's set, add if not present */
    int find(int name);
    
    /* Union two sets */
    int unionSets(int name, int other);
    
    /* Check if two items are in different sets */
    bool different(int itemA, int itemB);
    
    /* Get all items in same set as name (O(n) operation) */
    std::vector<int> getList(int name);
    
    /* Convert to list of lists, one per set (O(n) operation) */
    std::vector<std::vector<int>> toLists();
    
//private:
//    std::unordered_map<int, int> ranks;
//    std::unordered_map<int, int> parents;
};

/*
 * UnionFind with attached data carried with each union.
 * Unions data together using set union.
 */
class UnionFindAttach : public UnionFind {
public:
    UnionFindAttach();
    
    /* Find with optional data attachment */
    int find(int name, const std::set<int>* attachData = nullptr);
    
    /* Union two sets and merge their attached data */
    int unionSets(int name, int other);
    
    /* Get attached set for name */
    std::set<int> getAttached(int name);
    
    /* Clear attached set for name */
    void clearAttached(int name);
    
private:
    std::unordered_map<int, std::set<int>> attached;
};

} // namespace mol2db2

#endif // UNIONFIND_H
