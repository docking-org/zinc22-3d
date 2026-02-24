#include "unionfind.h"
#include <iostream>

namespace mol2db2 {

/* UnionFind implementation */

UnionFind::UnionFind() {
}

void UnionFind::printPar() const {
    std::cout << "Parents: ";
    for (auto it = parents.begin(); it != parents.end(); ++it) {
        std::cout << it->first << "->" << it->second << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Ranks: ";
    for (auto it = ranks.begin(); it != ranks.end(); ++it) {
        std::cout << it->first << ":" << it->second << " ";
    }
    std::cout << std::endl;
}

bool UnionFind::check(int name) const {
    return (parents.find(name) != parents.end());
}

int UnionFind::find(int name) {
    /* If not seen before, create singleton set */
    if (parents.find(name) == parents.end()) {
        parents[name] = name;
        ranks[name] = 0;
        return name;
    }
    
    /* Path compression: find root and compress path */
    std::vector<int> path;
    path.push_back(name);
    int parent = parents[name];
    
    while (parent != path.back()) {
        path.push_back(parent);
        parent = parents[parent];
    }
    
    /* Compress path */
    for (size_t i = 0; i < path.size() - 1; i++) {
        parents[path[i]] = parent;
    }
    
    return parent;
}

int UnionFind::unionSets(int name, int other) {
    int onePar = find(name);
    int otherPar = find(other);
    
    if (onePar == otherPar) {
        return onePar;  /* Already in same set */
    }
    
    /* Union by rank */
    if (ranks[onePar] < ranks[otherPar]) {
        parents[onePar] = otherPar;
        return otherPar;
    } else {
        parents[otherPar] = onePar;
        if (ranks[onePar] == ranks[otherPar]) {
            ranks[onePar]++;
        }
        return onePar;
    }
}

bool UnionFind::different(int itemA, int itemB) {
    int parA = find(itemA);
    int parB = find(itemB);
    return (parA != parB);
}

std::vector<int> UnionFind::getList(int name) {
    int parent = find(name);
    std::vector<int> returnList;
    
    /* Make all pointers direct */
    for (auto it = parents.begin(); it != parents.end(); ++it) {
        find(it->first);
    }
    
    /* Collect all items with same parent */
    for (auto it = parents.begin(); it != parents.end(); ++it) {
        if (it->second == parent) {
            returnList.push_back(it->first);
        }
    }
    
    return returnList;
}

std::vector<std::vector<int>> UnionFind::toLists() {
    std::unordered_map<int, std::vector<int>> lists;
    
    /* Make all pointers direct */
    for (auto it = parents.begin(); it != parents.end(); ++it) {
        find(it->first);
    }
    
    /* Group by parent */
    for (auto it = parents.begin(); it != parents.end(); ++it) {
        lists[it->second].push_back(it->first);
    }
    
    /* Convert to vector of vectors */
    std::vector<std::vector<int>> listOfLists;
    for (auto it = lists.begin(); it != lists.end(); ++it) {
        listOfLists.push_back(it->second);
    }
    
    return listOfLists;
}

/* UnionFindAttach implementation */

UnionFindAttach::UnionFindAttach() : UnionFind() {
}

int UnionFindAttach::find(int name, const std::set<int>* attachData) {
    int parent = UnionFind::find(name);
    
    /* Initialize attached set if needed */
    if (attached.find(parent) == attached.end()) {
        attached[parent] = std::set<int>();
    }
    
    /* Add attached data if provided */
    if (attachData != nullptr) {
        attached[parent].insert(attachData->begin(), attachData->end());
    }
    
    return parent;
}

int UnionFindAttach::unionSets(int name, int other) {
    int onePar = find(name, nullptr);
    int otherPar = find(other, nullptr);
    
    if (onePar == otherPar) {
        return onePar;
    }
    
    /* Perform union using parent class logic */
    int result;
    if (ranks[onePar] < ranks[otherPar]) {
        parents[onePar] = otherPar;
        /* Merge attached data */
        attached[otherPar].insert(attached[onePar].begin(), attached[onePar].end());
        result = otherPar;
    } else {
        parents[otherPar] = onePar;
        if (ranks[onePar] == ranks[otherPar]) {
            ranks[onePar]++;
        }
        /* Merge attached data */
        attached[onePar].insert(attached[otherPar].begin(), attached[otherPar].end());
        result = onePar;
    }
    
    return result;
}

std::set<int> UnionFindAttach::getAttached(int name) {
    int parent = find(name, nullptr);
    return attached[parent];
}

void UnionFindAttach::clearAttached(int name) {
    int parent = find(name, nullptr);
    attached[parent].clear();
}

} // namespace mol2db2
