#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <vector>

namespace mol2db2 {

/* 
 * Takes a list of lists. Returns all possible ways of picking one element
 * from each list, as a new list.
 * Example: [[1,2], [3,4]] -> [[1,3], [1,4], [2,3], [2,4]]
 */
template<typename T>
std::vector<std::vector<T>> allCombinations(const std::vector<std::vector<T>>& inputLists);

} // namespace mol2db2

/* Template implementation must be in header */
namespace mol2db2 {

template<typename T>
std::vector<std::vector<T>> allCombinations(const std::vector<std::vector<T>>& inputLists) {
    std::vector<std::vector<T>> stack;
    stack.push_back(std::vector<T>());  /* Start with empty list */
    
    for (size_t i = 0; i < inputLists.size(); i++) {
        std::vector<std::vector<T>> newStack;
        
        for (size_t j = 0; j < stack.size(); j++) {
            for (size_t k = 0; k < inputLists[i].size(); k++) {
                std::vector<T> oldCopy = stack[j];
                oldCopy.push_back(inputLists[i][k]);
                newStack.push_back(oldCopy);
            }
        }
        
        stack = newStack;
    }
    
    return stack;
}

} // namespace mol2db2

#endif // COMBINATORICS_H
