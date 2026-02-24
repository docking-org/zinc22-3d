#ifndef PRIODICT_H
#define PRIODICT_H

#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

namespace mol2db2 {

/*
 * Priority dictionary using binary heaps.
 * Based on David Eppstein's Python implementation.
 */
template<typename K, typename V>
class PriorityDictionary {
public:
    PriorityDictionary();
    
    /* Find smallest item after removing deleted items from heap */
    K smallest();
    
    /* Set/update value for key */
    void set(const K& key, const V& val);
    
    /* Get value for key */
    V get(const K& key) const;
    
    /* Check if key exists */
    bool contains(const K& key) const;
    
    /* Remove key */
    void remove(const K& key);
    
    /* Get size */
    size_t size() const { return dict.size(); }
    
    /* Check if empty */
    bool empty() const { return dict.empty(); }
    
    /* Pop smallest element */
    K pop();
    
private:
    struct HeapPair {
        V value;
        K key;
        
        bool operator<(const HeapPair& other) const {
            return value < other.value;
        }
        
        bool operator<=(const HeapPair& other) const {
            return value <= other.value;
        }
        
        bool operator>(const HeapPair& other) const {
            return value > other.value;
        }
    };
    
    std::vector<HeapPair> heap;
    std::unordered_map<K, V> dict;
    
    void rebuildHeap();
};

} // namespace mol2db2

/* Template implementation */
namespace mol2db2 {

template<typename K, typename V>
PriorityDictionary<K, V>::PriorityDictionary() {
}

template<typename K, typename V>
K PriorityDictionary<K, V>::smallest() {
    if (dict.empty()) {
        throw std::out_of_range("smallest of empty PriorityDictionary");
    }
    
    /* Remove deleted items from top of heap */
    while (!heap.empty()) {
        HeapPair& top = heap[0];
        
        /* Check if this heap entry is still valid */
        auto it = dict.find(top.key);
        if (it != dict.end() && it->second == top.value) {
            return top.key;
        }
        
        /* Remove invalid entry and re-heapify */
        HeapPair lastItem = heap.back();
        heap.pop_back();
        
        if (heap.empty()) {
            break;
        }
        
        size_t insertionPoint = 0;
        while (true) {
            size_t smallChild = 2 * insertionPoint + 1;
            
            if (smallChild + 1 < heap.size() && 
                heap[smallChild] > heap[smallChild + 1]) {
                smallChild++;
            }
            
            if (smallChild >= heap.size() || lastItem <= heap[smallChild]) {
                heap[insertionPoint] = lastItem;
                break;
            }
            
            heap[insertionPoint] = heap[smallChild];
            insertionPoint = smallChild;
        }
    }
    
    if (heap.empty() && !dict.empty()) {
        /* Heap is empty but dict is not - rebuild */
        rebuildHeap();
        if (!heap.empty()) {
            return heap[0].key;
        }
    }
    
    throw std::out_of_range("PriorityDictionary corrupted");
}

template<typename K, typename V>
void PriorityDictionary<K, V>::set(const K& key, const V& val) {
    dict[key] = val;
    
    /* Rebuild heap if too many deleted items */
    if (heap.size() > 2 * dict.size()) {
        rebuildHeap();
    } else {
        /* Add new pair to heap */
        HeapPair newPair = {val, key};
        heap.push_back(newPair);
        
        size_t insertionPoint = heap.size() - 1;
        while (insertionPoint > 0) {
            size_t parent = (insertionPoint - 1) / 2;
            if (!(newPair < heap[parent])) {
                break;
            }
            heap[insertionPoint] = heap[parent];
            insertionPoint = parent;
        }
        heap[insertionPoint] = newPair;
    }
}

template<typename K, typename V>
V PriorityDictionary<K, V>::get(const K& key) const {
    auto it = dict.find(key);
    if (it == dict.end()) {
        throw std::out_of_range("Key not found in PriorityDictionary");
    }
    return it->second;
}

template<typename K, typename V>
bool PriorityDictionary<K, V>::contains(const K& key) const {
    return dict.find(key) != dict.end();
}

template<typename K, typename V>
void PriorityDictionary<K, V>::remove(const K& key) {
    dict.erase(key);
}

template<typename K, typename V>
K PriorityDictionary<K, V>::pop() {
    K key = smallest();
    dict.erase(key);
    return key;
}

template<typename K, typename V>
void PriorityDictionary<K, V>::rebuildHeap() {
    heap.clear();
    for (auto it = dict.begin(); it != dict.end(); ++it) {
        HeapPair pair = {it->second, it->first};
        heap.push_back(pair);
    }
    
    /* Sort to create heap */
    std::sort(heap.begin(), heap.end());
}

} // namespace mol2db2

#endif // PRIODICT_H
