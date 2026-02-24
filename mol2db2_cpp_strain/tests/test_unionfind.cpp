#include "unionfind.h"
#include <iostream>
#include <cassert>
#include <algorithm>

using namespace mol2db2;

void test_basic_union_find() {
    UnionFind uf;
    
    // Create elements
    uf.find(0);
    uf.find(1);
    uf.find(2);
    
    // Initially all different
    assert(uf.different(0, 1));
    assert(uf.different(1, 2));
    assert(uf.different(0, 2));
    
    // Union 0 and 1
    uf.unionSets(0, 1);
    assert(!uf.different(0, 1));
    assert(uf.different(1, 2));
    
    // Union 1 and 2 (should connect all three)
    uf.unionSets(1, 2);
    assert(!uf.different(0, 1));
    assert(!uf.different(1, 2));
    assert(!uf.different(0, 2));
    
    std::cout << "✓ test_basic_union_find passed" << std::endl;
}

void test_getList() {
    UnionFind uf;
    
    // Create two groups: {0,1,2} and {3,4}
    uf.find(0);
    uf.find(1);
    uf.find(2);
    uf.find(3);
    uf.find(4);
    
    uf.unionSets(0, 1);
    uf.unionSets(1, 2);
    uf.unionSets(3, 4);
    
    // Get list for element 0
    std::vector<int> list0 = uf.getList(0);
    assert(list0.size() == 3);
    assert(std::find(list0.begin(), list0.end(), 0) != list0.end());
    assert(std::find(list0.begin(), list0.end(), 1) != list0.end());
    assert(std::find(list0.begin(), list0.end(), 2) != list0.end());
    
    // Get list for element 3
    std::vector<int> list3 = uf.getList(3);
    assert(list3.size() == 2);
    assert(std::find(list3.begin(), list3.end(), 3) != list3.end());
    assert(std::find(list3.begin(), list3.end(), 4) != list3.end());
    
    std::cout << "✓ test_getList passed" << std::endl;
}

void test_toLists() {
    UnionFind uf;
    
    // Create three groups: {0,1}, {2,3,4}, {5}
    for (int i = 0; i < 6; i++) {
        uf.find(i);
    }
    
    uf.unionSets(0, 1);
    uf.unionSets(2, 3);
    uf.unionSets(3, 4);
    
    std::vector<std::vector<int>> lists = uf.toLists();
    
    assert(lists.size() == 3);
    
    // Check sizes
    std::vector<int> sizes;
    for (const auto& list : lists) {
        sizes.push_back(list.size());
    }
    std::sort(sizes.begin(), sizes.end());
    assert(sizes[0] == 1);  // {5}
    assert(sizes[1] == 2);  // {0,1}
    assert(sizes[2] == 3);  // {2,3,4}
    
    std::cout << "✓ test_toLists passed" << std::endl;
}

void test_path_compression() {
    UnionFind uf;
    
    // Create a long chain: 0-1-2-3-4
    for (int i = 0; i < 5; i++) {
        uf.find(i);
    }
    
    uf.unionSets(0, 1);
    uf.unionSets(1, 2);
    uf.unionSets(2, 3);
    uf.unionSets(3, 4);
    
    // Path compression should happen on find
    int root = uf.find(4);
    
    // All should now point directly to root
    assert(uf.find(0) == root);
    assert(uf.find(1) == root);
    assert(uf.find(2) == root);
    assert(uf.find(3) == root);
    assert(uf.find(4) == root);
    
    std::cout << "✓ test_path_compression passed" << std::endl;
}

void test_unionfind_attach() {
    UnionFindAttach ufa;
    
    // Create elements with attached data
    std::set<int> data1 = {10, 20};
    std::set<int> data2 = {30, 40};
    
    ufa.find(0, &data1);
    ufa.find(1, &data2);
    
    // Check attached data before union
    std::set<int> attached0 = ufa.getAttached(0);
    assert(attached0.size() == 2);
    assert(attached0.count(10) == 1);
    assert(attached0.count(20) == 1);
    
    // Union should merge attached data
    ufa.unionSets(0, 1);
    
    std::set<int> attachedMerged = ufa.getAttached(0);
    assert(attachedMerged.size() == 4);
    assert(attachedMerged.count(10) == 1);
    assert(attachedMerged.count(20) == 1);
    assert(attachedMerged.count(30) == 1);
    assert(attachedMerged.count(40) == 1);
    
    std::cout << "✓ test_unionfind_attach passed" << std::endl;
}

void test_large_dataset() {
    UnionFind uf;
    
    // Create 1000 elements
    for (int i = 0; i < 1000; i++) {
        uf.find(i);
    }
    
    // Union into 10 groups of 100
    for (int group = 0; group < 10; group++) {
        for (int i = 1; i < 100; i++) {
            uf.unionSets(group * 100, group * 100 + i);
        }
    }
    
    std::vector<std::vector<int>> lists = uf.toLists();
    assert(lists.size() == 10);
    
    for (const auto& list : lists) {
        assert(list.size() == 100);
    }
    
    std::cout << "✓ test_large_dataset passed" << std::endl;
}

int main() {
    std::cout << "Running UnionFind tests..." << std::endl;
    
    test_basic_union_find();
    test_getList();
    test_toLists();
    test_path_compression();
    test_unionfind_attach();
    test_large_dataset();
    
    std::cout << "\n✓✓✓ All UnionFind tests passed! ✓✓✓\n" << std::endl;
    return 0;
}
