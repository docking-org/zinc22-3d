#include "mol2.h"
#include <iostream>
#include <cassert>
#include <sstream>
#include <cmath>

using namespace mol2db2;

void test_mol2_creation() {
    Mol2 mol;
    assert(mol.name == "fake");
    assert(mol.xyzCount == -1);
    assert(mol.atomNum.empty());
    
    std::cout << "✓ test_mol2_creation passed" << std::endl;
}

void test_mol2_copy() {
    Mol2 mol1;
    mol1.name = "test";
    mol1.protName = "protein";
    mol1.atomNum.push_back(1);
    mol1.atomNum.push_back(2);
    
    Mol2 mol2 = mol1.copy();
    
    assert(mol2.name == "test");
    assert(mol2.protName == "protein");
    assert(mol2.atomNum.size() == 2);
    assert(mol2.atomNum[0] == 1);
    assert(mol2.atomNum[1] == 2);
    
    // Verify deep copy
    mol1.atomNum[0] = 99;
    assert(mol2.atomNum[0] == 1);
    
    std::cout << "✓ test_mol2_copy passed" << std::endl;
}

void test_atom_converter() {
    AtomConverter converter(nullptr);  // Use defaults
    
    Mol2 mol;
    mol.atomNum.push_back(1);
    mol.atomName.push_back("C1");
    mol.atomType.push_back("C.ar");
    mol.atomBonds.push_back(std::vector<std::pair<int, std::string>>());
    
    int dockType = converter.convertMol2atomNum(mol, 1);
    assert(dockType == 1);  // C.ar -> 1
    
    std::cout << "✓ test_atom_converter passed" << std::endl;
}

void test_color_table() {
    ColorTable colorTable(nullptr);  // Use defaults
    
    assert(colorTable.defaultColor == "neutral");
    assert(colorTable.colorInts["positive"] == 1);
    assert(colorTable.colorInts["negative"] == 2);
    assert(colorTable.colorInts["neutral"] == 7);
    
    std::cout << "✓ test_color_table passed" << std::endl;
}

void test_solv_data() {
    SolvData solv;
    
    assert(solv.totalAtoms == 0);
    assert(solv.totalCharge == 0.0);
    assert(solv.charge.empty());
    
    // Add some data
    solv.totalAtoms = 2;
    solv.charge.push_back(0.5);
    solv.charge.push_back(-0.5);
    
    assert(solv.charge.size() == 2);
    assert(std::fabs(solv.charge[0] - 0.5) < 0.0001);
    
    std::cout << "✓ test_solv_data passed" << std::endl;
}

void test_mol2_from_text() {
    std::vector<std::string> mol2text = {
        "@<TRIPOS>MOLECULE",
        "TEST",
        "2 1 0 0 0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
        "1 C1 0.0 0.0 0.0 C.ar 1 <0> -0.15",
        "2 C2 1.4 0.0 0.0 C.ar 1 <0> -0.15",
        "@<TRIPOS>BOND",
        "1 1 2 ar"
    };
    
    Mol2 mol(nullptr, nullptr, &mol2text);
    
    assert(mol.name == "TEST");
    assert(mol.atomNum.size() == 2);
    assert(mol.atomNum[0] == 1);
    assert(mol.atomNum[1] == 2);
    assert(mol.bondNum.size() == 1);
    assert(mol.xyzCount == 0);  // One conformation
    
    std::cout << "✓ test_mol2_from_text passed" << std::endl;
}

void test_bonded_to() {
    std::vector<std::string> mol2text = {
        "@<TRIPOS>MOLECULE",
        "ETHANE",
        "8 7 0 0 0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
        "1 C1 0.0 0.0 0.0 C.3 1 <0> 0.0",
        "2 C2 1.5 0.0 0.0 C.3 1 <0> 0.0",
        "3 H1 -0.5 0.5 0.5 H 1 <0> 0.0",
        "4 H2 -0.5 -0.5 0.5 H 1 <0> 0.0",
        "5 H3 -0.5 0.5 -0.5 H 1 <0> 0.0",
        "6 H4 2.0 0.5 0.5 H 1 <0> 0.0",
        "7 H5 2.0 -0.5 0.5 H 1 <0> 0.0",
        "8 H6 2.0 0.5 -0.5 H 1 <0> 0.0",
        "@<TRIPOS>BOND",
        "1 1 2 1",
        "2 1 3 1",
        "3 1 4 1",
        "4 1 5 1",
        "5 2 6 1",
        "6 2 7 1",
        "7 2 8 1"
    };
    
    Mol2 mol(nullptr, nullptr, &mol2text);
    
    // C1 (atom 1) should be bonded to C2 (atom 2)
    assert(mol.bondedTo(1, "C", 1, nullptr, false, nullptr));
    
    // C1 should be bonded to H at distance 1
    assert(mol.bondedTo(1, "H", 1, nullptr, false, nullptr));
    
    // C1 should be bonded to C at distance 1
    assert(mol.bondedTo(1, "C.3", 1, nullptr, false, nullptr));
    
    // H3 (atom 3) should be bonded to C at distance 1
    assert(mol.bondedTo(3, "C", 1, nullptr, false, nullptr));
    
    // H3 should be bonded to H at distance 2 (through C1)
    assert(mol.bondedTo(3, "H", 2, nullptr, false, nullptr));
    
    std::cout << "✓ test_bonded_to passed" << std::endl;
}

void test_rmsd() {
    std::vector<std::string> mol2text = {
        "@<TRIPOS>MOLECULE",
        "TEST",
        "2 1 0 0 0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
        "1 C1 0.0 0.0 0.0 C.ar 1 <0> 0.0",
        "2 C2 1.0 0.0 0.0 C.ar 1 <0> 0.0",
        "@<TRIPOS>BOND",
        "1 1 2 ar"
    };
    
    Mol2 mol(nullptr, nullptr, &mol2text);
    
    // Add a second conformation
    std::vector<std::array<double, 3>> conf2 = {
        {0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}  // Moved second atom
    };
    mol.atomXyz.push_back(conf2);
    mol.xyzCount = 2;
    
    double rmsd = mol.getRMSD(0, 1);
    
    // RMSD should be sqrt((0^2 + 1^2) / 2) = sqrt(0.5) ≈ 0.707
    assert(std::fabs(rmsd - 0.707) < 0.01);
    
    std::cout << "✓ test_rmsd passed" << std::endl;
}

void test_keep_confs_only() {
    Mol2 mol;
    mol.xyzCount = 5;
    
    for (int i = 0; i < 5; i++) {
        mol.atomXyz.push_back({{0.0, 0.0, 0.0}});
        mol.inputEnergy.push_back(i * 10.0);
        mol.inputHydrogens.push_back(0);
    }
    
    mol.keepConfsOnly(1, 4);  // Keep conformations 1, 2, 3
    
    assert(mol.xyzCount == 3);
    assert(mol.atomXyz.size() == 3);
    assert(mol.inputEnergy.size() == 3);
    assert(std::fabs(mol.inputEnergy[0] - 10.0) < 0.001);
    assert(std::fabs(mol.inputEnergy[2] - 30.0) < 0.001);
    
    std::cout << "✓ test_keep_confs_only passed" << std::endl;
}

int main() {
    std::cout << "Running MOL2 tests..." << std::endl;
    
    test_mol2_creation();
    test_mol2_copy();
    test_atom_converter();
    test_color_table();
    test_solv_data();
    test_mol2_from_text();
    test_bonded_to();
    test_rmsd();
    test_keep_confs_only();
    
    std::cout << "\n✓✓✓ All MOL2 tests passed! ✓✓✓\n" << std::endl;
    return 0;
}
