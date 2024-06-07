#!/bin/bash

#set install path to current directory
INSTALL_PATH=$(pwd)
export CONDA_ROOT=which conda

#if conda is not installed, exit    
if [ -z $CONDA_ROOT ]; then
    echo "Conda is not installed. Please install conda and try again."
    exit 1
fi

mkdir -p $INSTALL_PATH/soft/DOCK-latest
mkdir licenses

#move all folders except submit
mv 3D $INSTALL_PATH/soft/DOCK-latest/ligand
mv amsol $INSTALL_PATH/soft/DOCK-latest/ligand
mv db22mol2 $INSTALL_PATH/soft/DOCK-latest/ligand
mv finish $INSTALL_PATH/soft/DOCK-latest/ligand
mv generate $INSTALL_PATH/soft/DOCK-latest/ligand
mv mol2db $INSTALL_PATH/soft/DOCK-latest/ligand
mv mol2db2 $INSTALL_PATH/soft/DOCK-latest/ligand
mv mol2db2_py3 $INSTALL_PATH/soft/DOCK-latest/ligand
mv mol2db2_py3_strain $INSTALL_PATH/soft/DOCK-latest/ligand
mv omega $INSTALL_PATH/soft/DOCK-latest/ligand/ligand
mv protonate $INSTALL_PATH/soft/DOCK-latest/ligand
mv rdkit $INSTALL_PATH/soft/DOCK-latest/ligand
mv reactor $INSTALL_PATH/soft/DOCK-latest/ligand
mv strain $INSTALL_PATH/soft/DOCK-latest/ligand1
