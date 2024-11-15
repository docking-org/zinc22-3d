#base=/change/this
#source /path/to/conda/bin/activate env name

export BINDIR=$base/submit
export SHRTCACHE=//scratch
export LONGCACHE=/scratch
export SOFT_HOME=$base/soft
export LICENSE_HOME=$base/licenses
export PATH=$PATH:$base/submit
export SKIP_DELETE=false

export DOCK_VERSION=DOCK-latest
export CORINA_VERSION=corina-2025
# export PYENV_VERSION=lig_build_py3
export JCHEM_VERSION=jchem-19.15.r4
export OPENBABEL_VERSION=openbabel-3.1.1
export EXTRALIBS_VERSION=extralibs

export SUBMIT_MODE=TEST_LOCAL
export INPUT_FILE=/path/to/in/smiles
export OUTPUT_DEST=/path/to/out/directory
export BUILD_MOL2="false"

export LD_LIBRARY_PATH=/path/to/miniconda3/envs/envname/lib:$LD_LIBRARY_PATH
export CHEMAXON_PATH=$LONGCACHE/build_3d_common_$(whoami)/$JCHEM_VERSION
export CHEMAXON_LICENSE_PATH=$LICENSE_HOME/jchem-license.cxl
export PATH=$PATH:$CHEMAXON_PATH/bin

export COMMON_DIR=$LONGCACHE/build_3d_common_$(whoami)
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$COMMON_DIR/$OPENBABEL_VERSION/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$COMMON_DIR/$EXTRALIBS_VERSION
export CORINA_ROOT=$COMMON_DIR/build_3d_common_$(whoami)/$CORINA_VERSION
export PATH="$COMMON_DIR/$CORINA_VERSION:${PATH}"

export OBABELBASE=$COMMON_DIR/$OPENBABEL_VERSION
OB_VER=$(echo $OPENBABEL_VERSION | cut -d'-' -f2-)
export BABEL_LIBDIR=$COMMON_DIR/$OPENBABEL_VERSION/lib/openbabel/$OB_VER
export BABEL_DATADIR=$COMMON_DIR/$OPENBABEL_VERSION/share/openbabel/$OB_VER
export PATH="${PATH}:${OBABELBASE}/bin"

export OE_LICENSE=$LICENSE_HOME/oe-license.txt
export CHEMAXON_PATH=$COMMON_DIR/$JCHEM_VERSION
export CHEMAXON_LICENSE_URL=$LICENSE_HOME/jchem-license.cxl
export PATH="$PATH:$CHEMAXON_PATH/bin"