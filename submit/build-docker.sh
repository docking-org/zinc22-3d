source /pyenv/bin/activate
cd /tmp
$TAUOMERIZE_PROTONATE_EXE -H 7.4 < /data/input.smi | $PROTOMER_STEREOCENTERS_EXE | python ${DOCKBASE}/ligand/generate/build_ligands_docker.py