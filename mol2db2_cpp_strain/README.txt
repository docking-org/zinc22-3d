================================================================================
mol2db2 C++ — Conversion Documentation
================================================================================

OVERVIEW
================================================================================

mol2db2 converts MOL2 molecular structure files to hierarchical DB2 format
for use with DOCK molecular docking software.

This C++ implementation is a faithful port of the Python reference:
mol2db2_py3_strain (Shoichet Lab, UCSF), including Benjamin Tingle's
_getRigidStructures / _countPositions rewrite and UCSF strain energy support.
Output is functionally equivalent to the Python version.


ORIGINAL PYTHON IMPLEMENTATION
================================================================================

Author:       Ryan G. Coleman
Affiliation:  Brian K. Shoichet Lab, University of California, San Francisco
New algorithm (_getRigidStructures, _countPositions, buckets2):
              Benjamin Tingle, February 2021
Strain energy integration:
              J.K. Lyu (jklyu), May 2020

Citation: If you use this software, please cite the original Python work.


C++ CONVERSION
================================================================================

Translator:   Claude (Anthropic AI Assistant)
              v1.0.0: Claude Sonnet 4  (February 2026)
              v1.1.0: Claude Sonnet 4.6 (February 2026)
Language:     C++17
Method:       Human-guided AI translation with iterative refinement
Human guidance: Trent Balius, Brendan Hall


ALGORITHM
================================================================================

The C++ implementation uses the same algorithm as mol2db2_py3_strain:

1. _getRigidStructures: identifies rigidly-bonded atom groups via cycle
   detection, rigid-bond virtual cycles, and intersection-graph DFS.

2. _countPositions (buckets2): hash-bucket spatial clustering of per-atom
   positions across conformers; splits output conformations at rotatable
   bond boundaries using the rigidStructures map.

3. decideDistanceRules: fast clash detection (same-type atom pairs only).

4. Strain fields: two values on S-lines (totalStrain, maxStrain).

5. No clouds (D-lines): divisive clustering is not used.

6. T-lines (color table) written only when non-default.

Known output difference vs Python:
  X-line ordering within same-sized conformer clusters may differ due to
  hash-map vs dict insertion order.  The same set of coordinates is always
  present.  


STRAIN ENERGY SUPPORT
================================================================================

MOL2 input — the parser recognises either format:
  mmff94s_NoEstat = <value>          (old, stored as inputEnergy)
  Strain = <totalStrain> <maxStrain> (new, stored as inputTotalStrain/inputMaxStrain)

S-line output:
  S <setNum> <lines> <confs> <broken> <hydro> <+totalStrain> <+maxStrain>
  Falls back to 9999.999 if no strain data is present.

Hydrogen rotation:
  Both rotateHydrogens() and rotateHydrogens_teb() propagate
  inputTotalStrain and inputMaxStrain to all generated conformers.

Python writeMol2File integration:
  mol2.py's writeMol2File() writes "Strain = total max" when
  inputTotalStrain is populated, so strain data survives piping through
  the binary via stdin with no post-hoc line replacement needed.


STDIN / STDOUT SUPPORT
================================================================================

  -m -    Read mol2 data from stdin (plain text, not gzip)
  -o -    Write db2 data to stdout (plain text, not gzip)

Or inline:
  import io, subprocess
  sio = io.StringIO()
  mol2_obj.writeMol2File(sio)
  db2_data = subprocess.run(
      [MOL2DB2_CPP, "-m", "-", "-s", solvfile, "-d", clashfile, "-o", "-"],
      input=sio.getvalue(), capture_output=True, text=True, check=True,
  ).stdout


USAGE
================================================================================

./bin/mol2db2 [options]

Options:
  -v, --verbose          Verbose output
  -m, --mol2 FILE        Input MOL2 file (use "-" to read from stdin)
  -s, --solv FILE        Solvation file
  -o, --db FILE          Output DB2.gz file (use "-" for plain-text stdout)
  -n, --namemol2 FILE    Name file
  -t, --atomtype FILE    Atom type conversion file
  -c, --colortable FILE  Color table file
  -d, --clash FILE       Clash parameter file
  -y, --hydrogens FILE   Hydrogen parameter file
  -r, --noreseth         Don't reset planar hydrogens
  -z, --norotateh        Don't rotate terminal hydrogens
  -x, --disttol FLOAT    Distance tolerance (default: 0.001)
  --limitset NUM         Maximum number of sets (default: 2500000)
  --limitconf NUM        Maximum number of conformations (default: 200000)
  --limitcoord NUM       Maximum number of coordinates (default: 1000000)
  --covalent             Process covalent ligands (SiH3 dummy removal)
  -a, --time             Show timing information
  --help                 Show help message

Examples:
  ./bin/mol2db2 -m db.mol2.gz -s db.solv -o db.db2.gz -v
  ./bin/mol2db2 -m db.mol2.gz -s db.solv -d clashfile.txt -o db.db2.gz


BUILD AND INSTALLATION
================================================================================

Dependencies:
  - C++17 compiler (GCC 7+ or Clang 5+)
  - Eigen3 (place eigen-3.4.0/ inside mol2db2_cpp_strain/)
  - zlib (sudo apt-get install zlib1g-dev)

Build:
  make
  # Binary is produced at bin/mol2db2

Install to your local bin (no sudo needed):
  make install
  # Copies to ~/.local/bin/mol2db2 by default.
  # Make sure ~/.local/bin is in your PATH:
  #   export PATH="$HOME/.local/bin:$PATH"

Install to a custom location:
  make install PREFIX=/path/to/your/dir
  # Copies to /path/to/your/dir/bin/mol2db2

Use without installing (reference by full path):
  /path/to/mol2db2_cpp_strain/bin/mol2db2 [options]

Debug build:
  CXXFLAGS="-std=c++17 -g -O0" make

Architecture-specific optimisation (faster, not portable):
  CXXFLAGS="-std=c++17 -O3 -march=native" make

Run tests:
  make test


FILE ORGANIZATION
================================================================================

src/core/         geometry, unionfind, combinatorics, priodict, floydwarshall
src/algorithms/   shortestpaths, buckets, buckets2, pca, divisive_clustering
src/molecular/    mol2, clash, hydrogens
src/hierarchy/    hierarchy
src/main/         mol2db2.cpp
mol2db2_cpp.py    Python subprocess wrapper for drop-in pipeline integration
tests/            Unit tests (geometry, unionfind, mol2)


PYTHON FILES CONVERTED
================================================================================

geometry.py              → geometry.h/cpp
combinatorics.py         → combinatorics.h (template)
unionfind2.py            → unionfind.h/cpp
priodict.py              → priodict.h (template)
floydwarshall.py         → floydwarshall.h/cpp
shortestpaths.py         → shortestpaths.h/cpp
buckets.py               → buckets.h/cpp
buckets2.py              → buckets2.h/cpp
pca.py                   → pca.h/cpp (using Eigen3)
divisive_clustering.py   → divisive_clustering.h/cpp
clash.py                 → clash.h/cpp
hydrogens.py             → hydrogens.h/cpp
hierarchy.py             → hierarchy.h/cpp
mol2db2.py               → mol2db2.cpp
mol2.py + solv.py + sybyl2dock.py + atom_color_table.py → mol2.h/cpp

mol2.py writeMol2File()  updated to write Strain lines
mol2db2_cpp.py           new Python subprocess wrapper
mol2hydroxyl.py          not converted (functionality in main program)


PERFORMANCE
================================================================================

Known future optimisation:
  confClusters uses std::map<vector<int>, ...> (O(log n) lookups).
  Switching to an unordered_map with a vector hash would give O(1).


KNOWN LIMITATIONS
================================================================================

1. X-line ordering differs from Python reference (cosmetic, no data loss)

2. No recursive subdivision on TooBigError
   The binary exits when size limits are exceeded.  Use --limitconf and
   --limitset to split large jobs upstream.

3. mol2hydroxyl standalone utility not converted


VERSION HISTORY
================================================================================

Version 1.1.0 (February 2026)
  - Brought C++ into full alignment with mol2db2_py3_strain Python reference
  - Replaced old cpp algorithm with the Python algorithm as the sole mode
    (removed --cpp / --python flags; now always uses Python algorithm)
  - Added buckets2.h/cpp: port of buckets2.py hash-bucket clustering
  - Added getRigidStructures(): cycle detection + rigid bond virtual cycles
  - Added countPositionsPython(): per-atom position clustering with conf
    splitting at rotatable bond boundaries
  - Added decideDistanceRules(): fast same-type-pair clash detection
  - Added strain energy fields throughout (inputTotalStrain, inputMaxStrain)
  - Fixed: std::showpos was not cleared after float fields, causing spurious
    + signs on integers in A, X, R, C, S lines
  - Fixed: posCount computed as per-atom local cluster count (not global
    confNum delta), fixing missing R-lines
  - Fixed: strain not propagated to hydrogen-rotation-generated conformers
  - Added: stdin support (-m -) and stdout support (-o -)
  - Added: mol2db2_cpp.py Python subprocess wrapper
  - Updated: mol2.py writeMol2File() writes "Strain = total max"
  - Updated: T-lines (color table) written only when non-default

Version 1.0.0 (February 2026)
  - Initial C++ conversion based on pre-strain Python version
  - Core utilities, algorithms, molecular structures: 100% faithful
  - Hierarchy: used deprecated findConformations algorithm only


LICENSING
================================================================================

This C++ implementation is a derivative work of the original Python code.
The original licensing terms of the Python implementation apply.


================================================================================
Last Updated: February 2026  |  Document Version: 1.2  |  Code Version: 1.1.0
================================================================================
