#!/bin/sh
USAGE="Usage: ${0} [OPTIONS] <SOURCE_FILE>

  Options:
    -h, --help - Display this message and exit
    -H, --pH <PH_LEVELS> - A quoted, space separated list of pH levels to 
                           build tautomers/protomers at
    -s, --single - Build a single db/db2 file instead of separate files
                   for each protomer
    -n, --name <NAME> - Override database name 
    -d, --dir <DIR> - Working directory
    -c, --covalent - Build a covalent library instead of standard
    -3, --3d - Use provided 3D structures (implies --pre-tautomerized)
    --no-limit-confs-by-hydrogens - Don't limit # conformations by # rotatable hydrogens
    --pre-tautomerized - Treat input file as pre-generated tautomers
    --permissive-taut-prot - Use lower tautomer and protomer cutoffs
    --no-conformations - Skip generating multiple rigid fragment conformations
    --no-db - Skip building db files
    --no-db2 - Skip building db2 files
    --no-solvation - Don't save solvation files
    --no-mol2 - Don't save mol2 files
    --save-table - Save the full protomer table
    --bad-charges <BAD_CHARGE_FILE> - List of bad protonation patterns to
                                      exclude
    --debug - Extra debugging output

  Overrideable Sub-programs:
    TAUOMERIZE_PROTONATE_EXE - Generate (multiple) tautomerized and protonated
                               variants of the input substances at a pH level
    PROTOMER_COALESE_EXE - Filter and merge protomers over pH levels
    PROTOMER_STEREOCENTERS_EXE - Expand any new stereocenters from protonation
    EMBED_PROTOMERS_3D_EXE - Create 3D mol2 files for each protomer 
                             (names should JUST be the line number of the 
                              protomer without any extension)
    PREPARE_NAME_EXE - Write the name.txt file to build a db2 file with
    SOLVATION_EXE - Calculation solvation for a given mol2 file
    GENERATE_CONFORMATIONS_EXE - Generate heirarchy conformations 
        GENERATE_RIGID_FRAGMENT_CONFORMATIONS_EXE - Generate standard 
                                                    heirarchy conformations
        GENERATE_COVALENT_CONFORMATIONS_EXE - Generate covalent heirarchy
                                              conformations
    BUILD_DB2_EXE - Generate a db2 file from conformations
    BUILD_DB_EXE - Genearte a db file from conformations
"

set -e

DOCKBASE="${DOCKBASE-$( dirname $( dirname $( dirname $BASH_SOURCE ) ) )}"

DEBUG="${DEBUG-}"
PROTOMER_PH_LEVELS=( ${PROTOMER_PH_LEVELS-7.4 6.4 8.4} )
COVALENT="${COVALENT-no}"
STORE_PROTOMERS="${STORE_PROTOMERS}"
USE_3D="${USE_3D-no}"
CREATE_TAUTOMERS="${CREATE_TAUTOMERS-yes}"
HYDROGEN_DEPENDENT_CONFORMATION_LIMITS="${HYDROGEN_DEPENDENT_CONFORMATION_LIMITS-yes}"
CREATE_CONFORMATIONS="${CREATE_CONFORMATIONS-yes}"
BUILD_DB_FILES="${BUILD_DB_FILES-yes}"
BUILD_DB2_FILES="${BUILD_DB2_FILES-yes}"
COPY_MOL2_FILES="${COPY_MOL2_FILES-yes}"
COPY_SOLV_FILES="${COPY_SOLV_FILES-yes}"
SINGLE_DATABASE="${SINGLE_DATABASE-no}"
WRITE_PROTOMER_TABLE="${WRITE_PROTOMER_TABLE-no}"
DBNAME="${NAME-}"
TASK_DIR="${TASK_DIR-$( pwd )}"
BAD_PROTOMER_CHARGES="${BAD_PROTOMER_CHARGES-${DOCKBASE}/ligand/protonate/data/bad-charges.txt}"

TAUOMERIZE_PROTONATE_EXE="${TAUOMERIZE_PROTONATE_EXE-${DOCKBASE}/ligand/protonate/tautprot_cxcalc.sh}"
#PROTOMER_COALESE_EXE="${PROTOMER_COALESE_EXE-${DOCKBASE}/ligand/protonate/coalese.py --sort --limits=1 --filter=${BAD_PROTOMER_CHARGES}}"
PROTOMER_COALESE_EXE="${PROTOMER_COALESE_EXE-${DOCKBASE}/ligand/protonate/coalese.py3.py --sort --limits=1 --filter=${BAD_PROTOMER_CHARGES}}"
#PROTOMER_STEREOCENTERS_EXE="${PROTOMER_STEREOCENTERS_EXE-${DOCKBASE}/ligand/protonate/expand-new-stereocenters.py}"
PROTOMER_STEREOCENTERS_EXE="${PROTOMER_STEREOCENTERS_EXE-${DOCKBASE}/ligand/protonate/expand-new-stereocenters.py3.py}"
UNEMBED_PROTOMERS_2D_EXE="${EMBED_PROTOMERS_3D_EXE-${DOCKBASE}/ligand/3D/unembed2d_molconvert.sh}"
EMBED_PROTOMERS_3D_EXE="${EMBED_PROTOMERS_3D_EXE-${DOCKBASE}/ligand/3D/embed3d_molconvert.sh}"
EXPAND_3D_EXE="${EXPAND_3D_EXE-${DOCKBASE}/ligand/3D/expand3d_molconvert.sh}"
#PREPARE_NAME_EXE="${PREPARE_NAME_EXE-${DOCKBASE}/ligand/generate/prepare.py}"
PREPARE_NAME_EXE="${PREPARE_NAME_EXE-${DOCKBASE}/ligand/generate/prepare.py3.py}"
SOLVATION_EXE="${SOLVATION_EXE-${DOCKBASE}/ligand/amsol/calc_solvation.csh}"
#COUNT_ROTATABLE_HYDROGENS_EXE="${COUNT_ROTATABLE_HYDROGENS_EXE-${DOCKBASE}/ligand/mol2db2/hydrogens.py}"
COUNT_ROTATABLE_HYDROGENS_EXE="${COUNT_ROTATABLE_HYDROGENS_EXE-${DOCKBASE}/ligand/mol2db2_py3/hydrogens.py}"
#GENERATE_RIGID_FRAGMENT_CONFORMATIONS_EXE="${GENERATE_RIGID_FRAGMENT_CONFORMATIONS_EXE-${DOCKBASE}/ligand/omega/omega_db2.e.py}"
GENERATE_RIGID_FRAGMENT_CONFORMATIONS_EXE="${GENERATE_RIGID_FRAGMENT_CONFORMATIONS_EXE-${DOCKBASE}/ligand/omega/omega_db2.e.py3.py}"
UCSF_STRAIN_EXE="${UCSF_STRAIN_EXE-${DOCKBASE}/ligand/strain/add_strain_noH_info.py}"
GENERATE_COVALENT_CONFORMATIONS_EXE="${COVALENT_CONFORMATIONS_EXE-${DOCKBASE}/ligand/omega/omega_warhead.py}"
#BUILD_DB2_STANDARD="${BUILD_DB2_STANDARD-${DOCKBASE}/ligand/mol2db2/mol2db2.py -v --solv=output.solv}"
#BUILD_DB2_STANDARD="${BUILD_DB2_STANDARD-${DOCKBASE}/ligand/mol2db2/mol2db2.py -v --solv=output.solv --clash=${DOCKBASE}/ligand/mol2db2/clashfile.txt}"
#BUILD_DB2_STANDARD="${BUILD_DB2_STANDARD-${DOCKBASE}/ligand/mol2db2_py3/mol2db2.py -v --solv=output.solv --clash=${DOCKBASE}/ligand/mol2db2/clashfile.txt}"
BUILD_DB2_STANDARD="${BUILD_DB2_STANDARD-${DOCKBASE}/ligand/mol2db2_py3_strain/mol2db2.py -v --solv=output.solv --clash=${DOCKBASE}/ligand/mol2db2/clashfile.txt}"
BUILD_DB2_COVALENT="${BUILD_DB2_COVALENT-${BUILD_DB2_STANDARD} --covalent}"
BUILD_DB_EXE="${BUILD_DB_EXE-${DOCKBASE}/ligand/mol2db/mol2db ${DOCKBASE}/ligand/mol2db/data/inhier}"

while [[ "$#" > 0 ]] ; do
    ARG="${1}"
    VALUE="${2}"
    shift 1
    case "${ARG}" in
        -h|--help)
            echo "${USAGE}" 1>&2
            exit -1
            ;;
        -H|--pH)
            PROTOMER_PH_LEVELS=( $( echo $VALUE | sed 's/,/ /g' ) )
            shift 1
            ;;
        -s|--single)
            SINGLE_DATABASE="yes"
            ;;
        -n|--name)
            DBNAME="${VALUE}"
            shift 1
            ;;
        -d|--dir)
            TASK_DIR="${VALUE}"
            shift 1
            ;;
        -C|--covalent)
            COVALENT="yes"
            ;;
        -3|--3d)
            USE_3D="yes"
            CREATE_TAUTOMERS="no"
            ;;
        --no-db)
            BUILD_DB_FILES="no"
            ;;
        --no-db2)
            BUILD_DB2_FILES="no"
            ;;
        --no-solv)
            COPY_SOLV_FILES="no"
            ;;
        --no-mol2)
            COPY_MOL2_FILES="no"
            ;;
        --save-table)
            SAVE_PROTOMER_TABLE="yes"
            ;;
        --bad-charges)
            BAD_PROTOMER_CHARGES="${VALUE}"
            shift 1
            ;;
        --no-conformations)
            CREATE_CONFORMATIONS="no"
            ;;
        --no-limit-confs-by-hydrogens)
            HYDROGEN_DEPENDENT_CONFORMATION_LIMITS="no"
            ;;
        --pre-tautomerized)
            CREATE_TAUTOMERS="no"
            ;;
        --permissive-taut-prot)
            echo "NOTICE: Using permissive tautomer and protomer thresholds. Could produce many protomers!" 1>&2
            export TAUT_PROT_CUTOFF=2
            export TAUTOMER_LIMIT=5
            export PROTOMER_LIMIT=5
            ;;
        --debug)
            DEBUG="yes"
            ;;
        *)
            if [ -z "${SOURCE_FILE}" ] ; then
                SOURCE_FILE="${ARG}"
            else
                echo "Invalid argument: ${ARG}" 1>&2
                echo "${USAGE}" 1>&2
            fi
    esac
done

if [ "${CREATE_CONFORMATIONS}" == "no" -a "${BUILD_DB_FILES}" == "yes" ] ; then
    echo "WARNING: Building DB files without conformations can hang!" 1>&2
    echo "Disabling DB preparation" 1>&2
    BUILD_DB_FILES="no"
fi

if [ -z "${SOURCE_FILE}" ] ; then
	echo "Error: source file not provided!" 1>&2
	exit -1
elif [ ! -e "${SOURCE_FILE}" ] ; then
	echo "Error: source file '${SOURCE_FILE}' does not exist!" 1>&2
	exit -1
else
	SOURCE_FILE="$( readlink -f ${SOURCE_FILE} )"
fi
if [ "${COVALENT}" == "yes" ] ; then
    GENERATE_CONFORMATIONS_EXE="${GENERATE_CONFORMATIONS_EXE-${GENERATE_COVALENT_CONFORMATIONS_EXE}}"
    BUILD_DB2_EXE="${BUILD_DB2_EXE-${BUILD_DB2_COVALENT}}"
    if [ "${BUILD_DB_FILES}" == "yes" ] ; then
        echo "Creation of covalent DB files not (yet) supported. DB files will not be build!" 1>&2
        BUILD_DB_FILES="no"
    fi
else
    GENERATE_CONFORMATIONS_EXE="${GENERATE_CONFORMATIONS_EXE-${GENERATE_RIGID_FRAGMENT_CONFORMATIONS_EXE}}"
    BUILD_DB2_EXE="${BUILD_DB2_EXE-${BUILD_DB2_STANDARD}}"
fi

BASE="$( basename "${SOURCE_FILE}" .gz )"

if [ "$( basename "${BASE}" .mol2 )" != "${BASE}" ] ; then
    echo "Input is a mol2 file using provided structures" 1>&2
    USE_3D="yes"
    CREATE_TAUTOMERS="no"
    if [ -z "${DBNAME}" ] ; then
        DBNAME="${BASE}"
        DBNAME="$( basename "${DBNAME}" .mol2 )"
    fi
elif [ -z "${DBNAME}" ] ; then
    DBNAME="${BASE}"
    DBNAME="$( basename "${DBNAME}" .smi )"
    DBNAME="$( basename "${DBNAME}" .ism )"
fi

if [ -z "$STORE_PROTOMERS" ]; then
    echo "STORE_PROTOMERS is not set! Will keep all results to finished directory" 1>&2
fi
BAD_PROTOMER_CHARGES="$( readlink -f "${BAD_PROTOMER_CHARGES}" )"

IN_PROGRESS="${IN_PROGRESS-${TASK_DIR}/working}"
FINISHED="${FINISHED-${TASK_DIR}/finished}"
FAILED="${FAILED-${TASK_DIR}/failed}"
ARCHIVE="${ARCHIVE-${TASK_DIR}/archive}"
UNLOADED="${UNLOADED-${TASK_DIR}/unloaded}"

PROTONATING_DIR="${PROTONATING_DIR-${IN_PROGRESS}/protonate}"
BUILD_DIR="${BUILD_DIR-${IN_PROGRESS}/building}"
CONFORMATIONS="${CONFORMATIONS-${IN_PROGRESS}/3D}"

CLEANED_SMILES="${IN_PROGRESS}/input-smiles.ism"
GENERATED_PROTOMER_SMILES="${GENERATED_PROTOMER_SMILES-${PROTONATING_DIR}/${DBNAME}-protomers.ism}"
EXPANDED_PROTOMER_SMILES="${EXPANDED_PROTOMER_SMILES-${PROTONATING_DIR}/${DBNAME}-protomers-expanded.ism}"
GENERATED_PROTOMER_SOLVATION="${EXPANDED_PROTOMER_SMILES-${IN_PROGRESS}/${DBNAME}-protomers.solv}"

PROTONATED_FILES=()

if ! [ -f $EXPANDED_PROTOMER_SMILES ]; then

    [ ! -z "${DEBUG}" ] && \
    echo "Creating $( basename "${GENERATED_PROTOMER_SMILES}" ) and $( basename "${GENERATED_PROTOMER_SOLVATION}" )" 1>&2
    mkdir -pv "${PROTONATING_DIR}" 1>&2
    echo "Storing results in ${FINISHED}" 1>&2
    echo "Working in ${IN_PROGRESS}" 1>&2

    pushd "${IN_PROGRESS}" 1>&2
    grep -v '^\s*$' "${SOURCE_FILE}" > "${CLEANED_SMILES}"

    echo "" 1>&2
    pushd "${PROTONATING_DIR}" 1>&2
    if [ "${CREATE_TAUTOMERS}" == "yes" ] ; then
        echo "Precomputing protomers for all compounds (pH: ${PROTOMER_PH_LEVELS[@]})" 1>&2
        INCOMMING_SMILES="${CLEANED_SMILES}"
        
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START all protonation | " && date
        for PH in "${PROTOMER_PH_LEVELS[@]}" ; do
            PROTONATED_FILE="${DBNAME}-protonated-${PH}.ism"
            if $TAUOMERIZE_PROTONATE_EXE -H "${PH}" > "${PROTONATED_FILE}" < "${INCOMMING_SMILES}" ; then
                echo "ph ${PH}:" $( cat "${PROTONATED_FILE}" | wc -l ) "protomers created" 1>&2
                PROTONATED_FILES=( "${PROTONATED_FILES[@]}" "${PROTONATED_FILE}" )
            else
                echo "pH ${PH} protomer generation failed!" 1>&2
                echo "${PH}" >> protonation-failures
                break
            fi
        done
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP all protonation | " && date
        
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START protonation post processing | " && date
        # Check if any protonation steps failed
        if [ -e "${PROTONATING_DIR}/protonation-failures" ] ; then
            echo "Protonation failures detected" | tee -a failure-reason 1>&2
            echo "Aborting!" 1>&2
            popd 1>&2
            popd 1>&2
            if [ -z "${DEBUG}" ] ; then
                mkdir -pv "${FAILED}" 1>&2
                mv -v ${PROTONATING_DIR}/* "${FAILED}" 1>&2
            fi
            exit -1
        fi
    else
        echo "Using inputs as protomers/tautomers. No processing done" 1>&2
        if [ "${USE_3D}" == "yes" ] ; then
            echo "Extracting 3D structures into protomers for future reference" 1>&2
            RAW_INCOMMING_SMILES="${PROTONATING_DIR}/${DBNAME}-raw.ism"
            INCOMMING_SMILES="${PROTONATING_DIR}/${DBNAME}-extracted.ism"
            if ! $UNEMBED_PROTOMERS_2D_EXE "${CLEANED_SMILES}" "${RAW_INCOMMING_SMILES}" ; then
                echo "Failed to extract smiles from provided file"
                popd 1>&2
                popd 1>&2
                if [ -z "${DEBUG}" ] ; then
                    mkdir -pv "${FAILED}" 1>&2
                    mv -v ${IN_PROGRESS}/* ${FAILED} 1>&2
                fi
                exit -1
            fi
            awk  -v NAME="${DBNAME}" '{
                if ( NF == 1 ) {
                    print $1, NAME "-" NR;
                } else {
                    print $0;
                }
        }' < "${RAW_INCOMMING_SMILES}" > "${INCOMMING_SMILES}" 
        else
            INCOMMING_SMILES="${SOURCE_FILE}"
        fi
        PROTONATED_FILE="${DBNAME}-protonated-manual.ism"
        PROTONATED_FILES=( "${PROTONATED_FILE}" )
        # Mock protomer selection
        awk '{ print $0, 100, 100, 100 }' < "${INCOMMING_SMILES}" > "${PROTONATED_FILE}"
        echo $( cut -d\   -f 2 "${PROTONATED_FILE}" | uniq | wc -l ) "substances and" \
            $( cut -d\   -f 2 "${PROTONATED_FILE}" | wc -l ) "protomers extracted" 1>&2
    fi
            
    echo "Coalesing and merging protomers" 1>&2
    # Pull reference Protomer (id: 0) from mid-pH protomer list
    if ! $PROTOMER_COALESE_EXE "${PROTONATED_FILES[0]}" "${PROTONATED_FILES[@]}" \
            | sed 's/\s\+/ /g' \
            > "${GENERATED_PROTOMER_SMILES}" ; then
        echo "Protomer coalese failed" | tee -a failure-reason 1>&2
        echo "Aborting!" 1>&2
        popd 1>&2
        popd 1>&2
        if [ -z "${DEBUG}" ] ; then
            mkdir -pv "${FAILED}" 1>&2
            mv -v ${IN_PROGRESS}/* ${FAILED} 1>&2
        fi
        exit -1
    fi
    if [ ! -s "${GENERATED_PROTOMER_SMILES}" ]; then 
        echo "Protomer coalese resulted in 0 protomers" | tee -a failure-reason 1>&2
        echo "Aborting!" 1>&2
        popd 1>&2
        popd 1>&2
        if [ -z "${DEBUG}" ] ; then
            mkdir -pv "${FAILED}" 1>&2
            mv -v ${IN_PROGRESS}/* "${FAILED}" 1>&2
        fi
        exit -1
    else
        echo $( cat "${GENERATED_PROTOMER_SMILES}" | wc -l ) "protomers generated for" $( cat "${INCOMMING_SMILES}" | wc -l ) "compounds" 1>&2
    fi

    if [ "${USE_3D}" == "yes" ] ; then
        echo "Using existing 3D embeddings. Not expanding stereochemistry" 1>&2
        cp -v "${GENERATED_PROTOMER_SMILES}" "${EXPANDED_PROTOMER_SMILES}" 1>&2
    else
        echo "Checking for new stereocenters and expanding" 1>&2
        if ! $PROTOMER_STEREOCENTERS_EXE "${GENERATED_PROTOMER_SMILES}" > "${EXPANDED_PROTOMER_SMILES}" ; then
            echo "Stereo expansion failed. Proeeding with original protomers" 1>&2
            cp -v "${GENERATED_PROTOMER_SMILES}" "${EXPANDED_PROTOMER_SMILES}" 1>&2
        else
            echo $( cat "${EXPANDED_PROTOMER_SMILES}" | wc -l ) "protomers after new stereo-center expansion" 1>&2
        fi
    fi

    popd 1>&2
    echo "" 1>&2
    [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP protonation post processing | " && date

else
    pushd "${IN_PROGRESS}" 1>&2
fi

[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START all building | " && date

# start rewrite here
# hello! if you are currently trying to muck around with the 3D scripts, you should know that the above code was written a generation before the below code
# the above bit of code was never a bottleneck and worked as far as I knew, so I kept it like it is. bash programming is not super readable, so good luck if you need to fix a bug :)

python ${DOCKBASE}/ligand/generate/build_ligands.py ${EXPANDED_PROTOMER_SMILES} &
buildpid=$!

function signal_build_ligands {
    kill -10 $buildpid
}
trap signal_build_ligands SIGUSR1
trap signal_build_ligands SIGINT

while [ -z "$(kill -0 $buildpid 2>&1)" ]; do
	sleep 5
done

[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP all building | " && date
