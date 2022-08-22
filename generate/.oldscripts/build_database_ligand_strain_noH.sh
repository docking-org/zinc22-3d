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

[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START 3D embedding | " && date
echo "Bulk generating 3D conformations all protomers in ${CONFORMATIONS}" 1>&2
mkdir -pv "${CONFORMATIONS}" 1>&2
if [ "${USE_3D}" == "yes" ] ; then
    EMBED_SOURCE="${SOURCE_FILE}"
    RUN_EMBEDDING="$EXPAND_3D_EXE "${EMBED_SOURCE}" ${CONFORMATIONS}/"
else
    EMBED_SOURCE="${EXPANDED_PROTOMER_SMILES}"
    RUN_EMBEDDING="$EMBED_PROTOMERS_3D_EXE  "${EMBED_SOURCE}" -o ${CONFORMATIONS}/"
fi
if ! $RUN_EMBEDDING 1>&2 ; then
    echo "Protomer conformation generation failed" | tee -a failure-reason 1>&2
    echo "Abortingy" 1>&2
    popd 1>&2
    if [ -z "${DEBUG}" ] ; then
        mkdir -pv "${FAILED}" 1>&2
        mv -v ${IN_PROGRESS}/* ${FAILED} 1>&2
    fi
    exit -1
else
    echo $( ls ${CONFORMATIONS} | wc -l ) "3D conformations generated for" $( cat ${INCOMMING_SMILES} | wc -l ) "compounds" 1>&2
    echo "" 1>&2
fi
[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP 3D embedding | " && date
       
[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START all building | " && date

LINE_NUM=0
cat "${INCOMMING_SMILES}" |
while read -r LINE
do
    LINE_NUM=`expr $LINE_NUM + 1`
    if [ -z "$LINE" ]; then
        continue
    fi

    TOKENS=( ${LINE//\s\+/} )
    SMILES="${TOKENS[0]}"
    NAME="${TOKENS[1]}"
    if [ -z "$SMILES" ] ; then
        echo "Missing information for line #${COMPOUND_NUM}" 1>&2
        echo "\t${NAME} ${SMILES}" 1>&2
        echo "Skipping" 1>&2
        continue
    elif [ -z "${NAME}" ] ; then
        echo "Missing name for line #${COMPOUND_NUM}" 1>&2
        NAME="${DBNAME}-${LINE_NUM}"
        _NAME_WAS_MISSING="yes"
    fi
        
    echo "Building $NAME" 1>&2

    MOL_DIR="${BUILD_DIR}/${NAME}"
    mkdir -pv "${MOL_DIR}" 1>&2
    pushd "${MOL_DIR}" 1>&2

    echo "${LINE}" > "${MOL_DIR}/input.ism"

    # Start output files
    NUMBERED_SMILES_FILE="${MOL_DIR}/${NAME}-numbered.ism"
    SMILES_FILE="${MOL_DIR}/${NAME}.ism"
    SOLV_FILE="${MOL_DIR}/${NAME}.solv"

    echo "Extracting previously generated protomers and correcting pH mod types" 1>&2
   
    grep -n "\s${NAME}\s" "${EXPANDED_PROTOMER_SMILES}" \
        | awk '{
        if (NR > 1 && $3 == "0") {
            print $1, $2, "1"
        } else {
            print $0
        }
    }' | sed -e 's/:/	/' -e 's/\s\+/	/g' > "${NUMBERED_SMILES_FILE}"
    cut -f2- "${NUMBERED_SMILES_FILE}" > "${SMILES_FILE}"
    if [ ! -s "${SMILES_FILE}" ]; then 
        echo "Protomer extracton resulted in 0 protomers" | tee -a failure-reason 1>&2
        echo "Marking $NAME as failed and skipping" 1>&2
        popd 1>&2
        mkdir -pv ${FAILED} 1>&2
        if [ -z "${DEBUG}" ] ; then
            mv -v ${MOL_DIR} ${FAILED}/${NAME} 1>&2
        else
            cp -rv ${MOL_DIR} ${FAILED}/${NAME} 1>&2
        fi
        continue
    else
        echo $( cat "${SMILES_FILE}" | wc -l ) "protomers extracted for" ${NAME} 1>&2
        echo "" 1>&2
    fi

    PROT_NUM=0
    cat "${NUMBERED_SMILES_FILE}" |
    while read -r PROTOMER_LINE
    do
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START protomer ${PROT_NUM} | " && date
        PROT_TOKENS=(${PROTOMER_LINE// / })          # Allow split on colon as well
        PROT_IDX=${PROT_TOKENS[0]}            # "Global" protomer index
        PROT_SMILES=${PROT_TOKENS[1]}            # SMILES
        COMPOUND_NAME=${PROT_TOKENS[2]}            # Compound Name
        PROT_TYPE=${PROT_TOKENS[3]}            # ph_mod_fk
        
        PROT_NAME="${PROT_NUM}"
        PROT_3D="${PROT_NUM}.mol2"
        PROT_3D_UNCHARGED="${PROT_3D}.original"
        PROT_3D_CHARGED='output.mol2'
        mkdir -pv $PROT_NUM 1>&2
        pushd ${PROT_NUM} 1>&2

        echo "Protomer ${PROT_NAME} (index: ${PROT_IDX})" 1>&2

        # Is this conformation actually for our compound? 
        # (Check that the compound name is exactly on the 2nd line)
	    PROT_EMBEDDED_NAME="$( grep -A 1 '@<TRIPOS>MOLECULE' "${CONFORMATIONS}/${PROT_IDX}" | tail -n+2 | sed 's/\s\+//g' )"
        if [ "${PROT_EMBEDDED_NAME}" == "${COMPOUND_NAME}" ] ; then
            echo "Found valid previously generated 3D confromation in ${CONFORMATIONS}/${PROT_IDX}" 1>&2
            mv -v "${CONFORMATIONS}/${PROT_IDX}" "${PROT_3D}" 1>&2
        elif [ "${_NAME_WAS_MISSING}" == "yes" -o "${USE_3D}" == "yes" ] ; then
            echo "No name was provided for compound #$LINE_NUM. Using ${CONFORMATIONS}/${PROT_IDX} by default for 3D embedding" 1>&2
            mv -v "${CONFORMATIONS}/${PROT_IDX}" "${PROT_3D}" 1>&2
       else
            echo "No valid previously generated conformations found!" 1>&2
            echo "Embedding SMILES in 3D and adding Hydrogens for $PROT_NAME" 1>&2
            if ! echo "${PROT_SMILES} ${COMPOUND_NAME}" | $EMBED_PROTOMERS_3D_EXE -s -o "${PROT_3D}" - 1>&2 ; then
                echo "Failed to generate 3D embedding in both batch and single-runs" | tee -a failure-reason 1>&2
                echo "Skipping $NAME $PROT_NAME" 1>&2
                popd 1>&2
                mv -v $PROT_NUM $PROT_NUM.failed 1>&2 
                break 
            elif [ ! -e "${PROT_3D}" -o ! -s "${PROT_3D}" ] ; then
                echo "Failed to generate 3D embedding in both batch and single-runs" | tee -a failure-reason 1>&2
                echo "Skipping $NAME $PROT_NAME" 1>&2
                popd 1>&2
                mv -v $PROT_NUM $PROT_NUM.failed 1>&2 
                break 
            fi
        fi

        echo "Preparing input files" 1>&2
        if ! $PREPARE_NAME_EXE \
        -n "${COMPOUND_NAME}" \
        -s "${PROT_SMILES}" \
        -t "${PROT_NUM}" \
        "${PROT_3D}" 1>&2 ; then
            echo "Protomer database file preparation failed" | tee -a failure-reason 1>&2
            echo "Skipping $NAME $PROT_NAME" 1>&2
            popd 1>&2
            mv -v $PROT_NUM $PROT_NUM.failed 1>&2 
            break 
        fi

        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START solvation ${PROT_NUM} | " && date
        echo "" 1>&2
        echo "Starting the preparation of the solvation calculations (AMSOL7.1)" 1>&2
        echo "  (SMILES: $PROT_SMILES)" 1>&2
        echo "" 1>&2
        ln -svfn "${PROT_3D}" "${COMPOUND_NAME}.mol2" 1>&2
        if ! $SOLVATION_EXE "${COMPOUND_NAME}.mol2" 1>&2 ; then
            echo "Solvation calculation failed" | tee -a failure-reason 1>&2
            echo "Skipping $NAME $PROT_NAME" 1>&2
            popd 1>&2
            if [ -z "${DEBUG}" ] ; then
                mv -v $PROT_NUM $PROT_NUM.failed 1>&2 
            else
                cp -rv $PROT_NUM $PROT_NUM.failed 1>&2
            fi
            break 1
        else
            mv -v "${PROT_3D}" "${PROT_3D_UNCHARGED}" 1>&2
            cp -v "${PROT_3D_CHARGED}" "${PROT_3D}" 1>&2
        fi
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP solvation ${PROT_NUM} | " && date

        if [ "${CREATE_CONFORMATIONS}" == "yes" ] ; then        # Generate and extract conformations
            [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START conformations ${PROT_NUM} | " && date
            if [ "${HYDROGEN_DEPENDENT_CONFORMATION_LIMITS}" != "no" ] ; then
                H_COUNT="$( ${COUNT_ROTATABLE_HYDROGENS_EXE} "${PROT_3D_CHARGED}" )"
            else
                H_COUNT=""
            fi
 
            if ! $GENERATE_CONFORMATIONS_EXE "${PROT_3D_CHARGED}" $H_COUNT 1>&2 ; then
                echo "Conformer generation failed" | tee -a failure-reason 1>&2
                echo "Skipping $NAME $PROT_NAME" 1>&2
                popd 1>&2
                if [ -z "${DEBUG}" ] ; then
                    mv -v $PROT_NUM $PROT_NUM.failed 1>&2
                else
                    cp -rv $PROT_NUM $PROT_NUM.failed 1>&2
                fi
                break 1
            fi
            [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP conformations ${PROT_NUM} | " && date
        else
            echo "No flexibility being generated!" 1>&2
            cp -v "${PROT_3D}" output.1.db2in.mol2 1>&2
        fi

        DBIN_FILES=( $( find . -maxdepth 1 -name '*.db2in.mol2' ) )
        if [[ ${#DBIN_FILES[@]} -lt 1 ]] ; then 
            echo "Conformer generation produced nothing" | tee -a failure-reason 1>&2
            echo "Skipping ${NAME} ${PROT_NAME}" 1>&2
            popd 1>&2
            if [ -z "${DEBUG}" ] ; then
                mv -v "${PROT_NUM}" "${PROT_NUM}.failed" 1>&2
            else
                cp -rv "${PROT_NUM}" "${PROT_NUM}.failed" 1>&2
            fi
            break 1
        fi

        # Build DB2 files from extracted conformations
        if [ "${BUILD_DB2_FILES}" == "yes" ] ; then
            DBIDX=0
            for DBIN in "${DBIN_FILES[@]}" ; do
                [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START db2 ${PROT_NUM} ${DBIDX} | " && date
                if ! $UCSF_STRAIN_EXE "${DBIN}" 1>&2 ; then
                    echo "UCSF strain failed" | tee -a failure-reason 1>&2
                    echo "Skipping $NAME $PROT_NAME" 1>&2
                    popd 1>&2
                    mv -v $PROT_NUM $PROT_NUM.failed  1>&2
                    break 2
                fi
                prefix=`basename ${DBIN} .mol2`
                #if ! $BUILD_DB2_EXE --mol2 "${DBIN}" --db "${DBIN}.db2.gz" 1>&2 ; then
                if ! $BUILD_DB2_EXE --mol2 "${prefix}.strain.mol2" --db "${DBIN}.db2.gz" 1>&2 ; then
                    echo "DB2 File generation failed" | tee -a failure-reason 1>&2
                    echo "Skipping $NAME $PROT_NAME" 1>&2
                    popd 1>&2
                    mv -v $PROT_NUM $PROT_NUM.failed  1>&2
                    break 2
                fi
                [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP db2 ${PROT_NUM} ${DBIDX} | " && date
                DBIDX=$(( $DBIDX + 1 ))
            done
            if ! ( zcat output.*.db2.gz > "${PROT_NAME}.db2" && gzip -9 "${PROT_NAME}.db2" ); then
                echo "DB2 Files not found as expected" | tee -a failure-reason 1>&2
                echo "Skipping $NAME $PROT_NAME" 1>&2
                popd 1>&2
                mv -v $PROT_NUM $PROT_NUM.failed 1>&2
                break 1
            fi
        fi

        if [ "${BUILD_DB_FILES}" == "yes" ] ; then
            DBIDX=0
            for DBIN in "${DBIN_FILES[@]}" ; do
                [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START db ${PROT_NUM} ${DBIDX} | " && date
                ln -sv "${DBIN}" db.mol2 1>&2
                ln -svfn output.solv mol.solv 1>&2
                if ! $BUILD_DB_EXE 1>&2 ; then  # Parameters from inheir
                    echo "DB File generation failed" | tee -a failure-reason 1>&2
                    echo "Skipping $NAME $PROT_NAME" 1>&2
                    popd 1>&2
                    mv -v $PROT_NUM $PROT_NUM.failed 1>&2
                    break 2
                fi

                mv -v db.db "${DBIN}.db" 1>&2
                rm -v db.mol2 1>&2
                gzip -9 "${DBIN}.db" 1>&2
                [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP db ${PROT_NUM} ${DBIDX} | " && date
                DBIDX=$(( $DBIDX + 1 ))
            done
            if ! ( zcat output*.db.gz > "${PROT_NAME}.db" && gzip -9 "${PROT_NAME}.db" ); then
                echo "DB Files not found as expected!" 1>&2
                echo "Proceeding" 1>&2
            fi
        fi
        
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START finalize ${PROT_NUM} | " && date
        gzip -9 "${PROT_3D}"
        cp output.solv ${PROT_NAME}.solv

        head -n1 output.solv |\
        	cut -f2- -d' ' |\
	        sed 's/^\s\+//' |\
        	sed "s/\s\+/    /g" >> ${SOLV_FILE}

        popd 1>&2
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP finalize ${PROT_NUM} | " && date
        [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP protomer ${PROT_NUM} | " && date
        PROT_NUM=`expr $PROT_NUM + 1`
    done

    echo "" 1>&2
    echo "Finished preparing $NAME" 1>&2
    echo "Recording results" 1>&2

    popd 1>&2
    if [ $( find ${MOL_DIR} -maxdepth 1 -name '*.failed' -print -quit ) ]; then 
        echo "Found failed protomers" | tee -a ${MOL_DIR}/failure-reason 1>&2
        echo "Marking $NAME as failed" 1>&2
        mkdir -pv ${FAILED} 1>&2
        mv -v ${MOL_DIR} ${FAILED}/${NAME} 1>&2
    elif [ "${SINGLE_DATABASE}" == "yes" ] ; then
        mkdir -pv "${FINISHED}" 1>&2
        if [ "${SAVE_PROTOMER_TABLE}" == "yes" ] ; then
            [ ! -e "${FINISHED}/${DBNAME}.xls" ] && \
                echo "smiles name type num_atoms total_charge polar_solv surface apolar_solv total_solv" \
                    | sed "s/\s\+/\t/g" > "${FINISHED}/${DBNAME}.xls"
            paste ${SMILES_FILE} ${SOLV_FILE} | sed 's/\s\+/	/g' >> "${FINISHED}/${DBNAME}.xls"
        fi
        echo "Appending to ${FINISHED}/${DBNAME}.*" 1>&2
        for PROT_DIR in $( find ${MOL_DIR} -mindepth 1 -maxdepth 1 -type d ); do
            PROT_NUM=$( basename $PROT_DIR )
            echo "$PROT_NUM: $PROT_DIR.*" 1>&2
            [ "${COPY_MOL2_FILES}" == "yes" ] && \
                gzip -dc "${PROT_DIR}/${PROT_NUM}.mol2.gz" >> "${FINISHED}/${DBNAME}.mol2"
            [ "${BUILD_DB2_FILES}" == "yes" ] && \
                gzip -dc "${PROT_DIR}/${PROT_NUM}.db2.gz" >> "${FINISHED}/${DBNAME}.db2"
            [ "${BUILD_DB_FILES}" == "yes" ] && \
                gzip -dc "${PROT_DIR}/${PROT_NUM}.db.gz" >> "${FINISHED}/${DBNAME}.db"
            [ "${COPY_SOLV_FILES}" == "yes" ] && \
                cat "${PROT_DIR}/${PROT_NUM}.solv" >> "${FINISHED}/${DBNAME}.solv"
            PROT_NUM=$(( $PROT_NUM + 1 ))
        done
        echo "Removing working files in ${MOL_DIR}" 1>&2
        if [ -z "${DEBUG}" ] ; then
            rm -rf ${MOL_DIR} 1>&2
        fi
    else
        echo "smiles name type num_atoms total_charge polar_solv surface apolar_solv total_solv" > "${MOL_DIR}/protomers.xls"
        paste ${SMILES_FILE} ${SOLV_FILE} >> "${MOL_DIR}/protomers.xls"
        sed -i "s/\s\+/\t/g" "${MOL_DIR}/protomers.xls"
        mkdir -pv ${FINISHED}/${NAME} 1>&2
        find ${MOL_DIR} -mindepth 1 -maxdepth 1 -type f -exec cp -v '{}' ${FINISHED}/${NAME}/ \; 1>&2
        for PROT_DIR in $( find ${MOL_DIR} -mindepth 1 -maxdepth 1 -type d ); do
            PROT_NUM=$( basename $PROT_DIR )
            [ "${COPY_MOL2_FILES}" == "yes" ] && \
                mv -v "${PROT_DIR}/${PROT_NUM}.mol2.gz" "${FINISHED}/${NAME}/" 1>&2
            [ "${BUILD_DB2_FILES}" == "yes" ] && \
                mv -v "${PROT_DIR}/${PROT_NUM}.db2.gz" "${FINISHED}/${NAME}/" 1>&2
            [ "${BUILD_DB_FILES}" == "yes" ] && \
                mv -v "${PROT_DIR}/${PROT_NUM}.db.gz" "${FINISHED}/${NAME}/" 1>&2
            [ "${COPY_SOLV_FILES}" == "yes" ] && \
                mv -v "${PROT_DIR}/${PROT_NUM}.solv" "${FINISHED}/${NAME}/" 1>&2
        done
        echo "Removing working files in ${MOL_DIR}" 1>&2
        if [ -z "${DEBUG}" ] ; then
            rm -rf ${MOL_DIR} 1>&2
        fi
    fi
done
[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP all building | " && date

popd 1>&2
[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: START storing | " && date
if [ "${SINGLE_DATABASE}" == "yes" ] ; then
    echo "Compressing combined databse files" 1>&2
    [ -e "${FINISHED}/${DBNAME}.mol2" ] && gzip -9 "${FINISHED}/${DBNAME}.mol2" && echo "${FINISHED}/${DBNAME}.mol2.gz" 1>&2
    [ -e "${FINISHED}/${DBNAME}.db2" ] && gzip -9 "${FINISHED}/${DBNAME}.db2" && echo "${FINISHED}/${DBNAME}.db2.gz" 1>&2
    [ -e "${FINISHED}/${DBNAME}.db" ] && gzip -9 "${FINISHED}/${DBNAME}.db" && echo "${FINISHED}/${DBNAME}.db.gz" 1>&2
    [ -e "${FINISHED}/${DBNAME}.solv" ] && gzip -9 "${FINISHED}/${DBNAME}.solv" && echo "${FINISHED}/${DBNAME}.solv.gz" 1>&2
elif [ "${SAVE_PROTOMER_TABLE}" == "yes" ] ; then
    [ ! -e "${FINISHED}/${DBNAME}.xls" ] && \
        echo "smiles name type num_atoms total_charge polar_solv surface apolar_solv total_solv" \
            | sed 's/\s\+/	/g' > "${FINISHED}/${DBNAME}.xls"
    for PROT_XLS in $( find "${FINISHED}" -mindepth 1 -maxdepth 1 -name "protomers.xls" ); do
        tail -n+2 "${PROT_XLS}" >> "${FINISHED}/${DBNAME}.xls"
    done
fi
if [ -z "$( find "${FINISHED}" -mindepth 1 -maxdepth 1 -print -quit )" ] ; then
    echo 1>&2
    echo "=======================================================" 1>&2
    echo "ERROR: No protomers successfully built!" 1>&2
    echo "=======================================================" 1>&2
    echo 1>&2
elif [ -d "${STORE_PROTOMERS}" ] ; then
    echo "Moving results into ${STORE_PROTOMERS}" 1>&2
    mv -v ${FINISHED}/*/ "${STORE_PROTOMERS}" 1>&2
elif [ ! -z "${STORE_PROTOMERS}" ] ; then
    echo "Running external protomer storage: ${STORE_PROTOMERS}" 1>&2
    if ! $STORE_PROTOMERS "${FINISHED}" 1>&2 ; then
        echo "External protomer storage failed" | tee -a failure-reason 1>&2
        echo "Some finished compounds may not have been loaded!" 1>&2
        echo "Marking all remaining finished protomers as failed" 1>&2
        mkdir -pv "${UNLOADED}" 1>&2
        mv -v "${FINISHED}"/* "${UNLOADED}/" 1>&2
    fi
else
    echo 1>&2
    echo 1>&2
    echo "=======================================================" 1>&2
    echo "WARNING: STORE_PROTOMERS not executable or a directory!" 1>&2
    echo "All results left in place (${FINISHED})" 1>&2
    echo "=======================================================" 1>&2
    echo 1>&2
fi
[ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP storing | " && date

echo "Finalizing..." 1>&2
if [ -z "${DEBUG}" ] ; then
   [ ! -z "${DEBUG}" ] &&  echo -n "BENCHMARK: START cleanup | " && date
    rm -rfv ${CONFORMATIONS} 1>&2
    if [ $( find ${BUILD_DIR} -mindepth 1 -maxdepth 2 -name '*' -print -quit ) ]; then
        mkdir -pv "${ARCHIVE}" 1>&2
        mv -v ${BUILD_DIR}/* ${ARCHIVE} || echo "None moved" 1>&2
    fi
    if [ -d "${FAILED}" ] ; then
        mv -v "${PROTONATING_DIR}" "${ARCHIVE}" 1>&2
	mv -v "${CLEANED_SMILES}" "${ARCHIVE}" 1>&2
    else
        rm -rfv "${PROTONATING_DIR}"
	rm -fv "${CLEANED_SMILES}"
    fi
    rmdir -v ${BUILD_DIR} 1>&2
    rmdir -v ${IN_PROGRESS} 1>&2
    [ ! -z "${DEBUG}" ] && echo -n "BENCHMARK: STOP cleanup | " && date
fi
