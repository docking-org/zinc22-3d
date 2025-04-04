# syntax=docker/dockerfile:1
FROM ubuntu:22.04

# Install needed software
RUN apt-get update && \
    apt-get install -y --no-install-recommends\
    build-essential\
    screen\
    csh\
    sudo\
    zlib1g-dev\
    openjdk-8-jdk\
    wget\
    ncurses-bin\
    time &&\
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# Copy the required licenses
COPY --from=jchem_license_folder license.cxl /licenses/jchem-license.cxl
COPY --from=openeye_license_folder oe_license.txt /licenses/oe_license.txt

# Copy the required folders from the "base" build context (containing the actual software - /mnt/nfs/soft/dock/versions/dock38/pipeline_3D_ligands)
COPY --from=base /soft/jchem-latest /soft/jchem-latest
COPY --from=base /soft/dock-latest /soft/dock-latest
COPY --from=base /soft/corina-latest /soft/corina-latest
COPY --from=base /soft/openbabel-latest /soft/openbabel-latest
COPY --from=base /soft/extralibs-latest /soft/extralibs-latest
COPY --from=base /soft/pyenv-latest /soft/pyenv-latest

# Untar all the required software
RUN mkdir -p /jchem && tar -xzf /soft/jchem-latest -C /jchem --strip-components=1
RUN mkdir -p /dock && tar -xzf /soft/dock-latest -C /dock --strip-components=1
RUN mkdir -p /corina && tar -xzf /soft/corina-latest -C /corina --strip-components=1
RUN mkdir -p /openbabel && tar -xzf /soft/openbabel-latest -C /openbabel --strip-components=1
RUN mkdir -p /extralibs && tar -xzf /soft/extralibs-latest -C /extralibs --strip-components=1
RUN mkdir -p /pyenv && tar -xzf /soft/pyenv-latest -C /pyenv --strip-components=1

# Set the base paths
ENV CHEMAXONBASE="/jchem"
ENV DOCKBASE="/dock"
ENV CORINABASE="/corina"
ENV OBABELBASE="/openbabel"
ENV EXTRALIBSBASE="/extralibs"

# Set Licenses
ENV OE_LICENSE="/licenses/oe_license.txt"
ENV CHEMAXON_LICENSE_URL="/licenses/jchem-license.cxl"

# Set building pipeline default values
ENV OMEGA_RMSD="0.5"
ENV OMEGA_FF="MMFF94Smod"
ENV OMEGA_MAX_CONFS="600"
ENV pH_LEVEL="7.4"
ENV OMEGA_TORLIB="Original"
ENV OMEGA_ENERGY_WINDOW="12"
ENV CORINA_MAX_CONFS="1"


# Sets other obabel paths
# TODO: Don't pass this as an arg - extract it from openbabel-latest
ARG OPENBABEL_VERSION_NUMBER
ENV BABEL_LIBDIR="${OBABELBASE}/lib/openbabel/${OPENBABEL_VERSION_NUMBER}"
ENV BABEL_DATADIR="${OBABELBASE}/share/openbabel/${OPENBABEL_VERSION_NUMBER}"

# Points to all the correct software that is needed
ENV MOLCONVERTEXE="${CHEMAXONBASE}/bin/molconvert"
ENV CXCALCEXE="${CHEMAXONBASE}/bin/cxcalc"
ENV OBABELEXE="${OBABELBASE}/bin/obabel"
ENV CORENAEXE="${CORINABASE}/corina"
ENV SPLITONEXE="${DOCKBASE}/ligand/3D/spliton_docker.py"
ENV EMBED_PROTOMERS_3D_EXE="${DOCKBASE}/ligand/3D/embed3d_corina.sh"
ENV TAUOMERIZE_PROTONATE_EXE="${DOCKBASE}/ligand/protonate/tautprot_cxcalc.sh"
ENV PROTOMER_STEREOCENTERS_EXE="${DOCKBASE}/ligand/protonate/expand-new-stereocenters.py3.py"
ENV AMSOLEXE="${DOCKBASE}/ligand/amsol/amsol7.1"

# Set Java home
ENV JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64"

# Adds required software to paths
ENV LD_LIBRARY_PATH="${OBABELBASE}/lib:${EXTRALIBSBASE}"
ENV PATH="${CORINABASE}:${OBABELBASE}/bin:${CHEMAXONBASE}/bin:${DOCKBASE}/bin:${PATH}"

# Link sh to bash (default is dash) b/c all shoichet things are expecting this
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Make sure the python environment is activated when we run the container
RUN echo '. /pyenv/bin/activate' >> ~/.bashrc
