#!/bin/bash

SINGLE="no"
SOURCE_IS_ROOT="no"
SOURCES=()
STDIN_SOURCES="no"
DEST_DIR=""
FORMATS=( db2.gz db.gz )

while [ "$#" -gt 0 ] ; do
	case "${1}" in
		-s|--single)
			SINGLE="yes"
			shift 1
			;;
		-d|--dirs)
			SOURCE_IS_ROOT="yes"
			shift 1
			;;
		-o|--output)
			DEST="${2}"
			shift 2
			;;
		-f|--formats)
			FORMATS=( ${2} )
			shift 2
			;;
		-)
			STDIN_SOURCES="yes"
			shift 1
			;;
		*)
			SOURCES=( "${SOURCES[@]}" "${1}" )
			shift 1
			;;
	esac
done

if [ "${#SOURCES}" -eq 0 -a "${STDIN_SOURCES}" == "no" ] ; then
	SOURCES=( $( pwd ) )
fi

if [ -z "${DEST}" ] ; then
	DEST=$( pwd )
else
	mkdir -pv "${DEST}" 1>&2
fi	

for SOURCE in  "${SOURCES[@]}" ; do
	if [ "${SOURCE_IS_ROOT}" == "yes" ] ; then
		INNER_SOURCES=( $( ls -d ${SOURCE}/*/ ) )
	else
		INNER_SOURCES=( "${SOURCE}" )
	fi
	COUNTER=0
	for INNER in "${INNER_SOURCES[@]}" ; do
		NAME=$( basename "${INNER}" )
		INNER_DEST="${DEST}/${NAME}"
		COUNTER=$(( $COUNTER + 1 ))
		for FORMAT in "${FORMATS[@]}" ; do
			if [[ "${FORMAT}" == *.gz ]] ; then
				gzip -dc "${INNER}"/*."${FORMAT}" | gzip -9 > "${INNER_DEST}.${FORMAT}"
			else
				cat "${INNER}"/*."${FORMAT}" > "${INNER_DEST}.${FORMAT}"
			fi
		done
	done
	echo "${SOURCE}": "${COUNTER}" 1>&2
done

if [ "${STDIN_SOURCES}" == "yes" ] ; then
	while read SOURCE ; do
		if [ "${SOURCE_IS_ROOT}" == "yes" ] ; then
			INNER_SOURCES=( $( ls ${SOURCE}/*/ ) )
		else
			INNER_SOURCES=( "${SOURCE}" )
		fi
		
		COUNTER=0
		for INNER in "${INNER_SOURCES[@]}" ; do
			COUNTER=$(( $COUNTER + 1 ))
			NAME=$( basename "${INNER}" )
			INNER_DEST="${DEST}/${NAME}"
			for FORMAT in "${FORMATS[@]}" ; do
				if [ "${FORMAT}" == *.db2.gz ] ; then
					gzip -dc "${INNER}"/*."${FORMAT}" | gzip -9 > "${INNER_DEST}.${FORMAT}"
				else
					cat "${INNER}"/*."${FORMAT}" > "${INNER_DEST}.${FORMAT}"
				fi
			done
		done
		echo "${SOURCE}": "${COUNTER}" 1>&2
	done
fi

if [ "${SINGLE}" == "yes" ] ; then
	for FORMAT in "${FORMATS[@]}" ; do
		if [ "${FORMAT}" == *.db2.gz ] ; then
			gzip -dc "${DEST}"/*."${FORMAT}" | gzip -9 > "${DEST}.${FORMAT}"
		else
			cat "${DEST}"/*."${FORMAT}" > "${DEST}.${FORMAT}"
		fi
		echo "Wrote ${DEST}.${FORMAT}" 1>&2
	done
	if [ $( readlink -f "${DEST}" ) != $( readlink -f $( pwd ) ) ] ; then
		rm -rf "${DEST}"/
	fi
fi
