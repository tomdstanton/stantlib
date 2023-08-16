#!/usr/bin/env bash

PROG_NAME=$(basename $0 ".sh")
min=10
max=100
step=10
out=$PWD

function usage {
        echo 2>&1
        echo "Usage: $PROG_NAME <file> [<file> ...] [options]"
        echo "	<file>	.pdb or .faa files to cluster"
        echo "	-out	Output directory	(default: CWD)"
        echo "	-min	Minimum %ID to cluster	(default: $min)"
        echo "	-max	Maximum %ID to cluster	(default: $max)"
        echo "	-step	Steps between min-max	(default: $step)"
        echo
        echo "	Clutering $min%ID -> $max%ID in steps of $step"
        echo
        exit 1
}

while [[ $# -gt 0 ]]; do
	case $1 in
		-out) shift; out=$1;;
		-min) shift; min=$1;;
		-max) shift; max=$1;;
		-step) shift; step=$1 ;;
		*) break ;;
  	esac
  	shift
done

# Store files in an array and check if at least one was provided
files=("$@")
[ ${#files[@]} -eq 0 ] && usage
faa=()
pdb=()

# Loop through the files and process each one
for file in "${files[@]}"; do
	if [ ! -f "$file" ]; then
		echo 2>&1 "File does not exist: $file"; exit 1
	elif [[ "$file" == *.pdb ]]; then
		pdb+=($file)
	elif [[ "$file" == *.faa ]]; then
		faa+=($file)
	else [[ "$file" != *.pdb ]]
		echo "WARNING: File doesn't have .faa/.pdb extension: $file"
	fi
done
[[ ${#pdb[@]} -eq 0 && ${#faa[@]} -eq 0 ]] && usage

echo "${#faa[@]} FAA file(s) provided"
echo "${#pdb[@]} PDB file(s) provided"
echo
echo "Output directory: $out"
echo "Minimum %ID to cluster: $min"
echo "Maximum %ID to cluster: $max"
echo "Steps between min-max: $step"
echo

# Check foldseek binary
command -v foldseek >/dev/null 2>&1 || { echo 2>&1 "foldseek binary not found"; exit 1; }
command -v mmseqs >/dev/null 2>&1 || { echo 2>&1 "mmseqs binary not found"; exit 1; }

if [ ! ${#pdb[@]} -eq 0 ]; then
	# Make temp dir
	tempdir=$(mktemp -dq ${out}/foldseek_temp_XXXXXX)
	if [ $? -ne 0 ]; then
	   echo "Can't create temp db, exiting..."
	   exit 1
	fi
	# Create db
	db="$tempdir/DB"
	echo "	FOLDSEEK: Creating db: $db"
	foldseek createdb ${files[@]} $db -v 1
	[ -f $db ] || { echo 2>&1 "Failed to create db: $db"; rm -r $tempdir; exit 1; }
	clu_prefix="$tempdir/CLU"
	tmp_="$tempdir/tmp"

	# Clustering loop
	for id in $(seq $min $step $max)
	do 
		clu=${clu_prefix}_${id}
		scale="0$(echo "scale=1; $id / 100" | bc)"
		echo "	FOLDSEEK: Clustering at $id%"
		foldseek cluster $db $clu $tmp_ --min-seq-id $scale -v 2
		[ -f ${clu}.index ] || { echo 2>&1 "Could not find ${clu}.index"; rm -r $tempdir; exit 1; }
		tsv=${out}/foldseek_cluster_${id}.tsv
		mmseqs createtsv $db $db $clu $tsv -v 2
		[ -f $tsv ] || { echo 2>&1 "Could not find $tsv"; rm -r $tempdir; exit 1; }
	done
	rm -r $tempdir
fi

if [ ! ${#faa[@]} -eq 0 ]; then
	# Make temp dir
	tempdir=$(mktemp -dq ${out}/mmseqs_temp_XXXXXX)
	if [ $? -ne 0 ]; then
	   echo "Can't create temp db, exiting..."
	   exit 1
	fi
	# Create db
	db="$tempdir/DB"
	echo "	MMSEQS: Creating db: $db"
	mmseqs createdb ${files[@]} $db -v 1
	[ -f $db ] || { echo 2>&1 "Failed to create db: $db"; rm -r $tempdir; exit 1; }
	clu_prefix="$tempdir/CLU"
	tmp_="$tempdir/tmp"

	# Clustering loop
	for id in $(seq $min $step $max)
	do 
		clu=${clu_prefix}_${id}
		scale="0$(echo "scale=1; $id / 100" | bc)"
		echo "	MMMSEQS: Clustering at $id%"
		mmseqs cluster $db $clu $tmp_ --min-seq-id $scale -v 2
		[ -f ${clu}.index ] || { echo 2>&1 "Could not find ${clu}.index"; rm -r $tempdir; exit 1; }
		tsv=${out}/mmseqs_cluster_${id}.tsv
		mmseqs createtsv $db $db $clu $tsv -v 2
		[ -f $tsv ] || { echo 2>&1 "Could not find $tsv"; rm -r $tempdir; exit 1; }
	done
	rm -r $tempdir
fi
echo "Done!"
exit 0
