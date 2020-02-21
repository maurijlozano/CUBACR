#!/bin/bash
clear
VI=0.5
show_help()
{
echo "
        Usage: modelParam2tab.sh [-h]

        -h Show this Help

"
echo -e '\nGenerates a table of M0 model parameters for every protein.\n\n'

exit 1
}

while getopts ":h" option; do
    case "${option}" in
        h) show_help
            exit 1
            ;;
        \?) echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

WD=$(pwd)
file="${WD}/modelParam.txt"
echo "SP HEP/LEP ProtID Ns Ls TreeLength Omega kappa" > "${file}"
#
ls -d */*/*/ | grep 'CC' -v | grep 'PHE' -v | while read D; do
	SP=$(echo $D | sed -E 's/^([^/]*)\/([^/]*)\/([^/]*)\/$/\1/')
	HL=$(echo $D | sed -E 's/^([^/]*)\/([^/]*)\/([^/]*)\/$/\2/')
	PI=$(echo $D | sed -E 's/^([^/]*)\/([^/]*)\/([^/]*)\/$/\3/')
	echo -e "Generating table for $SP, set $HL, protein ${PI}."
	cd "${D}"
	#results.mlc
	if [ ! -f "results.mlc" ]; then
		echo -e "Skipping, results.mlc required..."
		echo "${SP} ${HL} ${PI} Skipped no results.mlc file" >> "${file}"
		continue
	else
		echo -n "${SP} ${HL} ${PI} " >> "${file}"
		line=$(grep 'ns =' results.mlc)
		echo -n $line | sed -E 's/^ns = ([0-9]{1,3}) ls = ([0-9]{1,5}$)/\1 \2 /' >> "${file}"
		line=$(grep 'tree length = ' results.mlc)
		echo -n "${line} "| sed -E 's/^tree length = //' >> "${file}"
		line=$(grep "omega (dN/dS) = " results.mlc)
		echo -n "${line} " | sed -E 's/^omega \(dN\/dS\) =  ?//' >> "${file}"
		line=$(grep "kappa (ts/tv) = " results.mlc)
		echo "${line} " | sed -E 's/^kappa \(ts\/tv\) =  ?//' >> "${file}"
	fi
	cd "$WD"
done



