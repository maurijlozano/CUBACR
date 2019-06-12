#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: CRdelta_heatmap.sh [-h]

        -h Show this Help
"
echo -e '\nCalculates RSCU difference between HEP-CR, HEP and PHE.\n\n'
exit 1
}

while getopts h option; do
    case "${option}" in
       h) show_help;;
    esac
done



ls -d */ | while read F; do 
    echo Processing "${F}"
	FE=$(echo "${F}" | sed -e 's/[_ ][_ ][_ ]/_/' | sed -e 's/[_ ][_ ]/_/')
    FS=$(echo "${FE}" | sed -e 's/\(^[A-Za-z0-9]*\)[_ ]\(.\)[A-Za-z0-9]*[_ ]\([A-Za-z0-9]*\).*/\1_\2.\3/')
    FN=$(echo "${FS}" | sed -e 's/\///')

    cd "./${F}/"
    
    if [[ ! -f fmdata.txt ]] ; then
		echo 'Please run calculate_dist_tree_heatmap.sh first.'
    fi
  
    cd ..
    Rscript --vanilla heatmap.r "./${F}" "${FN}"
    cd "./${F}/"
    if [[ -d LEP ]] ; then
    	cd ..
    	Rscript --vanilla heatmap2.r "./${F}" "${FN}"
    else
    	cd ..	
	fi
done
