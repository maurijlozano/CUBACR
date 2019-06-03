#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: tAi_Modal.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nComputes tAi Adaptaion Index of Modal frequencies.\nThis script scans all the child directories for C1 and singleton .modal_freq files and computes de adaptation index tAi.\nGenerates a tAi vs Distance plot.\n\n'

exit 1
}

while getopts yh option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help;;
    esac
done


if [ ${SKIP} != 1 ]
then
    echo -e 'Beginning calculations.\n'
else
    echo -e 'PRESS ANY KEY TO CONTINUE.\n'
    read -n 1 -s
fi


ls -d */ | while read F; do 
    echo Processing "${F}"
    cd "./${F}/"
    
    if [[ ! -f s_opts_DCBS_GNM.txt ]] ; then
        echo 'File "s_opts_DCBS_GNM.txt" is not there, skipping.'
        cd ..
        continue
    fi

    if [[ ! -f modals.counts ]] ; then
		cat $(ls */single*.modal_freq.fas) $(ls LEP/*.modal_freq.fas) $(ls HEP/*.modal_freq.fas) $(ls *LEP*.modal_freq.fas) $(ls *HEP*.modal_freq.fas) $(ls PHE/*PHE*.modal_freq.fas) > modal.counts
		compile_codon_counts  < modal.counts > modal.count
		sed -e 's/  / /g' modal.count > modal.counts
		cut -f 1 modal.counts > modal.count
		echo -e $(ls */single*.modal_freq.fas) $(ls LEP/*.modal_freq.fas) $(ls HEP/*.modal_freq.fas) $(ls *LEP*.modal_freq.fas) $(ls *HEP*.modal_freq.fas) $(ls PHE/*PHE*.modal_freq.fas) | tr ' ' '\n' | sed -e 's/[^\/]*\///' -e 's/.modal_freq.fas//' -e 's/[_ ][RrNn][Ee][WwFf]$//' -e 's/__/_/' > modal.rname
		paste -d' ' modal.rname modal.count | sed -e 's/\^\n//' > modal.counts
		rm *.count modal.rname
	fi

    cd ..

    Rscript --vanilla Tai_modal.R "./${F}"

done
