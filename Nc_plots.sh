#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: Nc_plots.sh [-h]

        -h Show this Help
"
echo -e '\nNc plots for CR and VR regions.\n\n'
exit 1
}

while getopts h option; do
    case "${option}" in
       h) show_help;;
    esac
done



ls -d */ | while read D; do 
    echo Processing "${D}"
	cd "./${D}/"
    
    cat $(ls CC/*.modal_freq.fa* PHE/*.modal_freq.fa* HEP/*.modal_freq.fa* LEP/*.modal_freq.fa* *.modal_freq.fa*) > modal_seqs.tmp
    sed -E 's/>[^\/]*\/[^\/]*\//>/' modal_seqs.tmp > modal_seqs.tmp2
    mv modal_seqs.tmp2 modal_seqs.tmp
    codonw modal_seqs.tmp "modal_seqs_NC_GC3s.tab" -enc -gc3s -nomenu -silent -nowarn &> /dev/null
 	rm modal_seqs.tmp modal_seqs.blk
 	sed -E 's/ //g' modal_seqs_NC_GC3s.tab > modal_seqs_NC_GC3s.tab2
 	mv modal_seqs_NC_GC3s.tab2 modal_seqs_NC_GC3s.tab
 	
 	ls CC/*.fa* PHE/PHE*.fa* HEP/HEP.fa* LEP/LEP.fa* *.fa* | grep -v 'modal_freq' | grep -v 'concat' | while read F; do 
		FN=$(echo "${F}" | sed -E 's/[^\/]*\///' | sed -E 's/[^_]*_([^_]*)_([^_]*).fas?$/\1-\2/' | sed -E 's/([^_]*).*/\1/' | sed -E 's/.fas?$//')
		codonw "${F}" "${FN}.tab" "${FN}.blk" -enc -gc3s -nomenu -silent -nowarn &> /dev/null
		sed -E 's/([^\t]*)/'"${FN}"'/'  "${FN}.tab" | tail -n +2 | grep '*' -v >> NC_GC3s.tmp
		rm "${FN}.blk" "${FN}.tab"
 	done
 	mv NC_GC3s.tmp NC_GC3s.tab
    
    cd ..
    Rscript --vanilla Nc_plots.R "./${D}"
  
done
