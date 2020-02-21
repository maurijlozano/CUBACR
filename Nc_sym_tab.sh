#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: Nc_sym_tab.sh [-h]

        -h Show this Help
"
echo -e '\nCalculates Nc table.\n\n'
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
    codonw modal_seqs.tmp "modal_seqs_NC_GC3s_control.tab" -enc -gc3s -nomenu -silent -nowarn &> /dev/null
 	rm modal_seqs.tmp modal_seqs.blk
 	sed -E 's/ //g' modal_seqs_NC_GC3s_control.tab | sed -E 's/_[Nn][Ee][Ww]//' > modal_seqs_NC_GC3s_control.tab2
 	mv modal_seqs_NC_GC3s_control.tab2 modal_seqs_NC_GC3s_control.tab
	cd ..
done
