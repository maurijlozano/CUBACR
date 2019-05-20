#!/bin/bash
clear

show_help()
{
echo "
        Usage: calculate_dist_tree_heatmap.sh [-h]]

        -h Show this Help

"

echo -e '\nGenerates a NJ distances tree and heatmap of CUF for Ci, Singletons, and expression sets.\n\n'

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


ls -d */ | while read D; do
	cd "$D"
	TN=$(echo "${D}" | sed 's/\/$//')
	cat $(ls *.modal_freq | sort -V) $(ls */*.modal_freq | sort -V) | sed -e 's/,/ /g' -e 's/\t/ /g' -e 's/|/ /g' > CC_modal_freqs.txt
	ls *.modal_freq | sort -V | sed -e 's/.modal_freq$//' -e 's/[^_]*_//' > CCindex.txt
	ls */*.modal_freq | sort -V | sed -e 's/.modal_freq$//' -e 's/\([^\/]*\)\///' -e 's/\([^_]*\)_.*/\1/' >> CCindex.txt
	paste -d' ' CCindex.txt CC_modal_freqs.txt > fmdata.txt
	rm CC_modal_freqs.txt 
	#calculate NJ distances tree
	cat $(ls *.modal_freq | sort -V) $(ls */*.modal_freq | sort -V) > fmdata2.txt
	freqs_2_nj_tree_linux -a -m "fmdata2.txt"  > "freqs.tree"
	cat "freqs.tree" | tail -1 > freqs_only.tree
	mv freqs_only.tree freqs.tree

	cd ..
	Rscript --vanilla tree_heatmap.r "./${D}"
	cd "$D"
	figtree -graphic SVG renamed.tree "${TN}.svg" &> /dev/null
	cd ..
done



