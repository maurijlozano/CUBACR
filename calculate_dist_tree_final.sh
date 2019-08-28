#!/bin/bash
clear

show_help()
{
echo "
        Usage: calculate_dist_tree_final.sh [-h]]

        -h Show this Help

"

echo -e '\nGenerates a NJ distances tree of CUF for Ci, Singletons, and expression sets.\n\n'

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
	#calculate NJ distances tree
	ls *.modal_freq | sort -V | sed -e 's/.modal_freq$//' -e 's/[^_]*_//' > CCindex2.txt
	ls */*.modal_freq | sort -V | grep '[Pp][Hh[Ee]' -v  | grep '[Gg][Nn[Mm]' -v | sed -e 's/.modal_freq$//' -e 's/\([^\/]*\)\///' -e 's/\([^_]*\)_.*/\1/' >> CCindex2.txt
	cat $(ls *.modal_freq | sort -V) $(ls */*.modal_freq | sort -V | grep 'PHE' -v | grep '[Gg][Nn[Mm]' -v ) > fmdata22.txt
	freqs_2_nj_tree_linux -a -m "fmdata22.txt"  > "freqs2.tree"
	cat "freqs2.tree" | tail -1 > freqs_only.tree
	mv freqs_only.tree freqs2.tree
	freqs_2_dists "fmdata22.txt"  > "freqs2.dist"
	cd ..
	Rscript --vanilla tree_final.r "./${D}"
	cd "$D"
	rm fmdata22.txt CCindex2.txt
	cd ..	
done



