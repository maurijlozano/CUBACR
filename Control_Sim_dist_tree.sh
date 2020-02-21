#!/bin/bash
clear

show_help()
{
echo "
        Usage: Control_Sim_dist_tree [-h]]

        -h Show this Help

"

echo -e '\nCalculate modal sequences for simulated data and generates a NJ distances tree of CUF for Ci, Singletons, expression sets and simulated data.\n\n'

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

#generate modals

# Calculating CR and VR modals...

ls -d */ | while read D; do
	#Processing each species!
	SN=$(echo "${D}" | sed 's/\/$//')
	echo -e "Processing ${SN}.\n"
	cd "${D}"
	
	echo -e "Generating multifasta and modal sequence of simulated Conserved/variable regions of both HEP and LEP. \n"
	# Calculates modals and concatenated sequences for CR VR and CR VR
	ls -d */ | grep 'CC' -v | grep 'PHE' -v | while read ED; do
		EFN=$(echo "${ED}" | sed 's/\/$//')
		echo -e "\nProcessing ${EFN}.\n"
		#CR
		cat $(ls ${ED}*/* | grep 'sim.aa100C.fas$') > "${SN}_${EFN}_sim_CR.fas"
		#VR
		cat $(ls ${ED}*/* | grep 'aaV.fas$') > "${SN}_${EFN}_sim_VR.fas"
		#calculate modals and modal seqs. for CR and VR
		compile_codon_counts  < "${SN}_${EFN}_sim_CR.fas" | modal_usage_from_counts -a $(nproc) > "${SN}_${EFN}_sim_CR.modal_freq"
		compile_codon_counts  < "${SN}_${EFN}_sim_VR.fas" | modal_usage_from_counts -a $(nproc) > "${SN}_${EFN}_sim_VR.modal_freq"
		cd ..
		perl seq2AvgAA.pl "${D}${SN}_${EFN}_sim_CR.fas"
		echo -e "\n"
		perl seq2AvgAA.pl "${D}${SN}_${EFN}_sim_VR.fas"
		echo -e "\n"
		perl modal2seq2.pl "${D}${SN}_${EFN}_sim_CR.modal_freq" "${D}${SN}_${EFN}_sim_CR.aaav"
		perl modal2seq2.pl "${D}${SN}_${EFN}_sim_VR.modal_freq" "${D}${SN}_${EFN}_sim_VR.aaav"
		sed 's/>.*/>'${EFN}'_sim_CR/' "${D}${SN}_${EFN}_sim_CR.modal_freq.fas" > "${D}${SN}_${EFN}_sim_CR.modal_freq.fased"
		sed 's/>.*/>'${EFN}'_sim_VR/' "${D}${SN}_${EFN}_sim_VR.modal_freq.fas" > "${D}${SN}_${EFN}_sim_VR.modal_freq.fased"
		mv "${D}${SN}_${EFN}_sim_CR.modal_freq.fased" "${D}${SN}_${EFN}_sim_CR.modal_freq.fas"
		mv "${D}${SN}_${EFN}_sim_VR.modal_freq.fased" "${D}${SN}_${EFN}_sim_VR.modal_freq.fas"
		cd "${D}"
	done
	cd ..
done
echo -e "\n"


#Generate simulated HEP and LEP modals

ls -d */ | while read D; do
	cd "${D}"
	ls -d */ | grep 'CC' -v | grep 'PHE' -v | while read E; do
		cd "${E}"
		file="$(echo ${E} | sed -E 's/\/$//')_sim"
		if [ ! -f "$file.fas" ]; then 
			touch "$file.fas"
		else
			rm "$file.fas"
			touch "$file.fas"
		fi
		ls -d */ | while read P; do
			protID=$(echo "$P" | sed -E 's/^([^/]*).*\/$/\1/')
			#protID=$(echo "$P" | sed -E 's/^([^/]{10}).*\/$/\1/')
			grep -A 1 "${protID}" ${P}${protID}.sim >> $file.fas
		done
		cd ../../
		compile_codon_counts  < "${D}${E}${file}.fas" | modal_usage_from_counts -a $(nproc) > "${D}${E}${file}.modal_freq"
		perl seq2AvgAA.pl "${D}${E}${file}.fas"
		perl modal2seq2.pl "${D}${E}${file}.modal_freq" "${D}${E}${file}.aaav"
		sed 's/>.*/>'${file}'/' "${D}${E}${file}.modal_freq.fas" > "${D}${E}${file}.modal_freq.fased"
		mv "${D}${E}${file}.modal_freq.fased" "${D}${E}${file}.modal_freq.fas" 
		echo -e "\n"
		cd "${D}"
	done
	cd ..
done

echo -e "\n"
#generate tree control
ls -d */ | while read D; do
	cd "$D"
	#calculate NJ distances tree
	ls *.modal_freq | sort -V | sed -e 's/.modal_freq$//' -e 's/[^_]*_//' > CCindex2.txt
	ls */*.modal_freq | sort -V | grep '[Pp][Hh[Ee]' -v  | grep '[Gg][Nn[Mm]' -v | sed -e 's/.modal_freq$//' -e 's/\([^\/]*\)\///' -e 's/\([^_]*\)_[^s^y^m].*/\1/' >> CCindex2.txt
	cat $(ls *.modal_freq | sort -V) $(ls */*.modal_freq | sort -V | grep 'PHE' -v | grep '[Gg][Nn[Mm]' -v ) > fmdata22.txt
	freqs_2_nj_tree_linux -a -m "fmdata22.txt"  > "freqs3.tree"
	cat "freqs3.tree" | tail -1 > freqs_only.tree
	mv freqs_only.tree freqs3.tree
	freqs_2_dists "fmdata22.txt"  > "freqs3.dist"
	cd ..
	Rscript --vanilla tree_control.r "./${D}"
	cd "$D"
	rm fmdata22.txt CCindex2.txt
	cd ..	
done

#control tree 2
ls -d */ | while read D; do
	cd "$D"
	#calculate NJ distances tree
	ls *.modal_freq | grep 'sim' | sort -V | sed -e 's/.modal_freq$//' -e 's/[^_]*_//' > CCindex2.txt
	ls */*.modal_freq | sort -V | grep '[Pp][Hh[Ee]' -v  | grep '[Gg][Nn[Mm]' -v | grep 'HEP.m' -v | grep 'LEP.m' -v | sed -e 's/.modal_freq$//' -e 's/\([^\/]*\)\///' -e 's/\([^_]*\)_[^s^y^m].*/\1/' >> CCindex2.txt
	cat $(ls *.modal_freq | grep 'sim' | sort -V) $(ls */*.modal_freq | sort -V | grep '[Pp][Hh[Ee]' -v  | grep '[Gg][Nn[Mm]' -v | grep 'HEP.m' -v | grep 'LEP.m' -v ) > fmdata22.txt
	freqs_2_nj_tree_linux -a -m "fmdata22.txt"  > "freqs4.tree"
	cat "freqs4.tree" | tail -1 > freqs_only.tree
	mv freqs_only.tree freqs4.tree
	freqs_2_dists "fmdata22.txt"  > "freqs4.dist"
	cd ..
	Rscript --vanilla tree_control2.r "./${D}"
	cd "$D"
	rm fmdata22.txt CCindex2.txt
	cd ..	
done

