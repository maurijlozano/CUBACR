#!/bin/bash
clear

SKIP=0
SKIPC=0
Ite=50
show_help()
{
echo "
        Usage: cCvsE.sh [-o] [-c] [-h] [-b]

        -o Skips G.Olsen Bootstrapped distances calculation
        -c Skips AA guided codon alingment calculation
        -h Show this Help
        -b n for bootstrap of the optimization algorithm (50 default)
"

echo -e '\nGenerates an amino acid guided codon alignment and calcualtes codon usage for conserved protein sequences, both on higly and lowly expressed proteins.\n
Then calculates the distance between modal and bootstraped sequences.\n\n'

exit 1
}

while getopts ":ochb:" option; do
    case "${option}" in
        o) SKIP=1;;
        c) SKIPC=1;;
        h) show_help
            exit 1
            ;;
        b) Ite=$OPTARG;;
        \?) echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done



if [[ $# -eq 0 ]] ; then
    echo 'Parameters required, run -h for help'
    exit 1
fi


WD=$(pwd)

#*************************************************************************************************
#************************************    1    ****************************************************
#*************************************************************************************************


#eliminate spaces from folder names!
echo -e "Eliminating spaces and undesscores from folder names!\n"
ls -d */ | while read D; do
	if [[ $D =~ ' ' || $D =~ '_' ]]; then DN=$(echo "${D}" | sed 's/ //g' | sed 's/_/-/g');	mv "${D}" "${DN}";D=$DN;fi
	cd "${D}"
	ls -d */ | while read E; do
		if [[ $E =~ ' ' || $E =~ '_' ]]; then EN=$(echo "${E}" | sed 's/ //g' | sed 's/_/-/g');	mv "${E}" "${EN}";E=$EN;fi
		cd "${E}"
		ls | while read F; do
				if [[ $F =~ ' ' ]]; then FN=$(echo "${F}" | sed 's/ //g');	mv "${F}" "${FN}";fi
		done
		cd ..
	done
	cd ..
done 2>/dev/null




if [[ $SKIPC == 0 ]]
then	
#*************************************************************************************************
#************************************    2    ****************************************************
#*************************************************************************************************



	#First, calculate aa guided codon alignments 
	ls -d */*/ | grep 'CC' -v | grep 'PHE' -v | while read D; do
		echo -e "Calculating aa guided codon alignments in ${D} \n"
		ls $D | grep '.fa[s]\?$' | grep '^HEP' -v | grep '^LEP' -v | while read F; do
			FN=$(echo "${F}" | sed 's/.fa.\?$//')
			FP="${D}${F}"
			#remove enters from seqs
			cat "${FP}" | perl -pe 'unless(/^>/){s/\n//g};s/>/\n>/' | sed '/^$/d' > "${FP}ed"
			mv "${FP}ed" "${FP}"
			FOP="${D}${FN}/${FN}"
			if [ ! -d "${D}${FN}" ]; then mkdir "${D}${FN}"; fi
			./translatorx_vLocal.pl -i "./${FP}" -o "./${FOP}" &> /dev/null
			echo -e "Calculating ${FN} conserved codon fasta file. \n"
			./codonAlnStats.pl -i "./${FOP}.nt_ali.fasta" -a "./${FOP}.aa_ali.fasta" -o "./${FOP}" > /dev/null
			Rscript --vanilla stats.R "./${FOP}.aa.table" "./${FOP}.nt.table" "./${FOP}" > /dev/null
		done
	done 

	echo -e "Generating conservation table. \n"

	#*************************************************************************************************
#************************************    3    ****************************************************
#*************************************************************************************************


	ls -d */ | while read D; do
		cd ${D}
		SN=$(echo "${D}" | sed 's/\/$//')
		echo "Expression/Protein/file AA.ID CODON.ID At.leats.5.AA" > "${SN}_CE.table"
		ls */*/*.stats | while read F; do
			echo "${F} $(cat "${F}")" >> "${SN}_CE.table"
		done
		cd ..
	done
else
    	echo -e 'Skipping AA guided codon alingment.\n'
fi




#*************************************************************************************************
#************************************    4    ****************************************************
#*************************************************************************************************

echo -e "Generating multifasta and modal sequence of Conserved/variable regions of both HEP and LEP. \n"

ls -d */ | while read D; do
	#Processing each species!
	SN=$(echo "${D}" | sed 's/\/$//')
	echo -e "Processing ${SN}.\n"
	cd "${D}"
	ls -d */ | grep 'CC' -v | grep 'PHE' -v | while read ED; do
		EFN=$(echo "${ED}" | sed 's/\/$//')
		echo -e "\nProcessing ${EFN}.\n"
		#CR
		cat $(ls ${ED}*/* | grep 'nt100C.fas$') > "${SN}_${EFN}_CR.fas"
		#make concatenated sequence..
		echo ">${EFN}_CR_concat" > "${SN}_${EFN}_CR_concat.fas" 
		grep '>' -v "${SN}_${EFN}_CR.fas" | sed -e 's/^ATG//' | sed -e 's/[T][AG][A]$//'| tr -d "\n" > seq.tmp
		echo "ATG"$(cat seq.tmp)"TAA"  >> "${SN}_${EFN}_CR_concat.fas"
		#VR
		cat $(ls ${ED}*/* | grep 'aaV.fas$') > "${SN}_${EFN}_VR.fas"
		#make concatenated sequence..
		echo ">${EFN}_VR_concat" > "${SN}_${EFN}_VR_concat.fas" 
		grep '>' -v "${SN}_${EFN}_VR.fas" | sed -e 's/^ATG//' | sed -e 's/[T][AG][A]$//'| tr -d "\n" > seq.tmp
		echo "ATG"$(cat seq.tmp)"TAA"  >> "${SN}_${EFN}_VR_concat.fas"
		rm seq.tmp
		
		#calculate modals and modal seqs. for CR and VR
		compile_codon_counts  < "${SN}_${EFN}_CR.fas" | modal_usage_from_counts -a $(nproc) > "${SN}_${EFN}_CR.modal_freq"
		compile_codon_counts  < "${SN}_${EFN}_VR.fas" | modal_usage_from_counts -a $(nproc) > "${SN}_${EFN}_VR.modal_freq"
		cd ..
		perl seq2AvgAA.pl "${D}${SN}_${EFN}_CR.fas"
		perl seq2AvgAA.pl "${D}${SN}_${EFN}_VR.fas"
		perl modal2seq2.pl "${D}${SN}_${EFN}_CR.modal_freq" "${D}${SN}_${EFN}_CR.aaav"
		perl modal2seq2.pl "${D}${SN}_${EFN}_VR.modal_freq" "${D}${SN}_${EFN}_VR.aaav"
		sed 's/>.*/>'${EFN}'_CR/' "${D}${SN}_${EFN}_CR.modal_freq.fas" > "${D}${SN}_${EFN}_CR.modal_freq.fased"
		sed 's/>.*/>'${EFN}'_VR/' "${D}${SN}_${EFN}_VR.modal_freq.fas" > "${D}${SN}_${EFN}_VR.modal_freq.fased"
		mv "${D}${SN}_${EFN}_CR.modal_freq.fased" "${D}${SN}_${EFN}_CR.modal_freq.fas"
		mv "${D}${SN}_${EFN}_VR.modal_freq.fased" "${D}${SN}_${EFN}_VR.modal_freq.fas"
		cd "${D}"
 	done
 		
 	#1. calculate olsen bootstrapped distances between multifasta sequences CR_HEP, CR_LEP and PHE.
 	if [[ ${SKIP} == 0 ]]
	then	 
		declare -a FILES #array!!
	 	FILES=($(echo $(ls *_CR.fas  *_VR.fas PHE/*.fa* | grep '.modal_freq.fas$' -v)))   #PHE
	 	len=${#FILES[@]}
	 	echo -e "\nCalculating distances.\n"
	 	for i in "${FILES[@]}";do compile_codon_counts  < "${i}" > "${i}.counts"; done
		 		 	
	 	evaluate_replicon_distance $(echo $(ls *.counts PHE/*.counts | grep 'eval_replicon_dist_temp_' -v)) > modal_dist1.dist
	 	evaluate_replicon_distance_2m $(echo $(ls *.counts PHE/*.counts))  > modal_dist2.dist
		rm *.counts PHE/*.counts
	else
    	echo -e '\nSkipping G. Olsen bootstrapped distance calculations.\n'
    fi
	

 	#2 Calculate bootstrapped sequences for the concatenated sequence and then, calcualte distances using Olsen distance.
 	#2.1. calcualte PHE concat.
 	
 	#tests if PHE file is present
 	NDS=2
	if ! ls PHE/*.fa*
	then 
		echo -e "A PHE file located on SN/PHE folder is required for PHE distances calculations.\n"
		PHEEXIST=0
	else	
		#PHE concat
		PHEEXIST=1
		NDS=3
		PHE=($(echo $(ls PHE/*.fa* | grep '.modal_freq.fas' -v))) 
		PHEN=$(echo "${PHE}" | sed 's/^PHE\///' | sed 's/.fa.\?$//')
	 	echo ">PHE_concat" > "${PHEN}_concat.fas" 
		grep '>' -v "${PHE}" | sed -e 's/^ATG//' | sed -e 's/[T][AG][AG]$//'| tr -d "\n" > seq.tmp
		echo "ATG"$(cat seq.tmp)"TAA"  >> "${PHEN}_concat.fas"
		rm seq.tmp
	fi 2>/dev/null
 	
 	declare -a FILESC
 	FILESC=($(echo $(ls *concat.fas | grep 'HEP')  $(ls *concat.fas | grep 'LEP')  $(ls *concat.fas | grep 'PHE')))
 	
 	if [ ${PHEEXIST} == 1 ]
 	then 
	 	#2.2. PHE exists
	 	Rscript --vanilla ../codonboots.R "${FILESC[0]}" "${FILESC[1]}" "${FILESC[2]}" "${FILESC[3]}" "${Ite}" "${FILESC[4]}"
	else
		#2.3. PHE does not exists
	 	Rscript --vanilla ../codonboots.R "${FILESC[0]}" "${FILESC[1]}" "${FILESC[2]}" "${FILESC[3]}" "${Ite}"
	fi
	
	ls *.boot.fas | while read BOOT; do
		BN=$(echo "${BOOT}" | sed 's/.fas$//')
		compile_codon_freqs  < "${BOOT}" > "${BN}.freqs"
	done

	#aca la idea es mandar al distproc.r un argumento de si van distancias entre 2 o 4 datasets
	freqs_2_dists $(echo $(ls *.boot.freqs | grep -v '^PHE') $(ls *.boot.freqs | grep '^PHE'))  > distances  #esto calcula la distancia tipo2 para todos contra todos-> como está puesto nos interesan los elementos i de 1:n y j de i+n a 2n. Esas serian todas las combinaciones de distancias entre 1i y 2i (11-21; 11-22; ....;12-21;...))
	Rscript --vanilla ../distproc.R "${Ite}" "${NDS}"


 	cd ..
done



#*************************************************************************************************
#************************************    5    ****************************************************
#*************************************************************************************************


#3. CA mapped to all proteins for each species - add PHE and CR-PHE! Here, look for GNM, PHE data and use it on CA 
#3.1 calcualte modal CC,single y PHE sequences!
	#calculate modals and modal seqs.
echo -e "Calculating CC and singleton modals\n"
ls -d */*/* | grep '[Ff][Aa][Ss]\?$' | grep 'HEP' -v | grep 'LEP' -v | grep '.modal_freq.fas' -v | while read D; do 		
	BN=$(echo "${D}" | sed 's/.fas\?$//')
	compile_codon_counts  < "${D}" | modal_usage_from_counts -a $(nproc) > "${BN}.modal_freq"

	perl seq2AvgAA.pl "${D}"
	perl modal2seq2.pl "${BN}.modal_freq" "${BN}.aaav"
	echo -e "\n"
done

#3.2 HEP/LEP multifasta and modalseq
echo -e "Calculating HEP and LEP modals\n"
ls -d */*/ | grep -E 'HEP|LEP' | while read D; do 		
	cd $D
	NAME=$(echo "${D}" | sed 's/^[^\/]*\///' | sed 's/\///')
	head -n2 $(ls *.fas | grep '.modal_freq.fas' -v | grep '^HEP' -v | grep '^LEP' -v) | grep -E '==>|^$' -v | sed 's/-//g' > "${NAME}.fas"
	compile_codon_counts < "${NAME}.fas" | modal_usage_from_counts -a $(nproc) > "${NAME}.modal_freq"
	cd $WD
	perl seq2AvgAA.pl "${D}${NAME}.fas"
	perl modal2seq2.pl "${D}${NAME}.modal_freq" "${D}${NAME}.aaav"
	echo -e "\n"
done


#3.3 Generate multifasta for CA analysis. Genomic_sequences, modals, expression, CR
#Name of the files must be 'set-id'_'any-text' with an under score...
echo -e "Runing CA"
ls -d */ | while read D; do
	cd ${D}
	cat $(ls -d */* | grep '[Ff][Aa][Ss]\?$' | grep -E 'CC/genoma|CC/GNM|CC/genome' | grep 'modal_freq.fas$' -v) | sed -e 's/[ \t_]/-/g' -e 's/\(>.*\)/\1_Gene/' > genes.fas
	cat genes.fas $(ls -d */*modal_freq.fas | grep -E 'CC/genoma|CC/GNM|CC/genome') $(ls -d */* -v | grep 'modal.freq.fas$' | grep '^CC/genoma' -v | grep -E '^CC/GNM' -v) $(ls | grep 'concat.fas$') $(ls | grep 'CR.modal_freq.fas$') $(ls | grep 'VR.modal_freq.fas$') > coa.fas
	sed -e 's/[^_]*_\(Gene\)$/>\1/' -e 's/^>[^\/]*\/[^\/]*\/\([^_]*\).*/>\1/' coa.fas > coa.fased
	mv coa.fased coa.fas
	codonw coa.fas -nomenu -nowarn -silent -coa_rscu -noblk
    rm hilo.coa fop.coa cusort.coa coa.out coa.blk cbi.coa cai.coa eigen.coa coa_raw coa.fas genes.fas
    sed -e 's/^[0-9]*_//' -e 's/ $//' genes.coa > genes.coaed
    mv genes.coaed genes.coa
    touch coa_var_axis1_2.txt
	echo 'varC1 varC2' > coa_var_axis1_2.txt
	sed '17q;d' summary.coa | sed -e 's/   / /g' | cut -f3,7 -d' ' >> coa_var_axis1_2.txt
	sed -e 's/\+//g' coa_var_axis1_2.txt > coa_var_axis1_2.txted
	mv coa_var_axis1_2.txted coa_var_axis1_2.txt
	cd ..
    Rscript --vanilla graph_coa.R "./${D}" genes.coa
    Rscript --vanilla graph_codons.R "./${D}" codon.coa
    echo -e "\n"
done
 	

# to do Comparison beetween species??

