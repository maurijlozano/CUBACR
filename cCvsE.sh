#!/bin/bash
clear

SKIP=0
SKIPC=0
SKIPB=0
SKIPM=0
SKIPCRVR=0
SKIPCOA=0
SKIPT=0
Ite=50
NDS=0
VI=0.2
show_help()
{
echo "
        Usage: cCvsE.sh [-o] [-c] [-h] [-i] [-b] [-m] [-n] [-a] [-t] [-v]

        -o Skips G.Olsen Bootstrapped distances calculation
        -c Skips AA guided codon alingment calculation
        -h Show this Help
        -i n for bootstrap of the optimization algorithm (50 default)
        -b Skips Bootstrapped codon distances calculation
        -m Skips modal sequences calculations
        -n Skips CR and VR modal sequences calculations
        -a Skips CA analysis
        -t Skips Conservation table generation
        -v Defines variability cutoff, must be 0<v<1
"

echo -e '\nGenerates an amino acid guided codon alignment and calcualtes codon usage for conserved protein sequences, both on higly and lowly expressed proteins.\n
Then calculates the distance between modal and bootstraped sequences.\n\n'

exit 1
}

while getopts ":ochbmnatv:i:" option; do
    case "${option}" in
        o) SKIP=1;;
        c) SKIPC=1;;
        h) show_help
            exit 1
            ;;
        i) Ite=$OPTARG;;
        v) VI=$OPTARG;;
        b) SKIPB=1;;
        m) SKIPM=1;;
        n) SKIPCRVR=1;; 
        a) SKIPCOA=1;;
        t) SKIPT=1;;      
        \?) echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done



#if [[ $# -eq 0 ]] ; then
#    echo 'Parameters required, run -h for help'
#    exit 1
#fi


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



#*************************************************************************************************
#************************************    2    ****************************************************
#*************************************************************************************************


if [[ $SKIPC == 0 ]]
then	

	#First, calculate aa guided codon alignments 
	ls -d */*/ | grep 'CC' -v | grep 'PHE' -v | while read D; do
		echo -e "Calculating aa guided codon alignments in ${D} \n"
		ls $D | grep '.fa[s]\?$' | grep '^HEP' -v | grep '^LEP' -v | while read F; do
			FN=$(echo "${F}" | sed 's/.fa.\?$//')
			FP="${D}${F}"
			#remove enters from seqs
			cat "${FP}" | perl -pe 'unless(/^>/){s/\n//g};s/>/\n>/' | sed -e 's/-//g' -e '/^$/d' > "${FP}ed"
			mv "${FP}ed" "${FP}"
			FOP="${D}${FN}/${FN}"
			if [ ! -d "${D}${FN}" ]; then mkdir "${D}${FN}"; fi
			./translatorx_vLocal.pl -i "./${FP}" -o "./${FOP}" &> /dev/null
			echo -e "Calculating ${FN} conserved codon fasta file. \n"
			#removing enters from alignments
			cat "./${FOP}.nt_ali.fasta" | perl -pe 'unless(/^>/){s/\n//g};s/>/\n>/' | sed '/^$/d' > "./${FOP}.nt_ali.fastaed"
			cat "./${FOP}.aa_ali.fasta" | perl -pe 'unless(/^>/){s/\n//g};s/>/\n>/' | sed '/^$/d' > "./${FOP}.aa_ali.fastaed" 
			mv "./${FOP}.nt_ali.fastaed" "./${FOP}.nt_ali.fasta"
			mv "./${FOP}.aa_ali.fastaed" "./${FOP}.aa_ali.fasta"
						
			./codonAlnStats.pl -i "./${FOP}.nt_ali.fasta" -a "./${FOP}.aa_ali.fasta" -o "./${FOP}" > /dev/null
			Rscript --vanilla stats.R "./${FOP}.aa.table" "./${FOP}.nt.table" "./${FOP}" "${VI}" "${FN}" > /dev/null
		done
	done 
	
else
    	echo -e 'Skipping AA guided codon alingment.\n'
fi

#*************************************************************************************************
#************************************    3    ****************************************************
#*************************************************************************************************
if [[ $SKIPT == 0 ]]
then	

	echo -e "Generating conservation table. \n"

	ls -d */ | while read D; do
		cd ${D}
		SN=$(echo "${D}" | sed 's/\/$//')
		echo "Expression/Protein/file AA.ID CODON.ID AA.VAR AA.Length" > "${SN}_CE.table"
		ls */*/*.stats | while read F; do
			echo "${F} $(cat "${F}")" >> "${SN}_CE.table"
		done
		cd ..
		Rscript --vanilla stats.plot.R "./${D}${SN}_CE.table" >/dev/null
	done
else
    	echo -e 'Skipping conservation table generation.\n'
fi


#*************************************************************************************************
#************************************    4    ****************************************************
#*************************************************************************************************


ls -d */ | while read D; do
	#Processing each species!
	SN=$(echo "${D}" | sed 's/\/$//')
	echo -e "Processing ${SN}.\n"
	cd "${D}"
	
	#*********************************************************************************************
	#*********************************************************************************************
	if [[ $SKIPCRVR == 0 ]]
	then	
		echo -e "Generating multifasta and modal sequence of Conserved/variable regions of both HEP and LEP. \n"
		# Calculates modals and concatenated sequences for HEP CR VR and LEP CR VR
		ls -d */ | grep 'CC' -v | grep 'PHE' -v | while read ED; do
			EFN=$(echo "${ED}" | sed 's/\/$//')
			echo -e "\nProcessing ${EFN}.\n"
			#CR
			cat $(ls ${ED}*/* | grep 'aa100C.fas$') > "${SN}_${EFN}_CR.fas"
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
			echo -e "\n"
			perl seq2AvgAA.pl "${D}${SN}_${EFN}_VR.fas"
			echo -e "\n"
			perl modal2seq2.pl "${D}${SN}_${EFN}_CR.modal_freq" "${D}${SN}_${EFN}_CR.aaav"
			perl modal2seq2.pl "${D}${SN}_${EFN}_VR.modal_freq" "${D}${SN}_${EFN}_VR.aaav"
			sed 's/>.*/>'${EFN}'_CR/' "${D}${SN}_${EFN}_CR.modal_freq.fas" > "${D}${SN}_${EFN}_CR.modal_freq.fased"
			sed 's/>.*/>'${EFN}'_VR/' "${D}${SN}_${EFN}_VR.modal_freq.fas" > "${D}${SN}_${EFN}_VR.modal_freq.fased"
			mv "${D}${SN}_${EFN}_CR.modal_freq.fased" "${D}${SN}_${EFN}_CR.modal_freq.fas"
			mv "${D}${SN}_${EFN}_VR.modal_freq.fased" "${D}${SN}_${EFN}_VR.modal_freq.fas"
			cd "${D}"
	 	done
	else
    	echo -e 'Skipping CR and VR modal sequences calculation.\n'
	fi
 	
	#*********************************************************************************************
	#*********************************************************************************************
 	#*************************************************************************************************
 	#*************************************************************************************************
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
    
	#*************************************************************************************************
	#*************************************************************************************************


	#*************************************************************************************************
	#*************************************************************************************************
 	#2 Calculate bootstrapped sequences for the concatenated sequence and then, calcualte distances using Olsen distance.
 	#2.1. calcualte PHE concat.
 	if [[ $SKIPB == 0 ]]
	then	
	 	#tests if PHE file is present
	 	NDS=0
		HEP=0
		LEP=0
			
		if ! ls PHE/*.fa*
		then 
			echo -e "A PHE file located on SN/PHE folder is required for PHE distances calculations.\n"
			PHEEXIST=0
		else	
			#PHE concat
			PHEEXIST=1
			PHE=($(echo $(ls PHE/*.fa* | grep '.modal_freq.fas' -v))) 
			PHEN=$(echo "${PHE}" | sed 's/^PHE\///' | sed 's/.fa.\?$//')
		 	echo ">PHE_concat" > "${PHEN}_concat.fas" 
			grep '>' -v "${PHE}" | sed -e 's/^ATG//' | sed -e 's/[T][AG][AG]$//'| tr -d "\n" > seq.tmp
			echo "ATG"$(cat seq.tmp)"TAA"  >> "${PHEN}_concat.fas"
			rm seq.tmp
		fi &>/dev/null
	 	
	 	declare -a FILESC
	 	FILESC=($(echo $(ls *concat.fas | grep 'HEP')  $(ls *concat.fas | grep 'LEP')  $(ls *concat.fas | grep 'PHE')))
	 	
	 	#counting conditions!!
	 	if ls HEP/*.fa*; then NDS=$(($NDS + 2)); HEP=1 ; echo -e " HEP proteins present.";fi &>/dev/null	
	 	if ls LEP/*.fa*; then NDS=$(($NDS + 2)); LEP=1 ; echo -e " LEP proteins present.";fi &>/dev/null
	 	if ls PHE/*.fa*; then NDS=$(($NDS + 1)); echo -e " PHE proteins present.";fi &>/dev/null	
	 	echo -e "\n"
	 	
	 	#2.2. PHE exists
		 Rscript --vanilla ../codonboots.R "${Ite}" "${FILESC[0]}" "${FILESC[1]}" "${FILESC[2]}" "${FILESC[3]}" "${FILESC[4]}" > /dev/null
		 	
		ls *.boot.fas | while read BOOT; do
			BN=$(echo "${BOOT}" | sed 's/.fas$//')
			compile_codon_freqs  < "${BOOT}" > "${BN}.freqs"
		done

		#NDS depends on weather HEP CR/VR, LEP CR/V or PHE exist. 
		freqs_2_dists $(echo $(ls *.boot.freqs | grep -v '^PHE') $(ls *.boot.freqs | grep '^PHE'))  > distances  #esto calcula la distancia tipo2 para todos contra todos-> como está puesto nos interesan los elementos i de 1:n y j de i+n a 2n. Esas serian todas las combinaciones de distancias entre 1i y 2i (11-21; 11-22; ....;12-21;...))
		Rscript --vanilla ../distproc.R "${Ite}" "${NDS}" "${HEP}" "${LEP}" "${PHEEXIST}" &>/dev/null
	else
    	echo -e 'Skipping Bootstrapped sequences distance calculation.\n'
	fi 
	#*************************************************************************************************
#*************************************************************************************************
	
	
 	cd ..
done



#*************************************************************************************************
#************************************    5    ****************************************************
#*************************************************************************************************


#3. CA mapped to all proteins for each species - add PHE and CR-PHE! Here, look for GNM, PHE data and use it on CA 
#3.1 calcualte modal CC,single y PHE sequences!
	#calculate modals and modal seqs.


	
#*************************************************************************************************
#*************************************************************************************************
if [[ $SKIPM == 0 ]]
	then	
		
		
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
else
   	echo -e 'Skipping CC, PHE, HEP and LEP modal sequences calculation.\n'
fi

#*************************************************************************************************
#*************************************************************************************************


if [[ $SKIPCOA == 0 ]]
then	

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
		Rscript --vanilla graph_coa.R "./${D}" genes.coa &>/dev/null
		Rscript --vanilla graph_codons.R "./${D}" codon.coa &>/dev/null
		echo -e "\n"
	done
else
    	echo -e 'Skipping CA analysis.\n'
fi 	

# to do Comparison beetween species??

