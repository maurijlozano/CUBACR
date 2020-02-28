#!/bin/bash
clear
VI=0.5
show_help()
{
echo "
        Usage: SimulCR_VR_IQTREE.sh [-h] [-v]

        -h Show this Help
        -v Defines variability cutoff, must be 0<v<1
"

echo -e '\nGenerates sequences under the same evolutionary model of HEP and LEP proteins both on highly and lowly expressed proteins, using IQTREE to infer phylogeny under the GY M0 model.\n\n'

exit 1
}

while getopts ":hv:" option; do
    case "${option}" in
        h) show_help
            exit 1
            ;;
        v) VI=$OPTARG;;
        \?) echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
        :)  echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

WD=$(pwd)
#First, converts fasta alignments to sequential phylip format
ls -d */*/*/ | grep 'CC' -v | grep 'PHE' -v | while read -r D; do
	cd "$WD"
	echo -e "Analizing ${D}."
	protID=$(echo $D | sed -E 's/^[^/]*\/[^/]*\/([^/]*)\/$/\1/' | sed -E 's/([^_]*_[^_]*).*/\1/')
	cd "${D}"

	#Modeltest
	if [ ! -f "nt.phy.treefile" ]; then
		#generating nt.ali
		input=$(ls *.nt_ali.fasta | grep 'sim' -v)
		ALNL=$(sed -n '2p' ${input} | wc -m)
		ALNL=$(echo $((${ALNL} - 4)))
		ALNS=$(grep -v '>' ${input} | wc -l)
		touch nt.phy
		echo ' '${ALNS}' '${ALNL} > nt.phy
		while read -r line || [ -n "$line" ]; do
			if [[ $line == '>'* ]]
			then
				echo -n $line | sed -E 's/>([^_]*_[^_]*).*/\1  /' >> nt.phy
			else
				echo $line | sed -E 's/(.{3})/\1 /g' | sed -e 's/TAG//g' -e 's/TAA//g' -e 's/TGA//g' >> nt.phy
			fi
		done < "$input"
		{
		iqtree -s nt.phy -st CODON -m GY+F3x4 -T 2 #AUTO
		} || {
			echo "An error occurred..."
			exit
		}
	else
		echo "Phylogenetic tree and parameters already estimated."
	fi

	if [ ! -f "nt.phy.treefile" ]; then
		echo -e "Something failed..."
		exit
	fi

	#Evolver
	#Generate MCcodon.dat for palm-evolver
	touch MCcodon.dat
	echo -e '0\n13147' > MCcodon.dat
	nnucl=$(head -n 1 nt.phy | sed -E 's/^ *([^ ]*) ([^ ]*).*/\2/') 
	ncodons=$(echo $(( $nnucl / 3 )))
	line=$(echo "$(head -n 1 nt.phy | sed -E 's/^ *([^ ]*) ([^ ]*).*/\1/') ${ncodons} 1")
	echo $line >> MCcodon.dat
	line=$(grep 'Total tree length' nt.phy.log)
	echo $line | sed -E 's/^[^:]*: *([0-9\.])[ \t]*/\1/' >> MCcodon.dat
	line=$(head -n 1 nt.phy.treefile)
	echo $line >> MCcodon.dat
	grep -A 11 'FINALIZING TREE SEARCH' nt.phy.log > nt.phy.log2
	line=$(grep "Nonsynonymous/synonymous ratio (omega):" nt.phy.log2)
	echo $line | sed -E 's/^Nonsynonymous\/synonymous ratio \(omega\):[ \t]*//' >> MCcodon.dat
	line=$(grep "Transition/transversion ratio (kappa): " nt.phy.log2)
	echo $line | sed -E 's/^Transition\/transversion ratio \(kappa\): [ \t]*//' >> MCcodon.dat
	
	grep 'State frequencies:' nt.phy.iqtree -A 16 | tail -n +2 > modelFreqs.txt
	sed -E 's/^ *pi\(([AGTC]{3})\) = ([0-9\.]*) *pi\(([AGTC]{3})\) = ([0-9\.]*) *pi\(([AGTC]{3})\) = ([0-9\.]*) *pi\(([AGTC]{3})\) = ([0-9\.]*).*$/\1 \2\n\3 \4\n\5 \6\n\7 \8/' modelFreqs.txt | sed -E 's/^ *pi\(([AGTC]{3})\) = ([0-9\.]*)/\1 \2/' > modelFreqs2.txt
	mv modelFreqs2.txt modelFreqs.txt
	cd "$WD"
	Rscript --vanilla formatCodonfreqs.r "./${D}"
	cd "$D"
	cat MCcodon.dat modelFreqs.txt > MCcodon2.dat
	mv MCcodon2.dat MCcodon.dat
	echo -e "0\n// end of file." >> MCcodon.dat

	if [ ! -f "mc.paml" ]; then
		if [[ $(command -v evolver) ]]
			then 
			echo -e "Generating simulated dataset."
			{
				evolver 6 MCcodon.dat &> /dev/null
			} || {
				echo "An error with Paml evolver occurred!"
				rm mc.paml &> /dev/null
				exit
			}
		elif [[ $(command -v paml-evolver) ]]
			then
			echo -e "Generating simulated dataset."
			{
				paml-evolver 6 MCcodon.dat &> /dev/null
			} || {
				echo "An error with Paml evolver occurred!"
				rm mc.paml &> /dev/null
				exit
			}
		else
			echo "Paml evolver not found!"
			exit
		fi
	else
		echo "simulation already run."
	fi


	#Extract C and VR
	if [[ ! -f "${protID}.sim.aa100C.fas" && -f "mc.paml" ]]; then
		echo "Extracting CR and VR for the simulated data."
		#phy to fasta
		tail -n +5 mc.paml | sed -E 's/^([^ ]*)  />\1\n/' | sed -e 's/ //g' > ${protID}.sim #nombre del archivo
		AD=$(pwd)
		cd "$WD"
		./translatorx_vLocal.pl -i "${AD}/${protID}.sim" -o "${AD}/${protID}.sim" &> /dev/null
		#Generates all the aa_ali.fasta and nt_ali.fasta
		cd "$AD"
		
		#removing enters from alignments
		cat "./${protID}.sim.nt_ali.fasta" | perl -pe 'unless(/^>/){s/\n//g};s/>/\n>/' | sed '/^$/d' > "./${protID}.sim.nt_ali.fastaed"
		cat "./${protID}.sim.aa_ali.fasta" | perl -pe 'unless(/^>/){s/\n//g};s/>/\n>/' | sed '/^$/d' > "./${protID}.sim.aa_ali.fastaed" 
		mv "./${protID}.sim.nt_ali.fastaed" "./${protID}.sim.nt_ali.fasta"
		mv "./${protID}.sim.aa_ali.fastaed" "./${protID}.sim.aa_ali.fasta"
		
		#extracting cr and vr
		cd "$WD"
		./codonAlnStats.pl -i "${AD}/${protID}.sim.nt_ali.fasta" -a "${AD}/${protID}.sim.aa_ali.fasta" -o "${AD}/${protID}.sim" &> /dev/null
		Rscript --vanilla stats.R "${AD}/${protID}.sim.aa.table" "${AD}/${protID}.sim.nt.table" "${AD}/${protID}.sim" "${VI}" "${protID}" &> /dev/null
		echo "CR and VR correctly extracted"
	elif [ ! -f "mc.paml" ]; then
		echo "mc.paml file missing!"
		exit 1
	else
		echo "CR and VR already extracted."
		cd "$WD"
	fi


	echo -e "\n"
done



