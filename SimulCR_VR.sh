#!/bin/bash
clear
VI=0.5
show_help()
{
echo "
        Usage: SimulCR_VR.sh [-h] [-v]

        -h Show this Help
        -v Defines variability cutoff, must be 0<v<1
"

echo -e '\nGenerates sequences under the same evolutionary model of HEP and LEP proteins both on highly and lowly expressed proteins.\n\n'

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
	echo -e "Testing evolutionary models for ${D}."
	protID=$(echo $D | sed -E 's/^[^/]*\/[^/]*\/([^/]*)\/$/\1/' | sed -E 's/([^_]*_[^_]*).*/\1/')
	cd "${D}"
	seqname=$(grep '>' $(ls *.aa_ali.fasta | grep 'sim' -v ) | head -n 1 | sed -E 's/>([^_]*_[^_]*).*/\1/')

	if [[ ! -f aa.phy ]]
	then
		#generating aa.ali
		input=$(ls *.aa_ali.fasta | grep 'sim' -v )
		ALNL=$(sed -n '2p' ${input} | wc -m)
		ALNL=$(echo $((${ALNL} - 1)))
		ALNS=$(grep -v '>' ${input} | wc -l)
		echo ' '${ALNS}' '${ALNL} > aa.phy
		while read -r line || [ -n "$line" ]; do
			if [[ $line == '>'* ]]
			then
				seqlabel="$(echo -n $line | sed -E 's/>([^_]*_[^_]*).*/\1/')"
				echo -n "$seqlabel  " >> aa.phy
			else
				echo $line | sed -E 's/(.{10})/\1 /g' | sed -E 's/X/-/g'  >> aa.phy
			fi
		done < "$input"

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
	else
		echo "Alignments files already generated..."
	fi

	#Modeltest
	if [ ! -f "modeltest.out" ]; then
		#Now, we have to run modeltest-ng and codonPhyml to make treefile
		{
			modeltest-ng -i aa.phy -d aa -o modeltest -m "DAYHOFF,LG,DCMUT,JTT,MTREV,WAG,RTREV,CPREV,VT,BLOSUM62,MTMAM,MTART,HIVB,HIVW" --force &> /dev/null 
		} || {
			echo "  > phyml  -i aa.phy -d aa -m LG -f m -v 0 -a e -c 4 -o tlr" > modeltest.out
			echo "Modeltest produced an error, using LG+G4 as model!."
		}
	else
		echo "Modeltest already run for this protein."
	fi


	#Generate tree with codonphyml
	if [ ! -f "aa.phy_codonphyml_tree.txt" ]; then
		line=$(grep '> phyml' modeltest.out  | sed -n 1p)
		line=$(echo $line | sed -E 's/^> phyml//' | sed -e 's/-a 0/-a e/')
		echo -e "Making phylogenetic tree."
		{
			codonphyml $line &> /dev/null
		} || {
			echo "Analysis not completed! Deleting aa.phy_codonphyml_tree.txt..."
			rm aa.phy_codonphyml_tree.txt &> /dev/null
			exit
		}
	else
		echo "Phylogenetic tree already constructed for this protein."
	fi

	#Codeml
	#The next step is to create codeml config filename
	touch codeml.ctl
	echo "      seqfile = ./nt.phy" > codeml.ctl
	echo "      treefile = ./aa.phy_codonphyml_tree.txt" >> codeml.ctl
	echo "      outfile = ./results.mlc" >> codeml.ctl
	echo "      noisy = 0" >> codeml.ctl
	echo "      verbose = 1" >> codeml.ctl
	echo "      runmode = 0" >> codeml.ctl
	echo "      seqtype = 1" >> codeml.ctl
	echo "      CodonFreq = 2" >> codeml.ctl
	echo "      ndata = 1" >> codeml.ctl
	echo "      clock = 0" >> codeml.ctl
	echo "      aaDist = 0" >> codeml.ctl
	echo "      model = 0" >> codeml.ctl
	echo "      NSsites = 0" >> codeml.ctl
	echo "      icode = 0" >> codeml.ctl
	echo "      fix_kappa = 0" >> codeml.ctl
	echo "      kappa = 2.05154" >> codeml.ctl
	echo "      fix_omega = 0" >> codeml.ctl
	echo "      omega = 1" >> codeml.ctl
	echo "      getSE = 0" >> codeml.ctl
	echo "      RateAncestor = 0" >> codeml.ctl
	echo "      Small_Diff = .5e-6" >> codeml.ctl
	echo "      cleandata = 0" >> codeml.ctl
	echo "      fix_blength = 1" >> codeml.ctl
	echo "      method = 0" >> codeml.ctl

	# Optimizing codon model M0.

	if [ ! -f "results.mlc" ]; then
		echo -e "Optimizing codon model M0."
		{ 
			codeml codeml.ctl &> /dev/null

		} || { 
			echo "Analysis not completed! Deleting results.mlc..."
			rm results.mlc
			exit
		}
	else			
		echo "Codon model already optimized."
	fi



	#Evolver
	#Generate MCcodon.dat for palm-evolver
	touch MCcodon.dat
	echo -e '0\n13147' > MCcodon.dat
	line=$(grep 'ns =' results.mlc)
	echo $line | sed -E 's/^ns = ([0-9]{1,3}) ls = ([0-9]{1,5}$)/\1 \2 1/' >> MCcodon.dat
	line=$(grep 'tree length = ' results.mlc)
	echo $line | sed -E 's/^tree length = //' >> MCcodon.dat
	line=$(grep "$seqname:" results.mlc)
	echo $line >> MCcodon.dat
	line=$(grep "omega (dN/dS) = " results.mlc)
	echo $line | sed -E 's/^omega \(dN\/dS\) = //' >> MCcodon.dat
	line=$(grep "kappa (ts/tv) = " results.mlc)
	echo $line | sed -E 's/^kappa \(ts\/tv\) = //' >> MCcodon.dat
	echo "$(grep -A 16 'evolver' results.mlc | tail -n +2 )" >> MCcodon.dat
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



