#!/bin/bash
clear
echo -e "Runing CA"
ls -d */ | while read D; do
	cd ${D}
	cat $(ls -d */* | grep '[Ff][Aa][Ss]\?$' | grep -E 'CC/genoma|CC/GNM|CC/genome' | grep 'modal_freq.fas$' -v) | sed -e 's/[ \t_]/-/g' -e 's/\(>.*\)/\1_Gene/' > genes.fas
	cat genes.fas $(ls -d */*modal_freq.fas | grep -E 'CC/genoma|CC/GNM|CC/genome') $(ls -d */* -v | grep 'modal.freq.fas$' | grep '^CC/genoma' -v | grep -E '^CC/GNM' -v ) $(ls | grep 'CR.modal_freq.fas$') $(ls | grep 'VR.modal_freq.fas$') > coa.fas
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
	Rscript --vanilla graph_coa_final.R "./${D}" genes.coa &>/dev/null
	Rscript --vanilla graph_codons_final.R "./${D}" codon.coa &>/dev/null
	cd ${D}
	rm coa_var_axis1_2.txt
	cd ..
	echo -e "\n"
done
