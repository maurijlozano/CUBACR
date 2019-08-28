#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: Nc_table.sh [-h]

        -h Show this Help
"
echo -e '\nGenerates modal sequences Nc report table.\n\n'
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
    if [ ! -f NC_GC3s.tab ]
    then
    	echo -e "Please run Nc_plot.sh first.\n"
    	exit
    fi
     
    cd ..
    Rscript --vanilla Nc_table.R "./${D}"
  
done

cat $(ls */Nc_means.txt) | sed -e 's/.SP.*//' -e 's/"//g' -e 's/\.\///' | cut -d' ' -f2-6 | sed '/^$/d' > tmp.txt
echo "SP title Nc SD SEM" > Nc_means_table.txt
cat tmp.txt >> Nc_means_table.txt 
rm tmp.txt
