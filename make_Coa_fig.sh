#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: make_Coa_fig.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
echo -e '\nThis script Compiles the most important figures for visual inspection.\nGenerates a single .svg file containting the CA plots.\n\n'
exit 1
}

while getopts yh option; do
    case "${option}" in
        y) SKIP=0;;
        h) show_help;;
    esac
done



if [ ${SKIP} != 1 ]
then
    echo -e 'Beginning calculations.\n'
else
    echo -e 'PRESS ANY KEY TO CONTINUE.\n'
    read -n 1 -s
fi


ls -d */ | while read F; do 
    echo Processing "${F}"
    if [[ ! -f "./${F}/"CA_genes.svg || ! -f "./${F}/"CA_codons.svg ]] ; then
        echo 'You need to run CA analysis first.'
        continue
    fi

    FE=$(echo "${F}" | sed -e 's/[_ ][_ ][_ ]/_/' | sed -e 's/[_ ][_ ]/_/')
    FS=$(echo "${FE}" | sed -e 's/\(^[A-Za-z0-9]*\)[_ ]\(.\)[A-Za-z0-9]*[_ ]\([A-Za-z0-9]*\).*/\1_\2.\3/')
    FN=$(echo "${FS}" | sed -e 's/\///')

    svg_stack.py --direction=h --margin=30 "./${F}/"CA_genes.svg "./${F}/"CA_codons.svg  > "${FN}"_COA.svg
    
done

