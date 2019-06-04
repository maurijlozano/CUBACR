#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: EXPR_vs_CNSV.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"

echo -e '\nGenerates protein Expression vs percent of totally conserved amino acids.\n\n'

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

touch CNSV.table.temp
ls */*CE.table | while read F; do
	sed '1d' ${F} | sed -E 's/([^\/]*)\/[^\/]*\//\1 /' | sed -E 's/.stats//' >> CNSV.table.temp
done
mv CNSV.table.temp CNSV.table
Rscript --vanilla EXPR_vs_CNSV.R
 
