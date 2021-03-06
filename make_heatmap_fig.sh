#!/bin/bash

clear


SKIP=1

show_help()
{
echo "
        Usage: make_heatmap_fig.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
"
echo -e '\nmake_heatmap_fig.\n\nThis script Compiles the most important figures for visual inspection.\nGenerates a single .svg file containting the heatmap plots.\n\n'
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

cp f4.svg fig4.svg

ls -d */ | while read F; do 
    echo Processing "${F}"

    if [[ ! -f "./${F}/"delta-heatmap.svg ]] ; then
        echo ':('
        continue
    fi

    svg_stack.py --direction=h --margin=0 fig4.svg "./${F}"delta-heatmap.svg > fig4a.svg
	rm fig4.svg
    mv fig4a.svg fig4.svg
done

svg_stack.py --direction=h --margin=0 fig4.svg f4s.svg  > Heatmap.svg

rm fig4.svg 


cp f4b.svg fig4.svg

ls -d */ | while read F; do 
    echo Processing "${F}"

    if [[ ! -f "./${F}/"delta-heatmap2.svg ]] ; then
        echo ':('
        continue
    fi
    
	if [[ -d "./${F}/"LEP ]] ; then
    svg_stack.py --direction=h --margin=0 fig4.svg "./${F}"delta-heatmap2.svg > fig4a.svg
	rm fig4.svg
    mv fig4a.svg fig4.svg
    fi
done
svg_stack.py --direction=h --margin=0 fig4.svg f4s.svg  > Heatmap2.svg
rm fig4.svg 


cp f4c.svg fig4.svg

ls -d */ | while read F; do 
    echo Processing "${F}"

    if [[ ! -f "./${F}/"delta-heatmap3.svg ]] ; then
        echo ':('
        continue
    fi
    
	if [[ -d "./${F}/"LEP ]] ; then
    svg_stack.py --direction=h --margin=0 fig4.svg "./${F}"delta-heatmap3.svg > fig4a.svg
	rm fig4.svg
    mv fig4a.svg fig4.svg
    fi
done
svg_stack.py --direction=h --margin=0 fig4.svg f4s2.svg  > Heatmap3.svg
rm fig4.svg 



