#!/bin/bash
clear

show_help()
{
echo "
        Usage: bulk_rename.sh [-h]]

        -h Show this Help

"

echo -e '\nReads the first fasta ID and renames the alignment file.'

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


ls -d */ | while read D; do
	cd "${D}"
	ls * | while read F; do
		FIRST=$(sed -n '1p' "${F}")
		SEC=$(echo "${FIRST}" | sed -E 's/>([^ ]*).*/\1/')
		mv "${F}" "${SEC}".fas
	done
	cd ..
done



