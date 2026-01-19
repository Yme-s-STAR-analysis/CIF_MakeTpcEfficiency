#!/bin/bash

if [ $# -ne 2 ]; then
    echo "[WARNING] - $0: Usage: $0 <energy_tag> <action>"
    echo "==> An Example: $0 7 show"
    exit 1
fi

energy=$1
action=$2
xml_file1="Xsubmit_${energy}_14.xml"
xml_file2="Xsubmit_${energy}_15.xml"

if [ ! -f "$xml_file1" ]; then
    echo "[WARNING] - $0: Proton file $xml_file1 not found!"
    exit 1
fi
if [ ! -f "$xml_file2" ]; then
    echo "[WARNING] - $0: Antiproton file $xml_file2 not found!"
    exit 1
fi

target_path1=$(grep -o 'file:[^"]*/out/' "$xml_file1" | sed 's|file:\(.*\)/out/|\1|')
target_path2=$(grep -o 'file:[^"]*/out/' "$xml_file2" | sed 's|file:\(.*\)/out/|\1|')

if [ -z "$target_path1" ]; then
    echo "[WARNING] - $0: No valid path pattern found in $xml_file1"
    exit 1
fi
if [ -z "$target_path2" ]; then
    echo "[WARNING] - $0: No valid path pattern found in $xml_file2"
    exit 1
fi

case $action in
    show)
        echo "[LOG] - $0: proton output path: $target_path1"
        echo "[LOG] - $0: antiproton output path: $target_path2"
        ;;
    merge)
        output_file="${energy}GeV.hadd.root"
        echo "[LOG]- $0: merging files to $output_file"
        hadd "$output_file" "$target_path1"/out/*root "$target_path2"/out/*root
        ;;
    fmerge)
        output_file="${energy}GeV.hadd.root"
        echo "[LOG]- $0: merging files to $output_file"
        hadd -f "$output_file" "$target_path1"/out/*root "$target_path2"/out/*root
        ;;
    *)
        echo "[WARNING]: invalid action ('$action'), please try: show, merge, fmerge"
        exit 1
        ;;
esac

