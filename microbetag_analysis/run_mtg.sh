#!/usr/bin/bash

eval "$(conda shell.bash hook)"
conda activate microbetag

types=("treatment" "day")
ko_merged_gz="mgg_prec/KEGG_annotations/ko_merged.txt.gz"


run_microbetag_for_type() {
    local TYPE=$1
    local ko_merged_gz="mgg_prec/KEGG_annotations/ko_merged.txt.gz"

    if [ "$TYPE" == "day" ]; then
        files=(config_day_*.yaml)
    else
        files=(config_TG_*.yaml)
    fi

    for config_file in "${files[@]}"; do
        [ -e "$config_file" ] || continue

        base=$(basename "$config_file")
        part="${base#*_}"
        CASE="${part%%.*}"

        echo "microbetag_nets/${CASE}.cx2"

        if [ -f "microbetag_nets/${CASE}.cx2" ]; then
            echo -e "A microbetag annotated network is already built for $TYPE $CASE."
        else
            echo -e "Building microbetag-annotated network for $TYPE $CASE"

            #if [ -f "$ko_merged_gz" ]; then
            #    echo -e "Decompressing ko_merged.txt file.."
            #    gunzip -k "$ko_merged_gz"   # keep gz file to avoid re-downloading
            #fi

            echo -e "Run microbetag"
            microbetag --config "$config_file" > log 2>&1

            if compgen -G "mgg_prec/mtag_net_*" > /dev/null; then
                rm -f mgg_prec/network_output.edgelist
                if ls mgg_prec/*.cyjs >/dev/null 2>&1; then
                    rm mgg_prec/*.cyjs
                fi
                mv mgg_prec/mtag_net_* "microbetag_nets/${CASE}.cx2"
                echo -e "\n----  microbetag for $TYPE $CASE is complete. ---  \n"
            else
                echo -e "\n!! Something went wrong, no microbetag network was built for $TYPE $CASE\n"
            fi
        fi

        # ensure clean environment between runs
        rm -rf __pycache__ .mypy_cache 2>/dev/null
        unset CASE base part config_file
    done
}



for type in day treatment; do
    run_microbetag_for_type "$type"
done
