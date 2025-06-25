
# conda source 
conda deactivate
conda activate microbetag


types=("treatment" "day")
# TYPE=$1

# if [[ "$1" != "day" && "$1" != "treatment" ]]; then
#     echo "Invalid argument: $1"
#     echo "Please use 'day' or 'treattmentt'."
#     exit 1
# fi

ko_merged_gz="mgg_prec/KEGG_annotations/ko_merged.txt.gz"


types=("treatment" "day")


run_microbetag_for_type() {
    local TYPE=$1
    local ko_merged_gz="mgg_prec/KEGG_annotations/ko_merged.txt.gz"  # ko_merged_gz="ko_merged.txt.gz"

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

            if [ -f "$ko_merged_gz" ]; then
                echo -e "Decompressing ko_merged.txt file.."
                gunzip "$ko_merged_gz"
            fi

            echo -e "Run microbetag"
            microbetag --config "$config_file" > log 2>&1

            if compgen -G "mgg_prec/mtag_net_*" > /dev/null; then
                rm -f mgg_prec/network_output.edgelist
                mv mgg_prec/mtag_net_* "microbetag_nets/${CASE}.cx2"
                echo -e "\n----  microbetag for $TYPE $CASE is complete. ---  \n"
            else
                echo -e "\n!! Something went wrong, no microbetag network was built for $TYPE $CASE\n"
            fi
        fi
    done
}


for type in day treatment; do
    run_microbetag_for_type "$type"
done



# for config_file in *.yaml; do

#     base=$(basename "$config_file")
#     part="${base#*_}"
#     CASE="${part%%.*}"

#     if [ -f "microbetag_nets/${TYPE}_${CASE}.cx2" ]; then
#         echo -e "A microbetag annotated network is already built for $TYPE $CASE."

#     else
#         echo -e "Building microbetag-annotated network for $TYPE $CASE"

#         # microbetag compresses ko_merged.txt always, so we make sure it's decompressed
#         if [ -f $ko_merged_gz ]; then
#             echo -e "Decompressing ko_merged.txt file.."
#             gunzip $ko_merged_gz
#         fi

#         # # Run microbetag
#         echo -e "Run microbetag"
#         microbetag --config "$config_file" > log 2>&1

#         # Remove the FlashWeave network, otherwise microbetag will use the previous one and you'll end up with the same network across all your days/treatments.
#         if compgen -G "mgg_prec/mtag_net_*" > /dev/null; then
    
#             rm mgg_prec/network_output.edgelist
#             mv mgg_prec/mtag_net_* microbetag_nets/${TYPE}_${CASE}.cx2
#             echo -e "\n----  microbetag for $TYPE $CASE is complete. ---  \n"
#         else
#             echo -e "\n!! Something went wrong, no microbetag network was built for $TYPE $CASE\n"
#         fi
#     fi
# done
