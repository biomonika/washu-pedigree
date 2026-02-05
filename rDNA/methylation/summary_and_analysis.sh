#!/bin/bash

SCRIPT_DIR=$1
out_dir=$2
doing_subregion=$3
doing_tr=$4

rm -f $out_dir/final_analysis/*.png

bash $SCRIPT_DIR/get_read_summary.sh $SCRIPT_DIR $out_dir $doing_subregion $doing_tr

# if [[ $# -gt 2 ]]; then
# 	bash $SCRIPT_DIR/get_read_summary.sh $SCRIPT_DIR $out_dir 1
# 	python $SCRIPT_DIR/run_ordering_analysis.py $out_dir/final_analysis/read_summary.txt > $out_dir/final_analysis/ordering_analysis.txt
# 	python $SCRIPT_DIR/violin_plot.py $out_dir/final_analysis 1
# else
# 	bash $SCRIPT_DIR/get_group_summary.sh $SCRIPT_DIR $out_dir
# fi
