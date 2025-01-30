

# `microbetag` analysis 

In this directory, we share:

- `sample_metadata.xlsx`: metabolites of interest and samples metadata in original format
- `split_data_per_day.R`: script to split data and metadata per day in a microbetag-friendly format
- `microbetag_input`: the per-day and the overall abundance and metadata files provided as input. Also, the basic `config.yml` file (only the abundance and metadata files change across the runs)
- `microbetag_nets`: the per-day and overall microbetag annotated networks (`.cx2`) files
- `per_day_cx_stats.ipynb`: a notebook with some first statistics on the annotated networks


The rows of `sample` and `individual` were removed from the per-day metadata files, before using them on microbetag.

