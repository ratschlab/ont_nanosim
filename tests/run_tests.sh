#!/usr/bin/env bash

# Test that the NanoSim simualator runs without errors
set -eu

# todo: remove this file

mamba activate nanosim

set -x

reads_dir=example_nanosim_reads
reads_prefix="$reads_dir/reads"
mkdir -p "$reads_dir"

# checks fa.gz file
/usr/bin/time -v python "src/simulator.py" genome \
            --model_prefix "/cluster/customapps/biomed/grlab/users/mmordig/selective_sequencing/nanosim_trained_models/human_NA12878_DNA_FAB49712_guppy/training" \
            --ref_g "/cluster/work/grlab/mmordig/selective_sequencing/genomes/CHM13/chr1.fa.gz" \
            -dna_type linear \
            --output "$reads_prefix" \
            --number 100 \
            --seed 1 \
            --strandness 0.5 \
            --basecaller guppy \
            --aligned_rate "100%" \
            --num_threads 4 \
            --no_flanking \
            --no_error_profile

rm -rf "$reads_dir"