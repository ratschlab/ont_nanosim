# Checking basic things about the simulation

We check that running it twice results in identical output (also up to interleaving of processes!).
## Metagenomic reads

```{bash}
rm -rf simulated_metagenomic_reads simulated_metagenomic_reads2

python /Users/maximilianmordig/Desktop/sequencing/NanoSim/src/simulator.py metagenome \
    -gl "metagenome_sim_configs/metagenome_list_for_simulation" -a "metagenome_sim_configs/abundance_for_simulation_multi_sample.tsv" \
    -dl "metagenome_sim_configs/dna_type_list.tsv" \
    -c "metagenome_ERR3152366_Log/training" \
    -o "simulated_metagenomic_reads/noflanking_new_reads" -max 1000 --seed 2 -b guppy --strandness 0.5 -t 4 --aligned_rate "100%"

mv simulated_metagenomic_reads simulated_metagenomic_reads2

python /Users/maximilianmordig/Desktop/sequencing/NanoSim/src/simulator.py metagenome \
    -gl "metagenome_sim_configs/metagenome_list_for_simulation" -a "metagenome_sim_configs/abundance_for_simulation_multi_sample.tsv" \
    -dl "metagenome_sim_configs/dna_type_list.tsv" \
    -c "metagenome_ERR3152366_Log/training" \
    -o "simulated_metagenomic_reads/noflanking_new_reads" -max 1000 --seed 2 -b guppy --strandness 0.5 -t 4 --aligned_rate "100%"
```

## Genomic reads

```{bash}
rm -rf simulated_reads simulated_reads2
python /Users/maximilianmordig/Desktop/sequencing/NanoSim/src/simulator.py genome \
    -rg "ref_data/GCF_000001405_normalized.40_GRCh38.p14_genomic_beginning_reduced.fna" \
    -c "human_NA12878_DNA_FAB49712_guppy/training" \
    -o "simulated_reads/noflanking_new_reads" -n 2 -max 1000 --seed 2 -b guppy -dna_type linear -t 4 --aligned_rate "100%"

mv simulated_reads simulated_reads2

python /Users/maximilianmordig/Desktop/sequencing/NanoSim/src/simulator.py genome \
    -rg "ref_data/GCF_000001405_normalized.40_GRCh38.p14_genomic_beginning_reduced.fna" \
    -c "human_NA12878_DNA_FAB49712_guppy/training" \
    -o "simulated_reads/noflanking_new_reads" -n 2 -max 1000 --seed 2 -b guppy -dna_type linear -t 4 --aligned_rate "100%"
```
