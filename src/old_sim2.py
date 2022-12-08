# copied from old_simulator.py

def simulation_aligned_metagenome(min_l, max_l, median_l, sd_l, out_reads, out_error, kmer_bias, basecaller,
                                  read_type, fastq, num_simulate, per=False, chimeric=False):
    
    # Simulate aligned reads
    out_reads = open(out_reads, "w")
    out_error = open(out_error, "w")

    id_begin = '@' if fastq else '>'

    remaining_reads = num_simulate
    if chimeric:
        num_segment = np.random.geometric(1/segment_mean, num_simulate)
    else:
        num_segment = np.ones(num_simulate, dtype=int)
    remaining_segments = num_segment
    remaining_gaps = remaining_segments - 1
    passed = 0
    current_species_bases = {species: 0 for species in dict_abun.keys()} #todo10
    while remaining_reads > 0:
        if per:
            ref_lengths = get_length_kde(kde_aligned, sum(remaining_segments)) if median_l is None else \
                np.random.lognormal(np.log(median_l), sd_l, remaining_segments)
            ref_lengths = [x for x in ref_lengths if min_l <= x <= max_l]
        else:
            remainder_lengths = get_length_kde(kde_ht, int(remaining_reads * 1.3), True)
            remainder_lengths = [x for x in remainder_lengths if x >= 0]
            head_vs_ht_ratio_list = get_length_kde(kde_ht_ratio, int(remaining_reads * 1.5))
            head_vs_ht_ratio_list = [x for x in head_vs_ht_ratio_list if 0 <= x <= 1]
            if median_l is None:
                ref_lengths = get_length_kde(kde_aligned, sum(remaining_segments))
            else:
                total_lengths = np.random.lognormal(np.log(median_l + sd_l ** 2 / 2), sd_l, remaining_reads)
                num_current_loop = min(remaining_reads, len(remainder_lengths), len(head_vs_ht_ratio_list))
                ref_lengths = total_lengths[:num_current_loop] - remainder_lengths[:num_current_loop]
            ref_lengths = [x for x in ref_lengths if 0 < x <= max_l]

        gap_lengths = get_length_kde(kde_gap, sum(remaining_gaps), True) if sum(remaining_gaps) > 0 else []
        gap_lengths = [max(0, int(x)) for x in gap_lengths]

        # Select strain/species to simulate
        species_pool, ref_lengths, remaining_segments = \
            assign_species(ref_lengths, remaining_segments, current_species_bases)

        is_reversed = random.random() > strandness_rate #todo10: why here??

        seg_pointer = 0
        gap_pointer = 0
        species_pointer = 0
        for each_read in xrange(len(remaining_segments)):
            
            segments = remaining_segments[each_read]
            # In case too many ref length was filtered previously
            if (not per and each_read >= min(len(head_vs_ht_ratio_list), len(remainder_lengths))) or \
                seg_pointer + segments > len(ref_lengths):
                break
            ref_length_list = [int(round(ref_lengths[seg_pointer + x])) for x in range(segments)]
            gap_length_list = [int(round(gap_lengths[gap_pointer + x])) for x in range(segments - 1)]
            species_list = [species_pool[species_pointer + x] for x in range(segments)]


            if per:
                seg_pointer += 1
                gap_pointer += 1
                species_pointer += 1
                with total_simulated.get_lock():
                    sequence_index = total_simulated.value
                    total_simulated.value += 1

                # Extract middle region from reference genome
                new_read = ""
                new_read_name = ""
                base_quals = []
                for seg_idx in range(len(ref_length_list)):
                    new_seg, new_seg_name = extract_read("metagenome", ref_length_list[seg_idx], species_list[seg_idx])
                    new_read += new_seg
                    new_read_name += new_seg_name
                    if fastq:
                        base_quals.extend(mm.trunc_lognorm_rvs("match", read_type, basecaller,
                                                               ref_length_list[seg_idx]).tolist())

                new_read_name = new_read_name + "_perfect_" + str(sequence_index)
                read_mutated = case_convert(new_read)  # not mutated actually, just to be consistent with per == False

                head = 0
                tail = 0

                if is_reversed:
                    new_read_name += "_R"
                else:
                    new_read_name += "_F"
                
                new_read_name += "_0_" + str(sum(ref_length_list)) + "_0"

            else:
                gap_list = []
                gap_base_qual_list = []
                seg_length_list = []
                seg_error_dict_list = []
                seg_error_count_list = []
                remainder = int(round(remainder_lengths[each_read]))
                head_vs_ht_ratio = head_vs_ht_ratio_list[each_read]

                total = remainder
                for each_gap in gap_length_list:
                    mutated_gap, gap_base_quals = simulation_gap(each_gap, basecaller, read_type, "metagenome", fastq)
                    gap_list.append(mutated_gap)
                    gap_base_qual_list.append(gap_base_quals)
                    total += len(mutated_gap)
                for each_ref in ref_length_list:
                    middle, middle_ref, error_dict, error_count = \
                        error_list(each_ref, match_markov_model, match_ht_list, error_par, trans_error_pr, fastq)
                    total += middle
                    seg_length_list.append(middle_ref)
                    seg_error_dict_list.append(error_dict)
                    seg_error_count_list.append(error_count)

                if total < min_l or total > max_l:
                    continue

                seg_pointer += segments
                gap_pointer += segments - 1
                species_pointer += segments

                with total_simulated.get_lock():
                    sequence_index = total_simulated.value
                    total_simulated.value += 1

                if remainder == 0:
                    head = 0
                    tail = 0
                else:
                    head = int(round(remainder * head_vs_ht_ratio))
                    tail = remainder - head

                # Extract middle region from reference genome
                num_seg = len(seg_length_list)
                new_seg_list = [None] * num_seg
                read_name_components = list()
                for seg_idx in range(num_seg):
                    new_seg_list[seg_idx], new_seg_name = extract_read("metagenome", seg_length_list[seg_idx], species_list[seg_idx])
                    read_name_components.append(new_seg_name)
                    if seg_idx < len(gap_list):
                        read_name_components.append("gap_" + str(len(gap_list[seg_idx])))
                new_read_name = ';'.join(read_name_components) + "_aligned_" + str(sequence_index)
                
                if num_seg > 1:
                    new_read_name += "_chimeric"
                
                if is_reversed:
                    new_read_name += "_R"
                else:
                    new_read_name += "_F"
                
                new_read_name += "_" + str(head) + \
                                 "_" + ";".join(str(x) for x in seg_length_list) + \
                                 "_" + str(tail)
                    
                read_mutated = ""
                base_quals = []
                for seg_idx in range(num_seg):
                    # Mutate read
                    new_seg = case_convert(new_seg_list[seg_idx])
                    seg_mutated, seg_base_quals = \
                        mutate_read(new_seg, new_read_name, out_error, seg_error_dict_list[seg_idx],
                                    seg_error_count_list[seg_idx], basecaller, read_type, fastq, kmer_bias)

                    if kmer_bias:
                        seg_mutated, seg_base_quals = mutate_homo(seg_mutated, seg_base_quals, kmer_bias, basecaller,
                                                                  None)
                    read_mutated += seg_mutated
                    base_quals.extend(seg_base_quals)
                    if seg_idx < len(gap_list):
                        read_mutated += gap_list[seg_idx]
                        base_quals.extend(gap_base_qual_list[seg_idx])

                    # Update base level abundance info
                    current_species_bases[species_list[seg_idx]] += len(new_seg)

                if fastq:  # Get head/tail qualities and add to base_quals
                    ht_quals = mm.trunc_lognorm_rvs("ht", read_type, basecaller, head + tail).tolist()
                    base_quals = ht_quals[:head] + base_quals + ht_quals[head:]

            # Add head and tail region
            read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + \
                           ''.join(np.random.choice(BASES, tail))

            # Reverse complement half of the reads
            if is_reversed:
                read_mutated = reverse_complement(read_mutated)
                base_quals.reverse()

            out_reads.write(id_begin + new_read_name + '\n')
            out_reads.write(read_mutated + '\n')

            if fastq:
                out_reads.write("+\n")
                out_quals = "".join([chr(qual + 33) for qual in base_quals])
                out_reads.write(out_quals + "\n")

            check_print_progress(sequence_index)

            passed += 1

        remaining_reads = num_simulate - passed
        remaining_segments = num_segment[passed:]
        remaining_gaps = remaining_segments - 1

    sys.stdout.write('\n')
    out_reads.close()
    out_error.close()