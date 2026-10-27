# validate_marker_map snapshot for ok case

    Code
      print(result)
    Message
      
      -- galick_result: OK -----------------------------------------------------------
      Errors: 0; Warnings: 0
      
      -- Summary 
      n_markers_in_map: 3
      inferred_chrom_style: chr_prefixed
      chromosome_summary: c("chr1", "chr2"), 2:1, and c(0, 0)
      n_duplicated_marker_ids: 0
      n_duplicated_chrom_pos: 0
      n_unsorted_within_chrom: 0
      markers_in_table_not_in_map:
      markers_in_map_not_in_table:
      n_markers_in_table:

# validate_marker_map snapshot for failing case

    Code
      print(result)
    Message
      
      -- galick_result: NOT OK -------------------------------------------------------
      Errors: 4; Warnings: 0
      
      -- Errors (top 3) 
      x marker_map_invalid_marker_id: 1 marker_id value(s) are NA, empty, or whitespace-only.
      x marker_map_invalid_pos: 1 pos value(s) are NA, non-finite, non-positive, or non-integer.
      x marker_map_duplicated_marker_ids: marker_id 'rs1' appears more than once.
      i ... and 1 more.
      
      -- Summary 
      n_markers_in_map: 3
      inferred_chrom_style: chr_prefixed
      chromosome_summary: chr1, 2, and 0
      n_duplicated_marker_ids: 1
      n_duplicated_chrom_pos: 1
      n_unsorted_within_chrom: 0
      markers_in_table_not_in_map:
      markers_in_map_not_in_table:
      n_markers_in_table:

