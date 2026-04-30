# validate_genotype_table snapshot for ok case

    Code
      print(result)
    Message
      
      -- beagler_result: OK ----------------------------------------------------------
      Errors: 0; Warnings: 0
      
      -- Summary 
      n_samples: 3
      n_markers: 2
      inferred_encoding: AHB
      inferred_ploidy: 2
      missing_rate_overall: 0
      missing_rate_by_sample: 0, 0, and 0
      missing_rate_by_marker: 0 and 0
      n_monomorphic_markers: 0
      n_duplicated_samples: 0
      n_duplicated_markers: 0

# validate_genotype_table snapshot for failing case

    Code
      print(result)
    Message
      
      -- beagler_result: NOT OK ------------------------------------------------------
      Errors: 1; Warnings: 1
      
      -- Errors (top 3) 
      x genotype_table_duplicated_samples: Sample ID 's1' appears more than once.
      
      -- Warnings (top 3) 
      ! genotype_table_monomorphic_markers: Marker 'rs2' is monomorphic.
      
      -- Summary 
      n_samples: 3
      n_markers: 2
      inferred_encoding: AHB
      inferred_ploidy: 2
      missing_rate_overall: 0
      missing_rate_by_sample: 0, 0, and 0
      missing_rate_by_marker: 0 and 0
      n_monomorphic_markers: 1
      n_duplicated_samples: 1
      n_duplicated_markers: 0

