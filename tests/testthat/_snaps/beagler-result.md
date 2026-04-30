# print.beagler_result snapshot for ok result

    Code
      print(r)
    Message
      
      -- beagler_result: OK ----------------------------------------------------------
      Errors: 0; Warnings: 0
      
      -- Summary 
      n_samples: 100
      n_markers: 50

# print.beagler_result snapshot for failing result

    Code
      print(r)
    Message
      
      -- beagler_result: NOT OK ------------------------------------------------------
      Errors: 2; Warnings: 1
      
      -- Errors (top 3) 
      x vcf_unsorted: Records out of order
      x vcf_dup_id: Marker ID rs1 appears twice
      
      -- Warnings (top 3) 
      ! vcf_low_maf: MAF below 0.01 on 5 markers
      
      -- Summary 
      n_records: 1000

