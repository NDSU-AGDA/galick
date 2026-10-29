# diagnostic reports have stable compact and detailed renderings

    Code
      print(x)
    Message
      
      -- galick_result: NOT OK -------------------------------------------------------
      Validation entries: 3 errors; 2 warnings.
      Entry counts are not counts of affected records or genotype calls.
      Full scan: 4 of 4 records; 5 samples.
      
      -- Impact in scanned records 
      Records with located findings: 3; with errors: 3.
      File/sample or unlocated entries: 1.
      Flagged genotype cells: 4; affected sample columns: 3.
      Cell counts cover genotype/FORMAT findings, not every call in a structurally
      flagged record.
      
      -- Rule impact 
    Output
      # A tibble: 5 x 6
        severity code                     n_entries n_records n_samples scope         
        <chr>    <chr>                        <int>     <int>     <int> <chr>         
      1 error    vcf_duplicated_chrom_pos         1         2        NA record        
      2 error    vcf_invalid_ref_alt              1         3        NA record        
      3 error    vcf_multiallelic                 1         3        NA record        
      4 warning  vcf_unsorted_records             1         2        NA record        
      5 warning  vcf_beagle_not_bgzf              1        NA        NA file_or_unloc~
    Message
      
      -- Overlapping record findings 
    Output
      # A tibble: 3 x 2
        n_records codes                                                               
            <int> <chr>                                                               
      1         1 vcf_duplicated_chrom_pos; vcf_invalid_ref_alt; vcf_multiallelic     
      2         1 vcf_duplicated_chrom_pos; vcf_invalid_ref_alt; vcf_multiallelic; vc~
      3         1 vcf_invalid_ref_alt; vcf_multiallelic; vcf_unsorted_records         
    Message
      
      -- Unsupported allele dependencies 
      Allele "-": 3 records; 4 dependent genotype cells.
      Example allele indices (record 1): 0=A, 1=G, 2=-
      ! Unsupported allele meaning is unresolved. Numeric FORMAT summaries are descriptive, not probabilities, universal quality thresholds, or a missing-call/deletion diagnosis.
      Removing an allele requires review of dependent calls and allele-numbered
      INFO/FORMAT annotations.
    Output
      # A tibble: 3 x 4
        genotype category           n_records n_genotype_cells
        <chr>    <chr>                  <int>            <dbl>
      1 2/2      unsupported_allele         2                2
      2 0/2      unsupported_allele         1                1
      3 2|.      unsupported_allele         1                1
    Message
      
      -- Numeric FORMAT evidence on problem records 
    Output
      # A tibble: 6 x 7
        field category           n_finite median n_missing n_nonfinite n_nonnumeric
        <chr> <chr>                 <int>  <dbl>     <int>       <int>        <int>
      1 Depth unsupported_allele        3   2            0           0            1
      2 Score unsupported_allele        2   0.75         0           2            0
      3 Depth missing                   0  NA            2           0            0
      4 Score missing                   0  NA            1           0            1
      5 Depth called_supported          9  20            0           0            0
      6 Score called_supported          9   0.9          0           0            0
    Message
      These pooled summaries are descriptive; within-record comparisons are reported
      separately.
      
      -- Action assessment 
      vcf_duplicated_chrom_pos [review_required]: Compare alleles and calls before
      choosing, merging or discarding duplicate records.
      vcf_invalid_ref_alt [review_required]: Confirm allele semantics. Review
      genotype indices and allele-numbered annotations before changing alleles.
      vcf_multiallelic [review_required]: Review target-specific representation;
      splitting or filtering requires genotype and annotation handling.

---

    Code
      summary(x)
    Message
      
      -- galick_result: NOT OK -------------------------------------------------------
      Validation entries: 3 errors; 2 warnings.
      Entry counts are not counts of affected records or genotype calls.
      Full scan: 4 of 4 records; 5 samples.
      
      -- Impact in scanned records 
      Records with located findings: 3; with errors: 3.
      File/sample or unlocated entries: 1.
      Flagged genotype cells: 4; affected sample columns: 3.
      Cell counts cover genotype/FORMAT findings, not every call in a structurally
      flagged record.
      
      -- Rule impact 
    Output
      # A tibble: 5 x 6
        severity code                     n_entries n_records n_samples scope         
        <chr>    <chr>                        <int>     <int>     <int> <chr>         
      1 error    vcf_duplicated_chrom_pos         1         2        NA record        
      2 error    vcf_invalid_ref_alt              1         3        NA record        
      3 error    vcf_multiallelic                 1         3        NA record        
      4 warning  vcf_unsorted_records             1         2        NA record        
      5 warning  vcf_beagle_not_bgzf              1        NA        NA file_or_unloc~
    Message
      
      -- Overlapping record findings 
    Output
      # A tibble: 3 x 2
        n_records codes                                                               
            <int> <chr>                                                               
      1         1 vcf_duplicated_chrom_pos; vcf_invalid_ref_alt; vcf_multiallelic     
      2         1 vcf_duplicated_chrom_pos; vcf_invalid_ref_alt; vcf_multiallelic; vc~
      3         1 vcf_invalid_ref_alt; vcf_multiallelic; vcf_unsorted_records         
    Message
      
      -- Unsupported allele dependencies 
      Allele "-": 3 records; 4 dependent genotype cells.
      Example allele indices (record 1): 0=A, 1=G, 2=-
      ! Unsupported allele meaning is unresolved. Numeric FORMAT summaries are descriptive, not probabilities, universal quality thresholds, or a missing-call/deletion diagnosis.
      Removing an allele requires review of dependent calls and allele-numbered
      INFO/FORMAT annotations.
    Output
      # A tibble: 3 x 4
        genotype category           n_records n_genotype_cells
        <chr>    <chr>                  <int>            <dbl>
      1 2/2      unsupported_allele         2                2
      2 0/2      unsupported_allele         1                1
      3 2|.      unsupported_allele         1                1
    Message
      
      -- Numeric FORMAT evidence on problem records 
    Output
      # A tibble: 6 x 7
        field category           n_finite median n_missing n_nonfinite n_nonnumeric
        <chr> <chr>                 <int>  <dbl>     <int>       <int>        <int>
      1 Depth unsupported_allele        3   2            0           0            1
      2 Score unsupported_allele        2   0.75         0           2            0
      3 Depth missing                   0  NA            2           0            0
      4 Score missing                   0  NA            1           0            1
      5 Depth called_supported          9  20            0           0            0
      6 Score called_supported          9   0.9          0           0            0
    Message
      These pooled summaries are descriptive; within-record comparisons are reported
      separately.
      
      -- FORMAT declarations 
    Output
      # A tibble: 5 x 5
        field number type    description                     profiled
        <chr> <chr>  <chr>   <chr>                           <lgl>   
      1 GT    1      String  Genotype                        FALSE   
      2 Score 1      Float   Custom score, not a probability TRUE    
      3 Depth 1      Integer Depth                           TRUE    
      4 AD    R      Integer Allele depths                   FALSE   
      5 Label 1      String  Text                            FALSE   
    Message
      Only unambiguous Number=1 Integer/Float declarations are profiled in extended
      mode.
      
      -- Within-record comparisons 
      Unsupported-allele calls versus fully called supported alleles at the same
      record.
    Output
      # A tibble: 2 x 6
        field n_comparable_records n_not_comparable_records n_lower n_equal n_higher
        <chr>                <int>                    <int>   <int>   <int>    <int>
      1 Score                    2                        0       1       0        1
      2 Depth                    2                        0       2       0        0
    Message
      Score: median difference 0; median ratio 1.028 from 2 positive baselines.
      Depth: median difference -13; median ratio 0.1361 from 2 positive baselines.
      Ratios use positive baseline medians only; numeric differences do not establish
      biological meaning.
      
      -- Genotype totals on problem records 
    Output
      # A tibble: 8 x 4
        genotype category           n_records n_genotype_cells
        <chr>    <chr>                  <int>            <dbl>
      1 2/2      unsupported_allele         2                2
      2 0/0      called_supported           2                6
      3 ./.      missing                    2                2
      4 0/2      unsupported_allele         1                1
      5 2|.      unsupported_allele         1                1
      6 1/1      called_supported           1                1
      7 0/1      called_supported           1                1
      8 1|1      called_supported           1                1
    Message
      
      -- Genotype categories on problem records 
    Output
      # A tibble: 11 x 4
         record genotype category               n
          <int> <chr>    <chr>              <int>
       1      1 2/2      unsupported_allele     1
       2      1 0/0      called_supported       1
       3      1 ./.      missing                1
       4      1 0/2      unsupported_allele     1
       5      1 2|.      unsupported_allele     1
       6      2 2/2      unsupported_allele     1
       7      2 1/1      called_supported       1
       8      2 ./.      missing                1
       9      2 0/1      called_supported       1
      10      2 1|1      called_supported       1
      11      3 0/0      called_supported       5
    Message
      
      -- Allele-index usage on problem records 
    Output
      # A tibble: 9 x 7
        record allele_index allele role  supported n_copies n_genotype_cells
         <int>        <int> <chr>  <chr> <lgl>        <dbl>            <int>
      1      1            0 A      REF   TRUE             3                2
      2      1            1 G      ALT   TRUE             0                0
      3      1            2 -      ALT   FALSE            4                3
      4      2            0 C      REF   TRUE             1                1
      5      2            1 T      ALT   TRUE             5                3
      6      2            2 -      ALT   FALSE            2                1
      7      3            0 G      REF   TRUE            10                5
      8      3            1 A      ALT   TRUE             0                0
      9      3            2 -      ALT   FALSE            0                0
    Message
      
      -- Record locations 
    Output
      # A tibble: 3 x 4
        record chrom pos   marker_id
         <int> <chr> <chr> <chr>    
      1      1 1     100   m1       
      2      2 1     100   m2       
      3      3 1     90    m3       
    Message
      
      -- Chromosomes 
    Output
      # A tibble: 2 x 3
        chrom n_records n_unsorted
        <chr>     <int>      <int>
      1 1             3          1
      2 2             1          0
    Message
      
      -- Action assessment 
      vcf_duplicated_chrom_pos [review_required]: Compare alleles and calls before
      choosing, merging or discarding duplicate records.
      vcf_invalid_ref_alt [review_required]: Confirm allele semantics. Review
      genotype indices and allele-numbered annotations before changing alleles.
      vcf_multiallelic [review_required]: Review target-specific representation;
      splitting or filtering requires genotype and annotation handling.
      vcf_unsorted_records [ordering_only]: Sort records after resolving invalid
      coordinates; do not change genotype calls.
      vcf_beagle_not_bgzf [storage_only]: Plain and compressed VCF are accepted as
      gt= or ref=; the runner does not require an input index.

