## Default Wright–Fisher Parameters

### Directory for Data Storage
- Manually create your own folder somewhere, and then set the file path to that folder	

### General Constants
- "population_size” : 1e4
- "release_size” : 100
- “recomb_rate” : 1e-6  (rate of recombination between target sites)

### Constants for Drive Locus
- “resistance_rate” : 0.01  (all functional)
- "conversion_rate” : 1  (aka cleavage)
- "drive_coeff” : -0.1  (selection coefficient)
- "drive_dom” : 0.5  (dominance)

### Constants for Distant Target Loci
- "num_target_loci” : 20  (number of distant target loci, whether linked or unlinked)
- "disruption_rate” : 0.2  (rate of non-functional target loci, aka broken alleles)
- "func_resist_rate” : 0.01  (rate of functional resistance)
- "broken_coeff” : -0.15  (selection coefficient)
