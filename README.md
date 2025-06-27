## Default Wright–Fisher Parameters

### Directory for Data Storage
- Manually create your own folder somewhere, and then set the file path to that folder	

### General Constants
- "population_size” : 1e4
- "release_size” : 100  (in addition to the population size)
- “recomb_rate” : 1e-6  (rate of recombination between target sites, drive is unlinked)

### Constants for Drive Locus
- “resistance_rate” : 0.01  (all functional alleles, no non-functional)
- "conversion_rate” : 1  (aka cleavage rate before resistance then conversion)
- "drive_coeff” : -0.1  (selection coefficient)
- "drive_dom” : 0.5  (dominance)

### Constants for Distant Target Loci
- "num_target_loci” : 20  (number of distant target loci, whether linked or unlinked)
- "disruption_rate” : 0.2  (cleavage rate of wild-type, creates non-functional/broken alleles)
- "func_resist_rate” : 0.01  (rate of functional resistance, on top of disruption rate)
- "broken_coeff” : -0.15  (selection coefficient, all fully recessive)
