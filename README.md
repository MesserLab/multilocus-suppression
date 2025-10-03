## Default Wright–Fisher Parameters

### Directory for Data Storage
- Manually create your own folder somewhere, and then set the file path to that folder	

### General Constants
- "population_size” : 1e4
- "release_size” : 100  (in addition to the population size)
- “recomb_rate” : 1e-5  (rate of recombination between target sites, drive is unlinked)

### Constants for Drive Locus
- “resistance_rate” : 0.001  (all functional alleles, no non-functional)
- "cleavage_rate” : 0.90  (cleavage rate before resistance then conversion)
- "drive_coeff” : -0.1  (selection coefficient)
- "drive_dom” : 0.5  (dominance)

### Constants for Distant Target Loci
- "num_target_loci” : 20  (number of distant target loci, whether linked or unlinked)
- "disruption_rate” : 0.2  (cleavage rate of wild-type, creates non-functional/broken alleles)
- "func_resist_rate” : 0.001  (rate of functional resistance, on top of disruption rate)
- "broken_coeff” : -0.15  (selection coefficient)
- "sd_broken_coeff" : 0  (standard deviation, normally distributed around broken_coeff)
- "broken_dom" : 0  (dominance)
