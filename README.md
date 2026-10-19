# A resistance-robust suppression gene drive leveraging the accumulation of recessive deleterious mutations

This repository contains all code and data needed to reproduce results from "A resistance-robust suppression gene drive leveraging the accumulation of recessive deleterious mutations".

## SLiM

SLiM can be installed at https://messerlab.org/slim/. See the SLiM manual and workshop material for more details on using SLiM.

The SLiM file [WF_multilocus_cluster.slim](WF_multilocus_cluster.slim) was used for all simulations in the paper.

### Default parameters

* `make_csv=1`: this tells SLiM to create a csv file tracking allele frequencies each generation
* `population_size=1e4`: we always kept the population size fixed at 10,000
* `release_size=200`: this is the number of drive heterozygotes introduced at generation 100. Setting this to 200 ensures an initial drive allele frequency of 1%.
* `shared_cas9=0`: this is an option to simulate Cas9 saturation across the drive and target loci. This was turned off by default.
* `baseline_cleavage_rate=0.5`: this is the baseline germline cleavage rate at the target loci when there's no Cas9 saturation. 
  + note: we used 0.8 before; now switching to 0.5
* `saturation_factor=0.0`: this is a parameter controlling how much the target cleavage rate decline with `num_target_loci` due to Cas9 saturation. Setting this to 0.0 ensures there's no saturation effect.
* `drive_cleavage_rate=1.0`: we assumed perfect germline cleavage at the drive locus
* `func_resist_rate=0.01`: this is the resistance rate at the target loci (only used when `num_target_loci` > 0)
* `sd_broken_coeff=0`: we implemented an option to vary the fitness cost of target site disruption, drawing from a Gaussian distribution with this standard deviation. Setting this to 0 ensures the fitness cost is constant.
* `broken_dom=0.0`: this is the dominance coefficient of the disrupted target site. We always assumed recessive mutations and set this to 0.
* `recovery_val=0.2`: we tracked the time until the genetic load imposed by the drive fell below this value.
* `end_simulation_at_genetic_load_cutoff=0`: this is an option to end the simulation when the genetic load reached the state above. We didn't use this by default, since we wanted to observe later states of the population, particularly for the log-jam scenarios.
* `PLOT_HOMZ_COUNTS=F`: this is an option to liveplot the number of individuals with 0, 1, ..., `num_target_loci` homozygous disrupted target sites. It is helpful for observing log-jam scenarios but was not used by default.
* `PRINT_COMMON_HAPLOTYPES=F`: this is an option to print the most common haplotypes each generation to the SLiM console. This makes the simulation much slower, so I recommend turning this off and instead calling the `reportHaplotypeCounts` function from the Eidos console at desired time points.
* `NUM_HAPLOTYPES_PRINT=10`: when `PRINT_COMMON_HAPLOTYPES=T`, this parameter would control the number of haplotypes printed to the console. 

### Varied parameters

* `recomb_rate`: this is the recombination rate between target sites. We varied this across 0.5 (unlinked), 1e-2, 1e3, 1e-4, 1e-5, 1e-6, 1e-7
* `resistance_rate`: this is the functional resistance rate at the drive locus. We varied this for the single-locus section across 1e-5, 1e-4, and 1e-3. For the multilocus drives, we set this parameter to 1e-2.
* `drive_coeff`: this is the fitness cost of the drive allele. We set this to 1 or 0.2 for the single-locus section. For the multilocus drives, we set this parameter to 0.1
* `drive_dom`: this is the dominance coefficient at the drive locus. We set this to 0.0 for the single-locus section and 0.5 for the multilocus section.
* `num_target_loci`: this was varied throughout the paper. For the single-locus drives, we set this to 0 and for the multilocus drives, we explored 1,2,...,10.
* `broken_coeff`: this is the fitness cost of target site disruption. We varied this based on the `num_target_loci` and the target potential genetic load.

* `func_resist_rate`: for the multilocus unlinked drives, we will vary this in (0.01, 1/3)
* Cas9 saturation:
    + `saturation_factor`: exploring 0.0 (no saturation), 0.5 (some saturation), 1 (strong saturation)
    + `baseline_cleavage_rate`: 0.5 by default, with no saturation. When saturation is simulated, we set `baseline_cleavage_rate = 1`


## Data

* Data used in the single locus section are included in `data/single-locus.zip`
* Data used for the (s, n) panel are placed as follows:
  + Default parameters (50% flat cleavage rate at targets): 
      * Summary csv: `data/updated-s-vs-n-default.csv`
      * n = 10 trajectories (s = 0.1, 0.2, 0.3): `data/updated-s-vs-n-n10-trajectories.zip` 
      * s = 0.3 trajectories (n = 1, 5, 10): `data/updated-s-vs-n-s0.3-trajectories.zip`
      * n = 10, s = 0.25 trajectories: `data/updated-s-vs-n-n10-s0.25-trajectories.zip`
      * n = 1, s = 0.95 trajectories: `data/updated-s-vs-n-n1-s0.95-trajectories.zip`
  + r1 fraction of 1/3:
      * Summary csv: `updated-s-vs-n-r1-high.csv`
      * n = 10 trajectories (s = 0.1, 0.2, 0.3): `data/updated-s-vs-n-r1-high-n10-trajectories.zip`
      * s = 0.3 trajectories (n = 1, 5, 10): `data/updated-s-vs-n-r1-high-s0.3-trajectories.zip`
      * n = 10, s = 0.25 trajectories: `data/updated-s-vs-n-r1-high-n10-s0.25-trajectories.zip`
      * n = 1, s = 0.95 trajectories: `data/updated-s-vs-n-r1-high-n1-s0.95-trajectories.zip`
  + Cas9 saturated:
      * Summary csv: `updated-s-vs-n-Cas9-saturated.csv`
      * n = 10 trajectories (s = 0.1, 0.2, 0.3): `data/updated-s-vs-n-Cas9-saturated-n10-trajectories.zip`
      * s = 0.3 trajectories (n = 1, 5, 10): `data/updated-s-vs-n-Cas9-saturated-s0.3-trajectories.zip`
      * n = 10, s = 0.25 trajectories: `data/updated-s-vs-n-Cas9-saturated-n10-s0.25-trajectories.zip`
      * n = 1, s = 0.95 trajectories: `data/updated-s-vs-n-Cas9-saturated-n1-s0.95-trajectories.zip`
  + Cas9 saturated and r1 fraction of 1/3:
      * Summary csv: `updated-s-vs-n-Cas9-saturated-r1-high.csv`
      * n = 10 trajectories (s = 0.1, 0.2, 0.3): `data/updated-s-vs-n-Cas9-saturated-r1-high-n10-trajectories.zip`
      * s = 0.3 trajectories (n = 1, 5, 10): `data/updated-s-vs-n-Cas9-saturated-r1-high-s0.3-trajectories.zip`
      * n = 10, s = 0.25 trajectories: `data/updated-s-vs-n-Cas9-saturated-r1-high-n10-s0.25-trajectories.zip`
      * n = 1, s = 0.95 trajectories: `data/updated-s-vs-n-Cas9-saturated-r1-high-n1-s0.95-trajectories.zip`

* Data used for the varied recombination rates (s,n) heatmaps are placed as follows:
  + Recombination rate = 10^-7: `data/recomb10_7-s-vs-n.csv`
  + Recombination rate = 10^-5: `data/recomb10_5-s-vs-n.csv`
  + Recombination rate = 10^-4: `data/recomb10_4-s-vs-n.csv`
  + Recombination rate = 0.5: `data/updated-s-vs-n-default.csv`

## Figures

[figure-making.Rmd](figure-making.Rmd) shows how each plot in the paper was made. 

Plots were arranged in Adobe Illustrator 2025.
