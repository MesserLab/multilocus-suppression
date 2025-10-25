# Plan for figures

## Figure 1: Schematic diagram (Isabel)

* Compare the standard dsx drive to the multilocus drive with unlinked target sites

## Part 1: WF model

### Figure 2: Dynamics of standard homing drives with no distant-site targets (Rachel)

* First row: dsx drive (s = 1, h = 0, cleavage rate = 1, 20 replicates) varying the r1 rate
  + Left plot: mean fitness over time
      * Plot the mean trajectory and report the fraction of replicates with suppression (for the lowest r1 rate)
  + Right plot: corresponding allele frequencies

* Second row: homing drive (h = 0.5, cleavage rate = 1, r1 rate = 0.001, 20 replicates) varying s
  + Left plot: mean fitness over time
  + Right plot: corresponding allele frequencies
  
* Takeaway from this figure: dsx can impose a high genetic load (minimum fitness) on the population but once resistance arises, is rapidly outcompeted. A homing drive with a lower fitness cost (s) can persist in the population longer after resistance arises, but it can't impose as high of a genetic load (minimum fitness)

### Figure 3: Introduction to the multilocus system (Rachel)

* Parameters:
  + No Cas9 sharing between drive and targets
  + 20 replicates per parameter combo
  + Drive:
      * cleavage rate = 1
      * h = 0.5
      * s = 0.1
      * r1 rate = 0.001
  + Target sites:
      * num target sites = 5
      * fitness cost at sites = 0.2
      * recessive (h = 0)
      * baseline cleavage rate = 1
      * saturation factor = 0.5
      * r1 rate = 0.01
      * recombination rate = 0.5

* Show the mean fitness over time and drive and disrupted-target frequencies (in same plot or two subplots)

* Takeaway from this figure: decoupling the drive and the target can lead to a high genetic load (minimum fitness) and longer time until population rebound (back to fitness 1)

### Figure 4: Exploring different number of target sites and cumulative fitness costs (Isabel)

* Parameters:
  + No Cas9 sharing between drive and targets
  + 20 replicates per parameter combo
  + Drive:
      * cleavage rate = 1
      * h = 0.5
      * s = 0.1
      * r1 rate = 0.001
  + Target sites:
      * num target sites varies in 1,2,...,10
      * fitness cost at sites varies such that (1-s)^n is some constant
      * recessive (h = 0)
      * baseline cleavage rate = 1
      * saturation factor = 0.5
      * r1 rate = 0.01
      * recombination rate = 0.5
      * Constant (1-s)^n = 0.01, 0.1, 0.25, 0.5

* Plot:
  + x-axis = number of target sites
  + y-axis = time until all broken alleles are lost 
  + Hold (1-s)^n = constant. For higher n, the time to recovery increases -- it takes longer for all broken alleles to be lost (Cas9 saturation will only slow the time to minimum fitness; can show that in supplemental plot)
    * Different line types for different cumulative costs. If the cost is higher, the time to recovery will be faster, since r1 alleles have a greater selective advantage.


* In the supplement, can show the same plot for different saturations (drive site also saturated; no saturation)

* Transition: how can we further increase the time to recovery? By increasing the amount of linkage between the target sites.

### Figure 5: varying the target site recombination rate (Rachel)

* Parameters:
  + No Cas9 sharing between drive and targets
  + 20 replicates per parameter combo
  + Drive:
      * cleavage rate = 1
      * h = 0.5
      * s = 0.1
      * r1 rate = 0.001
  + Target sites:
      * num target sites is set at 10
      * fitness cost at each site is set at 0.2
      * recessive (h = 0)
      * baseline cleavage rate = 1
      * saturation factor = 0.5
      * r1 rate = 0.01
      * recombination rate varies in 0.5, 1e-3, 1e-5, 1e-7, 0

* Plot shows the mean fitness over time at different recombination rates
  + Zoom in on a plateaued line at the end:
      * Show a diagram of fixed haplotypes there (ex: r1-broken-r1-broken / broken-r1-broken-r1)

* Takeaway from this figure:
  + If the target sites are linked, then assortative overdominance (the log-jam) can allow broken alleles to stay in the population for longer -> drive imposes a genetic load on the population for longer
  + These log-jam individuals will produce fewer offspring (if targets affect viability) or many offspring with low fecundity

#### Supplemental figure comparing our system to the Faber system

* In the Faber system, there's just 1 target site (unlinked from drive) with a fitness cost of 1 (2 disrupted alleles = sterility)

* Parameters:
  + No Cas9 sharing between drive and targets
  + Drive:
      * cleavage rate = 1 (maybe lower?)
      * h = 0.5
      * s = 0.1
      * r1 rate = 0.001 (maybe lower?)
  + Target sites:
      * num target sites = 1
      * fitness cost at each site = 1
      * recessive (h = 0)
      * baseline cleavage rate = 1 (maybe lower?)
      * r1 rate = 0.01 (maybe lower?)

* Plot shows the mean fitness over time

* Takeaway from this figure: this drive can impose a high genetic load but is also vulnerable to resistance

## Part 2: nonWF panmictic model

* More realistic systems
* Explore haplolethal or haplosufficient homing-rescue drive (at the drive locus)
* Maybe compare 2 of our multilocus systems (recomb rate of 0.5 and recomb rate much lower?), the dsx drive, and the Faber system:
  + Suppression rate?
  + Time to suppression?
  + Effect of different density curves?
