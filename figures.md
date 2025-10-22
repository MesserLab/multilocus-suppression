# Plan for figures

## Figure 1: Schematic diagram (Isabel)

* Compare the standard dsx drive to the multilocus drive with unlinked target sites

## Part 1: WF model

### Figure 2: Dynamics of standard homing drives with no distant-site targets (Rachel)

* First row: dsx drive (s = 1, h = 0, cleavage rate = 1, 20 replicates) varying the r1 rate
  + Left plot: mean fitness over time
      * Plot the mean trajectory and report the fraction of replicates with suppression (for the lowest r1 rate)
  + Right plot: corresponding allele frequencies

* Second row: homing drive (h = 0.5, cleavage rate = 1, r1 rate = ?, 20 replicates) varying s
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
      * num target sites = ? (maybe 5)
      * fitness cost at sites = ? (maybe 0.25)
      * recessive (h = 0)
      * baseline cleavage rate = 1
      * saturation factor = 0.5
      * r1 rate = 0.01
      * recombination rate = 0.5

* Show the mean fitness over time and drive and disrupted-target frequencies (in same plot or two subplots)

* Takeaway from this figure: decoupling the drive and the target can lead to a high genetic load (minimum fitness) and longer time until population rebound (back to fitness 1)

### Figure 4: Exploring different number of target sites and fitness cost at each site (Isabel)

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
      * fitness cost at sites varies in 0.1, 0.2, ..., 1.0
      * recessive (h = 0)
      * baseline cleavage rate = 1
      * saturation factor = 0.5
      * r1 rate = 0.01
      * recombination rate = 0.5

* 2 heatmaps:
  1. Color shows minimum population fitness (aka max genetic load) 
  2. Color shows time-to-rebound
  3. Could combine the above and show the average area *above* the curve (between fitness = 1 and fitness trajectory)

* Takeaway from this figure:
  + More target sites and a larger fitness cost at each -> lowest population fitness but fastest time-to-recovery
  + Best to target a lot of sites with a more moderate cost -> slower suppression but also slower time-to-recovery

* Transition: how can we allow the multilocus system to remain in the population for longer? Vary linkage.

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
      * num target sites is set at ? (maybe 5 or 10)
      * fitness cost at each site is set at ?
      * recessive (h = 0)
      * baseline cleavage rate = 1
      * saturation factor = 0.5
      * r1 rate = 0.01
      * recombination rate varies in 0.5, 0.001, 0.00001, 0.0000001 ?

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
