# Stabilizing Selection Model
This folder contains a simulation of model of stabilizing selection acting on a quantitative trait, with explicit incorporation of genetic recombination. This model is conceptually inspired by classic study of natural Selection by **Hermon C. Bumpus (1899)**.

## introduction
This folder contains a simulation model of stabilizing selection acting 
In one of the earliest attempts to measure the natural selection, Bumpus (1899) measured the wing length of sparrows killed in a storm, and those that survived. He found that the survivors contianed and excess of birds with wings of average length, and a deficiency of birds with very long or very short wings. Many subsequent measurements of natural selection on quantitative traits have given the same picture: 
**Typical individuals do better than either extreme.** Such selection is referred to as **Normalizing** or **Stabilizing** selection.
In natural populations, however, phenotypic distributions are shaped not only by selection but also by genetic architecture, includingL 
- Multiple loci contributing to a trait
- Recombination breaking and reshuffling allelic combinations
This simulation aims to capture these processes in a simplified but biologically interpretable framework. 

## Model behind the simulation
### Trait architecture 
The focal trait is quantitative, determined by multiple genetic loci. Each locus contributes additively to the phenotype. ab /ab individuals has a wing length of 10 units, and every substitution of A for a or B for b, adds one unit, so that AB / AB is 14 units. Alleles have fixed effects, and environmental is ignored. 

### Selection 
Fitness is modeled as a Gaussian functiton centered on an optimal phenotype (12 units of length). Individuals with trait values closer to the optimum have higher survival or reproductive success. In this model selection acts on adults. The strength of stabilizing selection is controlled by the variance of the fitness function. 
##### Selection Procedure
Selection is implemented based on individual phenotypes. Fitness in calculated for each individual, fitness values are normalized relative to the population. Individuals contribute to the mating pool proportional to their fitness. 
For each individual, the phenotype $z$ is calculated as the sum of genetic contributions across loci, and its frequency after selection is calculated as follows: 

$$f(z)' = \frac{w_z \cdot f(z) } {\bar{w}}$$

where 
* $w_z$ is the average fitness of individuals with phenotype $z$.
* $f(z)$ is the frequency of individuals with phenotype $z$.
* $ \bar{w}$ is the average fitness of the population.
In this simulation individuals have different viability, but all the surviving individuals have the same fecunditiy. 
### Recombination 
In this model, sexual reproduction is assumed. Recombination occures between loci at a specified recombination rate. In this simulation we track phenotypic variance and LD. For a detailed review on the mathematics behind the recombination, visit any evolution text beook in reach.
