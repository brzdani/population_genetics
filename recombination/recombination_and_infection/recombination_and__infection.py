from itertools import combinations_with_replacement
from itertools import product
import pandas as pd
import numpy as np
from tqdm import tqdm

#------------------------------
#   Initializing the population
#------------------------------

def init_pop(initial_population_size):
    population = pd.DataFrame()
    genotype = []
    haplotype = [
        (i, j)
        for i, j in product(range(2), repeat=2)
    ]
    for h1, h2 in combinations_with_replacement(haplotype, 2):
        genotype.append((h1, h2))
        
    return population

def init_para (population):
    genotypes = []
    numbers = []
    for _, genotype in enumerate(population['Genotype']):
        numbers.append(np.random.rand())
        genotypes.append(genotype)
    numbers = np.array(numbers)
    return pd.DataFrame({
        ['Genotype'] : genotypes,
        ['Number']: numbers/sum(numbers)
    })

#----------------------
#   Population Dynamics
#----------------------

def infection_and_reproduction(pop, pop_size,  RATES):
    pass

def recombination():
    pass

#--------------------------
#   Calculating frequencies
#--------------------------

def genotype_freq():
    pass

#------------------
#   Main simulation
#------------------

def simulation(generations, initial_population_size, RATES):
    pop =  init_pop(initial_population_size)
    para = init_para(pop)
    popukation_size = [initial_population_size]

    for _ in tqdm(range(generations)):
         
        pop, para = infection_and_reproduction(pop, para, RATES)
        pop = recombination(pop)
        popukation_size.append(sum(pop))

#------------------------
#   Setting the paramters
#------------------------

INIT_POP_SIZE = None
BETA = None
ALPHA = None
B0 = None
SELEC_DIF = None
DENS_DEPEN = None
D = None
SEX_R = None
GENS = None 


RATES = {
    'beta': BETA,
    'Alpha': ALPHA,
    'Maximum Birth rate': B0, 
    'Selection differential': SELEC_DIF,
    'Density Dependance': DENS_DEPEN,
    'Death Rate': D,
    'Sex Ratio' : SEX_R
}


