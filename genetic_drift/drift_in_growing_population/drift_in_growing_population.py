"""
This code is written to simulate the change in heterrozygosity in a population
that is following the logostic growth. First it calculate the population size
based on that it will find the effective population size. Lastlym it will use
the effective population size to track the changes in the population heterozygosity. 

"""


import argparse
import matplotlib.pyplot as plt

def logistic_growth(population_number,
                    growth_rate,
                    carrying_capacity):
    '''
    This function calculate the next population size based on it parameters.

    Parameters
    ----------
    population_number   :   Float
                            Current population size 
    growth_rate         :   Float
                            Growth rate (r)
    carrying_capacity   :   Float
                            The carrying capacity of population (k)

    Return
    ------
    Float 
                            the population size of the next time unit. 
    '''
    popsize = population_number + population_number * growth_rate * \
    (1 - (population_number / carrying_capacity))

    return popsize


def effective_population_size (population_number,
                               generations,
                               growth_rate,
                               carrying_capacity):
    '''
    This function calculate the effective population size of a growing population. 
    
    Parameters
    ----------
    population_number   :   Float
                            Current population size 
    growth_rate         :   Float
                            Growth rate (r)
    carrying_capacity   :   Float
                            The carrying capacity of population (k)
    generations         :   Int
                            The length of time during which the Ne is calculated. 

    Return
    ------
    Float 
                            the effective population size. 
    '''
    ne = 0
    for _ in range(generations):
        popsize = logistic_growth(population_number,
                                  growth_rate,
                                  carrying_capacity)
        ne += 1/popsize
    ne = generations / ne

    return ne

def heterozygosity (**kwargs):
    ne = effective_population_size (
        kwargs['n'],
        kwargs['gen'],
        kwargs['r'],
        kwargs['k']
    )
    heteroz = [kwargs['h']]
    for _ in range(kwargs['gen']):
        heteroz.append((1-(1/(2 * ne))) * heteroz[-1])
    return heteroz

def diagram_generator (heteros):
    plt.rcParams['figure.dpi'] = 250
    plt.rcParams['font.family'] = 'Times New Roman'
    fig, ax = plt.subplots(figsize = (8, 6))
    ax.set_title('\n Change of Heterozygosity \n ')
    ax.grid(True, linestyle = '-.')
    ax.plot(heteros)
    ax.set(ylim = (0, 0.5), xlabel = 'Generations', ylabel = 'Heterozygosity')
    plt.show()

if __name__ == '__main':
    parser = argparse.ArgumentParser()
    parser.add_argument("--GEN", type = int)
    parser.add_argument("--N", type = float)
    parser.add_argument("---R", type=float)
    parser.add_argument("--K", type=float)
    parser.add_argument("--H", type=float)
    args = parser.parse_args()

    K = args.K
    R = args.R
    N = args.N
    GEN = args.GEN
    H = args.H

    diagram_generator (heterozygocity(n = N, r = R, k = K, gen = GEN, h = H))
