"""
A python code for simulating a stabilizing slection for wing length.

In this program a dataframe is made to hold the information about every
genotype, icluding its phenotype and fitness. Then a population dataframe 
is made to keep track of frequcies of genotypes through geneations. In this, 
wing length in controlled by two loci on one chromose. so, we have recombination. 
and also the selecetion is stabilizing.

Author: Daniel Borzoo
Feb. 2026 
"""

import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['font.family'] = 'Times New Roman'

def genotype_df_generator (length_fitness_dic, init_length):
    '''
    Generate the genotype dataframe (describesd in module doc.)

    Parameters
    ----------
    lenghth_fitness_dic :   Dictionary
                            Contians related fitness (as values) 
                            to every phenotype (as key). 
    init_length :           Int. (can be float, too)
                            Set the initial value of the wing.

    Returns
    -------
    Data frame 
                            A genotype dataframe. 
    '''
    genotype = []
    index = []
    counter = 0
    for i in range(2):
        for j in range(2):
            for x in range(2):
                for y in range(2):
                    gen_temp = f'A{i+1} B{j+1} / A{x+1} B{y+1}'
                    if f'A{x+1} B{y+1} / A{i+1} B{j+1}' not in genotype:
                        genotype.append(gen_temp)
                        index.append(counter+1)
                        if f'A{x+1} B{y+1} / A{i+1} B{j+1}' != gen_temp:
                            genotype.append( f'A{x+1} B{y+1} / A{i+1} B{j+1}')
                            index.append(counter+1)
                        counter +=1
    fitness = []
    phenotype = []
    for i in genotype:
        count_2 = i.count('2')
        phenotype.append(init_length+count_2)
        fitness.append(length_fitness_dic.get(count_2, 1000))
    return pd.DataFrame({'Phenotype': phenotype, 'Fitness': fitness},
                         index=genotype)

def init_haplotype_frequncy(a1_freuqncy, b1_frequncy):
    '''
    Generate the haplotype dataframe based on initial frequncies of alleles.

    Parameters
    ----------
    a1_freuqncy :   float
                    The initial frequency of the A1 allele. 
    b1_frequncy :   floar (can be float, too)
                    The initial frequency of the B1 allele. 

    Returns
    -------
    Data frame 
                            A haplotype dataframe contianing 
                            the frequncy of every haplotype. 
    '''
    frequcncy_dic ={'A1': a1_freuqncy, 'A2': (1- a1_freuqncy),
                    'B1': b1_frequncy, 'B2' : (1-b1_frequncy)}
    h_frequency = []
    htype = []
    for i in range(1, 3):
        for j in range(1, 3):
            h_frequency.append(frequcncy_dic[f'A{i}'] * frequcncy_dic[f'B{j}'])
            htype.append(f'A{i} B{j}')
    haplotype_dataframe = pd.DataFrame({htype[i]:[h_frequency[i]] \
                                        for i in range(len(htype))},
                                        index=[1])
    return haplotype_dataframe

def population_frequency (population, haplotype_dataframe):
    '''
    Generate the poplulation dataframe containing the next generation 
    based on the last row of haplotype dataframe. This function add the next 
    generation genotypic frequencies as a new row at the end of the population
    data frame. 

    Parameters
    ----------
    population          :   Data Frame
                            The population data frame contaning the frequency
                            of genotypes from the begining of the simulation. 
    haplotype_dataframe :   Data Frame
                            The data frame containing every haplotype frequency
                            since the begening. 

    Returns
    -------
    Data frame 
                            A population dataframe contianing the new genotype
                            frequency as the last row. One row is added as the 
                            new frequencies.
    '''
    genotypes = [*population.columns]
    haplotypes = [i.split(' / ') for i in genotypes]
    frequency = []

    for i in range(len(genotypes)):

        frequency.append(
                float(
                haplotype_dataframe[haplotypes[i][0]].iloc[len(population)]
                ) * \
                float(
                haplotype_dataframe[haplotypes[i][1]].iloc[len(population)]
                )
            )

    population.loc[len(population)] = frequency
    return population

def selection (population, genotype_dataframe):
    '''
    Generate the poplulation dataframe containing the freuquencies after the 
    selection based of fitnesses in genotype dataframe. This function changes
    the last entry of population data frame.  

    Parameters
    ----------
    population          :   Data Frame
                            The population data frame contaning the frequency
                            of genotypes from the begining of the simulation. 
    genotype_dataframe  :   Data Frame
                            The data frame containing fitness of every genotype. 

    Returns
    -------
    Data frame 
                            A population dataframe contianing the new genotype
                            frequency after the selection at the last row. The 
                            last row is replaced with new frequencies. 
    '''
    genotypes = [*population.columns]

    mean_fitness = sum(
        float(population[i].iloc[len(population)-1]) * \
         float(genotype_dataframe['Fitness'][i]) for i in genotypes
        )

    frequency = [
        (float(population[i].iloc[len(population)-1]) * \
         float(genotype_dataframe['Fitness'][i])) / \
            mean_fitness \
                for i in genotypes
        ]

    population.loc[len(population)-1] = frequency
    return population

def haplotype_frequncy (haplotype_dataframe, population, r):
    '''
    Generate the haplotype dataframe containing the freuquencies of every
    haplotype after gamete production. This function changesthe last entry
    of data frame.  

    Parameters
    ----------
    haplotype_dataframe :   Data Frame
                            The data frame containing every haplotype frequency
                            since the begening. 
    population          :   Data Frame
                            The population data frame contaning the frequency
                            of genotypes from the begining of the simulation. 
    r                   :   float
                            recombination rate. 

    Returns
    -------
    Data frame 
                            A haplotype dataframe contianing the new haplotype
                            frequency after recombination and gamete production. 
                            A new row row is added for the new frequencies. 
    '''
    haplotypes = [*haplotype_dataframe.columns]
    next_hap_freq = []
    for i in haplotypes:
        temp = haplotypes.copy()
        temp.remove(i)
        next_hap_freq.append(float(population[f'{i} / {i}'].iloc[len(population)-1]))
        for j in temp:
            allele_1 = i.split()
            allele_2 = j.split()
            being = [i in allele_2 for i in allele_1]
            if being.count(True) > 0:
                next_hap_freq[-1]  = next_hap_freq[-1] + \
                    (float(population[f'{i} / {j}'].iloc[len(population)-1])*\
                      0.5)
                next_hap_freq[-1]  = next_hap_freq[-1] + \
                    (float(population[f'{j} / {i}'].iloc[len(population)-1])*\
                      0.5)
            else:
                next_hap_freq[-1]  = next_hap_freq[-1] + \
                    (float(population[f'{i} / {j}'].iloc[len(population)-1])*\
                      0.5 * (1-r))
                next_hap_freq[-1]  = next_hap_freq[-1] + \
                    (float(population[f'{j} / {i}'].iloc[len(population)-1])*\
                      0.5 * (1-r))
                next_hap_freq[-1]  = next_hap_freq[-1] + \
                (float(population[
                f'{allele_1[0]} {allele_2[1]} / {allele_2[0]} {allele_1[1]}'
                ].iloc[len(population)-1])* 0.5 * (r))
                next_hap_freq[-1]  = next_hap_freq[-1] + \
                (float(population[
                f'{allele_2[0]} {allele_1[1]} / {allele_1[0]} {allele_2[1]}'
                ].iloc[len(population)-1])* 0.5 * (r))

    haplotype_dataframe.loc[len(haplotype_dataframe)+1] = next_hap_freq
    return haplotype_dataframe

def histogram_diagram (population, init_length, generation):
    '''
    Generate and show the histogram showing the frequency of every phenotype in
    a specific generation. 

    Parameters
    ---------- 
    population  :   Data Frame
                    The population data frame contaning the frequency
                    of genotypes from the begining of the simulation. 
    init_length :   Int. (can be float, too)
                    Set the initial value of the wing.
    generation  :   The histogram show the frequency of every haplotype
                    in this generation. 

    Returns
    -------
    None.  
    '''
    phenotype = []
    freq = []
    for i in [*population.columns]:
        phenotype.append(init_length+i.count('2'))
        freq.append(float(population[i].iloc[generation]))

    phen_freq = {f'{i}': 0 for i in range(10, 15)}

    for i, ptype in enumerate(phenotype):
        match ptype:
            case 10:
                phen_freq['10'] += freq[i]
            case 11:
                phen_freq['11'] += freq[i]
            case 12:
                phen_freq['12'] += freq[i]
            case 13:
                phen_freq['13'] += freq[i]
            case 14:
                phen_freq['14'] += freq[i]

    return phen_freq


def phenotype_dataframe_generator(population, init_length):
    '''
    Create the phenotype dataframe that shows the frequency of every 
    phenotype from the begining. 

    Parameters
    ---------- 
    population  :   Data Frame
                    The population data frame contaning the frequency
                    of genotypes from the begining of the simulation. 
    init_length :   Int. (can be float, too)
                    Set the initial value of the wing.

    Returns
    -------
    Data Frame
                    The phenotype dataframe.   
    '''
    phenotypes = [10, 11, 12, 13, 14]
    phenotype_dataframe = pd.DataFrame({i:[] for i in phenotypes})
    j = len(population)-1
    for j in range(len(population)):
        phenotype = []
        freq = []
        for i in [*population.columns]:
            phenotype.append(init_length+i.count('2'))
            freq.append(float(population[i].iloc[j]))

        phen_freq = {f'{i}': 0 for i in range(10, 15)}

        for i, ptype in enumerate(phenotype):
            match ptype:
                case 10:
                    phen_freq['10'] += freq[i]
                case 11:
                    phen_freq['11'] += freq[i]
                case 12:
                    phen_freq['12'] += freq[i]
                case 13:
                    phen_freq['13'] += freq[i]
                case 14:
                    phen_freq['14'] += freq[i]
        phenotype_dataframe.loc[j] = [*phen_freq.values()]
    return phenotype_dataframe

def variance_calc (phenotype_dataframe):
    '''
    Create a list of phenotypic variance in every generation.  

    Parameters
    ---------- 
    phenotype_dataframe  :   Data Frame
                            The data frame contianing frequency of every phenotype
                            from the begining. 

    Returns
    -------
    List
                            A list of phenotypic variance in every generation
    '''
    var = []
    for i in range(len(phenotype_dataframe)):
        var.append(0)
        for j in range(len(phenotype_dataframe.iloc[i])):
            diff = [*phenotype_dataframe.columns][j] - 12
            var[i] = var[i] + ((diff**2) * [*phenotype_dataframe.iloc[i]][j])
    return var


def simulation_run (**kwargs):
    genotype_dataframe = genotype_df_generator(kwargs['length_fitness_dic'],
                                               kwargs['init_length'])
    haplotype_dataframe = init_haplotype_frequncy(kwargs['init_p_a1'],
                                                  kwargs['init_p_b1'])
    population = pd.DataFrame({ i:[] \
                               for i in [*genotype_dataframe.index]})
    population = population_frequency (population, haplotype_dataframe)
    for i in range(kwargs['time']):
        selected_pop = population.copy(deep = True)
        selected_pop = selection(selected_pop, genotype_dataframe)
        haplotype_dataframe = haplotype_frequncy(haplotype_dataframe,
                                                 selected_pop, kwargs['rec_rate'])
        population = population_frequency(population, haplotype_dataframe)
    return population, haplotype_dataframe

def ld_variance_diagram(parameters_dataframe):
    plt.rcParams['figure.dpi'] = 150
    fig, ax = plt.subplots(nrows=1, ncols= 2, figsize=(11, 6))
    fig.suptitle('''\n Stabilizing Selection: LD and Variance \n''', fontsize = 16)
    ax[0].grid(True, linestyle = '-.')
    ax[1].grid(True, linestyle = '-.')
    ax[0].plot(parameters_dataframe['LD'])
    ax[0].set_title('LD', fontsize = 12)
    ax[0].set(ylim = (-0.0005, 0.1), xlabel = 'Generation', ylabel = 'LD')
    ax[1].plot(parameters_dataframe['Phen_var'])
    ax[1].set_title('Variance', fontsize = 12)
    ax[1].set(ylim = (-0.005, max(parameters_dataframe['Phen_var'])+0.1),
            xlabel = 'Generation', ylabel = 'Phenotypic Variance')
    fig.tight_layout()
    plt.show()

def phenotype_histogram(population,init_length):
    plt.rcParams['figure.dpi'] = 150
    fig, ax = plt.subplots(nrows=3, ncols= 3, figsize=(11, 10), sharex=True, sharey=True)
    fig.suptitle('''\n Stabilizing Selection: Phenotypes Histogram \n''', fontsize = 16)
    counter = 0
    for i in range(0, 3):
        for j in range(0, 3):
            f_dic = histogram_diagram(population, init_length, (counter*10))
            ax[i, j].bar([*f_dic.keys()] , [*f_dic.values()])
            ax[i, j].set_title(f'Generation {counter*10}')
            if j == 0:
                ax[i, j].set( ylabel = 'Frequency')
            if i == 2:
                ax[i, j].set( xlabel = 'Phenotype (wing length)')
            counter +=1
    fig.tight_layout()
    plt.show()

def genotype_hitmap(population):
    mid = len(population) //2 
    df_1 = population.iloc[0:mid+1]
    df_2 = population.iloc[mid:]
    plt.rcParams['figure.dpi'] = 150
    fig , ax = plt.subplots(1,2, figsize = (16, 11))
    fig.suptitle('\n Stablizing Selection: Genotype Frequency \n', fontsize = 16)

    im1 = ax[0].imshow(
        df_1.values,
        aspect = 'auto',
        vmin = 0,
        vmax = 1
    )
    ax[0].set_xticks(
        range(len([*df_1.columns])),
        labels = [*df_1.columns],
        rotation = 45,
        ha="right",
        rotation_mode="anchor"
    )
    ax[0].set(ylabel = 'Generations')
    im2 = ax[1].imshow(
        df_2.values,
        aspect = 'auto',
        vmin = 0,
        vmax = 1
    )
    ax[1].set_xticks(
        range(len([*df_2.columns])),
        labels = [*df_2.columns],
        rotation = 45,
        ha="right",
        rotation_mode="anchor"
    )
    ax[1].set_yticks(
        range(0, df_2.shape[0]+1, 10),
        labels = list(range(df_1.shape[0], population.shape[0]+1, 10))
    )
    fig.colorbar(im1, ax = ax)
    plt.show()


TIME = 100
INIT_P_A1   = float(input('Please enter the initial frequency of A1 \n'))
INIT_P_B1   = float(input('Please enter the initial frequency of B1 \n'))
REC_RATE    = float(input('Please enter the recombination rate \n'))
fitness_dic = {0: 0.7, 1:0.9, 2:1, 3:0.9, 4:0.7}
INIT_LENGTH = 10

pop, haplotype = simulation_run (time=TIME, init_p_a1=INIT_P_A1,
                                                  init_p_b1=INIT_P_B1, rec_rate=REC_RATE,
                                                  length_fitness_dic=fitness_dic,
                                                  init_length=INIT_LENGTH)

parameters_dataframe = pd.DataFrame(index =
                                    range(1, len(haplotype)+1))
parameters_dataframe['LD'] =  ( haplotype['A2 B1'] * \
                               haplotype['A1 B2']) - \
                                (haplotype['A2 B2'] * \
                                 haplotype['A1 B1'])
parameters_dataframe['Phen_var'] = variance_calc (
    phenotype_dataframe_generator(pop, INIT_LENGTH)
    )

ld_variance_diagram(parameters_dataframe)
phenotype_histogram(pop,INIT_LENGTH)
genotype_hitmap(pop)
