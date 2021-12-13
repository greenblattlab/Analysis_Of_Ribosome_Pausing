""" Functions used regularly during my analysis to integrate the equations of the inhomogenous l-TASEP model into my work. 
"""

from plastid import BAMGenomeArray, VariableFivePrimeMapFactory, \
                        GTF2_TranscriptAssembler, GFF3_TranscriptAssembler, \
                        Transcript, ThreePrimeMapFactory
import numpy as np
from Bio import SeqIO
import numpy
import pandas as pd
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
import matplotlib.pyplot as plt
import csv
from scipy.sparse.linalg import lsqr
import time
import random
import math
from tqdm import tqdm
import copy
from sympy import symbols, solve, sqrt

def low_density(lamb,a,I = 10):
    '''
    A function that calculates the particle density along a transcript from a set of elongation rates
    inferred from ribosome profiling. This function assumes that elongation is in a low density regime
    (e.g. initiation limiting)
    '''
    Jl = (a*(lamb[0]-a))/(lamb[0] + (I-1)*a)
    pl = 1/(2*I) + (Jl*(I-1))/(2*I*lamb) - np.sqrt((1/(2*I) + (Jl*(I-1))/(2*I*lamb))**2 - Jl/(I*lamb))
    return pl, Jl 

def high_density(lamb,B,I = 10):
    '''
    A function that calculates the particle density along a transcript from a set of elongation rates
    inferred from ribosome profiling. This function assumes that elongation is in a high density regime
    (e.g. termination limiting)
    '''
    JR = (B*(lamb[-1]-B))/(lamb[-1] + (I-1)*B)
    pR = 1/(2*I) + (JR*(I-1))/(2*I*lamb) + np.sqrt((1/(2*I) + (JR*(I-1))/(2*I*lamb))**2 - JR/(I*lamb))
    return pR, JR

def maximum_current(lamb,a,B,I = 10):
    '''
    A function that calculates the particle density along a transcript from a set of elongation rates
    inferred from ribosome profiling. This function assumes that elongation is in a maximum current regime
    (e.g. elongation limiting)
    '''
    Jmax = min(lamb)/((1+np.sqrt(I))**2)
    flip = np.where(lamb == np.amin(lamb))[0][0]
    pR = 1/(2*I) + (Jmax*(I-1))/(2*I*lamb[0:flip]) + np.sqrt((1/(2*I) + (Jmax*(I-1))/(
        2*I*lamb[0:flip]))**2 - Jmax/(I*lamb[0:flip]))
    pl = 1/(2*I) + (Jmax*(I-1))/(2*I*lamb[flip:]) - np.sqrt((1/(2*I) + (Jmax*(I-1))/(
        2*I*lamb[flip:]))**2 - Jmax/(I*lamb[flip:]))
    p = np.concatenate((pR,pl))
    return p, Jmax

def make_mc(arr_c, position, a, B, I = 10):
    '''
    This function purposefully induces elongation limitation at a certain point 
    '''
    lamb_c = copy.deepcopy(arr_c)
    Jmax = min(lamb_c)/((1+np.sqrt(I))**2)
    crit_a = ((lamb_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[0]*Jmax)/((lamb_c[0] - (I - 1)*Jmax)**2)))
    crit_B = ((lamb_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[-1]*Jmax)/((lamb_c[-1] - (I - 1)*Jmax)**2)))
    mut_min = position
    while True:
        lamb_c[mut_min] = lamb_c[mut_min]*0.9 # It keeps doing this every run through. 
        Jmax = min(lamb_c)/((1+np.sqrt(I))**2)
        crit_a = ((lamb_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[0]*Jmax)/((lamb_c[0] - (I - 1)*Jmax)**2)))
        crit_B = ((lamb_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[-1]*Jmax)/((lamb_c[-1] - (I - 1)*Jmax)**2)))
        if crit_a < a and crit_B < B:
            break
    p, J = maximum_current(lamb_c,a=a,B=B,I = 10)
    return p, J

####!!!!#### This function is unsure of what to do in the presence of a high density regime. 
def make_ld(lamb, a, B, I = 10):
    '''
    This function attempts to force everything to be in a low density regime regardless 
    of the initiation rate and elongation rate inputs. It does this by increasing the minimum elongation rate. 
    '''
    lamb_c = copy.deepcopy(lamb)# Create a copy of lamb to work on. 
    Jmax = min(lamb_c)/((1+np.sqrt(I))**2)
    crit_a = ((lamb_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[0]*Jmax)/((lamb_c[0] - (I - 1)*Jmax)**2)))
    crit_B = ((lamb_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[-1]*Jmax)/((lamb_c[-1] - (I - 1)*Jmax)**2)))
    if a < crit_a and B > crit_B:
        p, J = low_density(lamb_c, a, I)
    else:
        while True:
            min_l = np.where(lamb_c == np.amin(lamb_c))[0][0]
            lamb_c[min_l] = lamb_c[min_l]*1.1 # It keeps doing this every run through. 
            Jmax = min(lamb_c)/((1+np.sqrt(I))**2)
            crit_a = ((lamb_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[0]*Jmax)/((lamb_c[0] - (I - 1)*Jmax)**2)))
            crit_B = ((lamb_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[-1]*Jmax)/((lamb_c[-1] - (I - 1)*Jmax)**2)))
            if a < crit_a and B > crit_B:
                break
        p, J = low_density(lamb_c, a, I) 
    return p, J

def get_density(lamb, a, B, I = 10, intermediates = False):
    '''A function that determines the correct elongation  regime/phase and calculates the correct density 
       accordingly
       
       arg lamb: A set of elongation rates for each codon in a transcript
       type: list of floats
       
       arg a: An initiation rate (the rate at which ribosomes are added to the transcript)
       type: float
       
       arg B: A termination rate (the rate a which ribosomes leave the transcript)
       type: float
       
       arg I: The particle size to be used in the TASEP model (should equal the size of the ribosome in codons)
       type: int
       
       arg intermediates: If set to true than the function will output all intermediate values used to calculate densities
       type: bool 
    
    '''
    lamb_c = copy.deepcopy(lamb)# Create a copy of lamb to work on. 
    Jmax = min(lamb_c)/((1+np.sqrt(I))**2)
    crit_a = ((lamb_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[0]*Jmax)/((lamb_c[0] - (I - 1)*Jmax)**2)))
    crit_B = ((lamb_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[-1]*Jmax)/((lamb_c[-1] - (I - 1)*Jmax)**2)))
    if a < crit_a and B > crit_B:
        p, J = low_density(lamb_c, a, I)
        phase = "LD"
    elif a > crit_a and B < crit_B:
        p, J = high_density(lamb_c, B, I)
        phase = "HD"
    elif a < crit_a and B < crit_B:
        Jl = (a*(lamb[0]-a))/(lamb[0] + (I-1)*a)
        JR = (B*(lamb[-1]-B))/(lamb[-1] + (I-1)*B)
        sign = Jl - JR 
        if sign > 0:
            p, J = low_density(lamb_c, a, I)
            phase = "LD"
        elif sign < 0:
            p, J = high_density(lamb_c, B, I)
            phase = "HD"
    elif a > crit_a and B > crit_B:
        p, J = maximum_current(lamb_c, a, B, I)
        phase = "MC"
    if intermediates == False:
        return p, J, phase
    elif intermediates == True:
        return p, J, phase, a, B, crit_a, crit_B, min(lamb_c), lamb_c[0], lamb_c[-1]
    
def get_all_intermediates(mean_lambda = 4, sd = 3, min_lambda = 0.8, length = 400, a = 0.02, B = 2, read_density = 1, pause_N = 4, 
                    pause_str = 0.5, a_frac = 1.0, rng_a = False, rng_pause = False, rng_p_range = (0.1, 0.8),
                    rng_a_range = (0.25,2), elon_frac = 1):
    '''
    A function that simulates ribosome profiling data from a mutant and a control for a single gene. 
    
    mean_lambda: The mean elongation rate used for the simulated gene
    type: float
    
    sd: The standard deviation of the elongation rates for your simulated gene
    type: float
    
    length: The length of your simulated gene
    type: int
    
    a: The initiation rate of your simulated gene
    type: float
    
    B: The termination rate of your simulated gene
    type: float
    
    read_density: how many reads should be extracted from your simulated gene per codon. 
    type: int
    
    pause_N: An integer specifying the number of pauses you want to induce in the mutant. Note that all 
             pause sites are induced in the first half of the transcript.
    type: int
    
    pause_str: The fraction by which to change the elongation rates at the pause sites in the mutant.
    type: float
    
    a_frac: The fraction by which to change the initiation rates in the mutant
    
    rng_a: If this parameter is set to True then the fraction by which to change the initiation rates in the mutant
           will be chosen randomly from a range of possible fractions specified in "rng_a_range" rather than being a single value. 
    type: boolean
    
    rng_pause: If this parameter is set to True then the fraction by which to change the elongation rates at the pause sites in the mutant.
               will be chosen randomly from a range of possible fractions specified in "rng_p_range" rather than being a single value. 
    type: boolean
    
    rng_p_range: A range of values from which the pause strength of the pause sites will be selected from
    type: tuple of length 2
    
    rng_a_range: A range of values from which the change to the alpha value for the mutant will be selected from
    type: tuple of length 2
    
    elon_frac: The fraction by which to change the global elongation rates in the mutant. 
    type: float
    
    '''
    mean_sqr = np.sqrt(mean_lambda)
    sample_size = int(read_density * length)
    lamb_c = np.random.gamma(mean_sqr, sd, length) + min_lambda 
    lamb_m = copy.deepcopy(lamb_c)
        
    # Create a set of random pause sights in your mutant. 
    pause_sites = np.random.randint(10, length/2, pause_N)
    if rng_pause == True:
        for ps in pause_sites:
            lamb_m[ps] = lamb_m[ps]* np.random.uniform(rng_p_range[0], rng_p_range[1], 1)[0]
    else:
        for ps in pause_sites:
            lamb_m[ps] = lamb_m[ps]*pause_str

    # alter the initiation rate of the mutant 
    if rng_a == True:
        mu_a = a* np.random.uniform(rng_a_range[0], rng_a_range[1], 1)[0]
    else:
        mu_a = a*a_frac
        
    # alter the elongation rate of the mutant
    lamb_m = lamb_m*elon_frac
    
    # simulate the reads for the control
    p_c, J_c, phase_c, a_c, B_c, crit_a_c, crit_B_c, min_l_c, ini_l_c, term_l_c = get_density(lamb_c, a, B, intermediates = True)
    prob_c = p_c/sum(p_c)
    reads_c = np.zeros(length)
    for i in range(sample_size):
        x = numpy.random.choice(np.arange(0, len(prob_c)), p = prob_c)
        reads_c[x] = reads_c[x]+1

    # simulate the reads for the mutant assuming LD
    mu_a = a* np.random.uniform(0.25, 2, 1)[0]
    p_m, J_m, phase_m, a_m, B_m, crit_a_m, crit_B_m, min_l_m, ini_l_m, term_l_m = get_density(lamb_m, mu_a, B, intermediates = True)
    prob_m = p_m/sum(p_m)
    reads_m = np.zeros(length)
    for i in range(sample_size):
        x = np.random.choice(np.arange(0, len(prob_m)), p = prob_m)
        reads_m[x] = reads_m[x]+1
    return(reads_c, J_c, phase_c, a_c, B_c, crit_a_c, crit_B_c, min_l_c, ini_l_c, term_l_c, reads_m, J_m, phase_m, a_m, B_m, 
           crit_a_m, crit_B_m, min_l_m, ini_l_m, term_l_m)

# Create a function to automate the simulation process
def simulate_profile(mean_lambda = 4, sd = 3, min_lambda = 0.8, length = 400, a = 0.02, B = 2, read_density = 1, pause_density = 0.01, 
                     pause_str = 0.5, a_frac = 1.0, rng_a = False, rng_pause = False, rng_p_range = (0.1, 0.8),
                     rng_a_range = (0.25,2), elon_frac = 1, return_min_lam = False):
    '''
    A function that simulates ribosome profiling data from a mutant and a control for a single gene. 
    
    mean_lambda: The mean elongation rate used for the simulated gene
    type: float
    
    sd: The standard deviation of the elongation rates for your simulated gene
    type: float
    
    length: The length of your simulated gene
    type: int
    
    a: The initiation rate of your simulated gene
    type: float
    
    B: The termination rate of your simulated gene
    type: float
    
    read_density: how many reads should be extracted from your simulated gene per codon. 
    type: int
    
    pause_N: An integer specifying the number of pauses you want to induce in the mutant. Note that all 
             pause sites are induced in the first half of the transcript.
    type: int
    
    pause_str: The fraction by which to change the elongation rates at the pause sites in the mutant.
    type: float
    
    a_frac: The fraction by which to change the initiation rates in the mutant
    
    rng_a: If this parameter is set to True then the fraction by which to change the initiation rates in the mutant
           will be chosen randomly from a range of possible fractions specified in "rng_a_range" rather than being a single value. 
    type: boolean
    
    rng_pause: If this parameter is set to True then the fraction by which to change the elongation rates at the pause sites in the mutant.
               will be chosen randomly from a range of possible fractions specified in "rng_p_range" rather than being a single value. 
    type: boolean
    
    rng_p_range: A range of values from which the pause strength of the pause sites will be selected from
    type: tuple of length 2
    
    rng_a_range: A range of values from which the change to the alpha value for the mutant will be selected from
    type: tuple of length 2
    
    elon_frac: The fraction by which to change the global elongation rates in the mutant. 
    type: float
    
    '''
    
    mean_sqr = np.sqrt(mean_lambda)
    sample_size = int(read_density * length)
    lamb_c = np.random.gamma(mean_sqr, sd, length)+ min_lambda 
    lamb_m = copy.deepcopy(lamb_c)
        
    # Create a set of random pause sights in your mutant. 
    pause_sites = []
    for i in range(round(length)):
        if random.random() < pause_density:
            pause_sites.append(i)
    if rng_pause == True:
        for ps in pause_sites:
            lamb_m[ps] = lamb_m[ps]* np.random.uniform(rng_p_range[0], rng_p_range[1], 1)[0]
    else:
        for ps in pause_sites:
            lamb_m[ps] = lamb_m[ps]*pause_str

    # alter the initiation rate of the mutant 
    if rng_a == True:
        mu_a = a* np.random.uniform(rng_a_range[0], rng_a_range[1], 1)[0]
    else:
        mu_a = a*a_frac
        
    # alter the elongation rate of the mutant
    lamb_m = lamb_m*elon_frac
        
    # Simulate the reads for the control
    p_c, J_c, phase_c = get_density(lamb_c, a, B)
    prob_c = p_c/sum(p_c)
    reads_c = np.zeros(length)
    for i in range(sample_size):
        x = numpy.random.choice(np.arange(0, len(prob_c)), p = prob_c)
        reads_c[x] = reads_c[x]+1
        
    # simulate the reads for the mutant
    p_m, J_m, phase_m = get_density(lamb_m, mu_a, B)
    prob_m = p_m/sum(p_m)
    reads_m = np.zeros(length)
    for i in range(sample_size):
        x = np.random.choice(np.arange(0, len(prob_m)), p = prob_m)
        reads_m[x] = reads_m[x]+1
    if return_min_lam == False:
        return reads_c, J_c, phase_c, reads_m, J_m, phase_m
    elif return_min_lam == True:
        return reads_c, J_c, phase_c, min(lamb_c), reads_m, J_m, phase_m, min(lamb_m)
    
def determine_sim_enrichment(ks_table, N_cats, max_ks):
    ratios = []
    all_ks = ks_table
    ks_MC = ks_table[ks_table["phase_mutant"] == "MC"]
    sections = split_equal(max_ks, N_cats)
    ratios.append(len(ks_MC.ks_stat[ks_MC.ks_stat < sections[0]])/len(all_ks.ks_stat[all_ks.ks_stat < sections[0]]))
    for sec, i in zip(sections, list(range(len(sections)))):
        try:
            ratios.append(len(ks_MC.ks_stat[(ks_MC.ks_stat > sec) & (ks_MC.ks_stat < sections[i+1])]
                )/len(all_ks.ks_stat[(all_ks.ks_stat > sec) & (all_ks.ks_stat < sections[i+1])]))
        except:
            pass
    ratios.append(len(ks_MC.ks_stat[ks_MC.ks_stat > sections[-1]])/len(all_ks.ks_stat[all_ks.ks_stat > sections[-1]]))
    return ratios, sections

def split_equal(value, parts):
    value = float(value)
    return [i*value/parts for i in range(1,parts+1)]

def big_dif_sim(diff_dist, data_mutant, data_control, figsize = (16,50), fontsize = 12, stat_name = "ks_stat ="):
    '''
    A function which creates a large graph showing the profile arrays for a list of transcripts
    
    returns a matplotlib axis object. 
    '''
    fig,ax = plt.subplots(len(diff_dist), 2, figsize = figsize)
    for axi, stat, gi in zip(ax, diff_dist, diff_dist.index):
            my_vec_mutant = data_mutant[gi]
            my_vec_control = data_control[gi]
            maxi = max([max(my_vec_mutant), max(my_vec_control)])*1.1

            axi[0].plot(my_vec_mutant)
            axi[0].text(len(my_vec_mutant)/2, maxi/1.2, stat_name + str(stat), fontsize = fontsize)
            axi[0].set_ylim([0,maxi])
            axi[0].set_ylabel("Read Counts", fontsize = fontsize)
            axi[0].set_xlabel("Codon Position", fontsize = fontsize)
            axi[0].set_title("mutant " + str(gi), fontsize = fontsize)
            axi[1].plot(my_vec_control)
            axi[1].set_ylim([0,maxi])
            axi[1].set_ylabel("Read Counts", fontsize = fontsize)
            axi[1].set_xlabel("Codon Position", fontsize = fontsize)
            axi[1].set_title("control " + str(gi), fontsize = fontsize)
    fig.tight_layout()
            
    return ax

# create a function that can determine the minimum elongation rate necessary to make alpha/crit_alpha equal to 1 (the minimum elongation rate necesary for a phase change)
def get_crit_lambda(alpha, l1, I = 10):
    '''
    This function calculates the minimum elongation rate that would be necessary to make the quotient of alpha divided by the critical alpha equal to one. 
    This means that this function outputs the minimum elongation rate that is necessary for a phase change to occur. 
    '''
    lmin = symbols('lmin', positive = True, real = True)
    expr = ((l1 - (I-1) * (lmin/((1+sqrt(I))**2))) / 2)*(1 - sqrt(1 - (4*l1*(lmin/((1+sqrt(I))**2)))/((l1 - (I - 1)*(lmin/((1+sqrt(I))**2)))**2))) - alpha
    sol = solve(expr)
    return sol