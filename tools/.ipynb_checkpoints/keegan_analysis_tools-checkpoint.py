""" Functions used regularly during my analysis to load and process data. 
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
import math
from tqdm import tqdm
import copy
from statsmodels.stats.proportion import proportions_ztest

def variable_threeprime_map_function(alignments,segment,p_offsets):
        '''
        This function is used to map read alignments to the location of the ribosomal p-site 
        from their 3' end. The offsets to use for each read length are specified by file
        generated using RiboWaltz.

        alignments:
            Information on the genome alignment of an individual read which is passed 
            to the function from a BamGenome array created by plastid. 

        segment:
            Information on the individual read segment which is passed 
            to the function from a BamGenome array created by plastid. 

        p_offsets:
            A pandas dataframe that has been loaded into the python environmemt.
            This dataframe should follow this template. 
                length          P_offsets
                 28              12
                 29              12
                 30              13
                ...             ...

        '''
        reads_out = []
        count_array = numpy.zeros(len(segment))
        for read in alignments: 
            for length, offset in zip(p_offsets["length"],p_offsets["p_offset"]): 
                if length != len(read.positions):
                    continue # skip read if it is not the length we are currently offsetting.

             # count offset 3' to 5' if the `segment` is on the plus-strand
             # or is unstranded
                if segment.strand == "+":
                    p_site = read.positions[-offset - 1]
                elif segment.strand == ".":
                    p_site = read.positions[-offset - 1]
             # count offset from other end if `segment` is on the minus-strand
                elif segment.strand == "-":
                    p_site = read.positions[offset]

                if p_site >= segment.start and p_site < segment.end:
                    reads_out.append(read)
                    count_array[p_site - segment.start] += 1
        return reads_out, count_array
    
def VariableThreePrimeMapFactory(p_offsets):
    '''
    BamGenome array objects will only be able to pass the alignments and segment
    arguments to the variable_threeprime_map_function. This wrapper all_gows me to
    also specify the offset that needs to be passed to the function. 
    '''
    def new_func(alignments,segment):
        return variable_threeprime_map_function(alignments,segment,p_offsets=p_offsets)

    return new_func

# Create a function that finds the proteins I need. 
def find_transcript(gene,transcripts, count_vectors):
    '''
    A function that takes the name of a gene as input and finds 
    the corresponding transcript from a transcript list. 
    
    returns both the transcript in question and the vector of counts for that transcript.
    
    This function is still a work in progress as for now it simply gives the last 
    transcript in the list that matches the gene ID. 
    '''
    
    for i in transcripts:
        if i.attr['transcript_biotype'] == 'protein_coding':
            if i.attr['gene_name'] == gene:
                my_transcript = i
                my_vector = count_vectors[transcripts.index(i)]
                index = transcripts.index(i)
                
    return my_transcript, my_vector, index

# Create a function that finds the proteins I need. 
def find_tran_mmus(gene,transcripts, count_vectors):
    '''
    A function that takes the name of a gene as input and finds 
    the corresponding transcript from a transcript list. 
    
    returns both the transcript in question and the vector of counts for that transcript.
    
    This function is still a work in progress as for now it simply gives the last 
    transcript in the list that matches the gene ID. 
    '''
    for i in transcripts:
        if i.attr['gene_name'] == gene:
            my_transcript = i
            my_vector = count_vectors[transcripts.index(i)]
            print(transcripts.index(i))
                
    return my_transcript, my_vector

def find_max_list(list):
    ''' 
    A function that finds the longest list/array in a list of lists. 
    '''
    list_len = [len(i) for i in list]
    return(max(list_len))

def tricubic(x):
    y = np.zeros_like(x)
    idx = (x >= -1) & (x <= 1)
    y[idx] = np.power(1.0 - np.power(np.abs(x[idx]), 3), 3)
    return y


class Loess(object):

    @staticmethod
    def normalize_array(array):
        min_val = np.min(array)
        max_val = np.max(array)
        if (max_val - min_val) == 0:
            return (array - min_val) / 1, min_val, max_val
        else:
            return (array - min_val) / (max_val - min_val), min_val, max_val

    def __init__(self, xx, yy, degree=1):
        self.n_xx, self.min_xx, self.max_xx = self.normalize_array(xx)
        self.n_yy, self.min_yy, self.max_yy = self.normalize_array(yy)
        self.degree = degree

    @staticmethod
    def get_min_range(distances, window):
        min_idx = np.argmin(distances)
        n = len(distances)
        if min_idx == 0:
            return np.arange(0, window)
        if min_idx == n-1:
            return np.arange(n - window, n)

        min_range = [min_idx]
        while len(min_range) < window:
            i0 = min_range[0]
            i1 = min_range[-1]
            if i0 == 0:
                min_range.append(i1 + 1)
            elif i1 == n-1:
                min_range.insert(0, i0 - 1)
            elif distances[i0-1] < distances[i1+1]:
                min_range.insert(0, i0 - 1)
            else:
                min_range.append(i1 + 1)
        return np.array(min_range)

    @staticmethod
    def get_weights(distances, min_range):
        max_distance = np.max(distances[min_range])
        if max_distance != 0:
            weights = tricubic(distances[min_range] / max_distance)
        else:
            weights = tricubic(distances[min_range] / 0.01)
        return weights

    def normalize_x(self, value):
        return (value - self.min_xx) / (self.max_xx - self.min_xx)

    def denormalize_y(self, value):
        return value * (self.max_yy - self.min_yy) + self.min_yy

    def estimate(self, x, window, use_matrix=False, degree=1):
        n_x = self.normalize_x(x)
        distances = np.abs(self.n_xx - n_x)
        min_range = self.get_min_range(distances, window)
        weights = self.get_weights(distances, min_range)

        if use_matrix or degree > 1:
            wm = np.multiply(np.eye(window), weights)
            xm = np.ones((window, degree + 1))

            xp = np.array([[math.pow(n_x, p)] for p in range(degree + 1)])
            for i in range(1, degree + 1):
                xm[:, i] = np.power(self.n_xx[min_range], i)

            ym = self.n_yy[min_range]
            xmt_wm = np.transpose(xm) @ wm
            beta = np.linalg.pinv(xmt_wm @ xm) @ xmt_wm @ ym
            y = (beta @ xp)[0]
        else:
            xx = self.n_xx[min_range]
            yy = self.n_yy[min_range]
            sum_weight = np.sum(weights)
            sum_weight_x = np.dot(xx, weights)
            sum_weight_y = np.dot(yy, weights)
            sum_weight_x2 = np.dot(np.multiply(xx, xx), weights)
            sum_weight_xy = np.dot(np.multiply(xx, yy), weights)
            if sum_weight != 0:
                mean_x = sum_weight_x / sum_weight
                mean_y = sum_weight_y / sum_weight
            else:
                mean_x = sum_weight_x / 0.01
                mean_y = sum_weight_y / 0.01
            denom = (sum_weight_x2 - mean_x * mean_x * sum_weight)
            if denom == 0:
                denom = 1e-20

            b = (sum_weight_xy - mean_x * mean_y * sum_weight) / \
                denom
            a = mean_y - b * mean_x
            y = a + b * n_x
        return self.denormalize_y(y)

def load_count_positions(csv_name, counts_path):
    # Create a list to hold the data and then fill it
    data = []
    gene_names = []
    with open(counts_path + csv_name, newline = '') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            data.append(row)
        
    # Remove the first row of the data (the header row)
    blank=data.pop(0)
            
    # Try to convert everything to a float if possible
    for i,ii in zip(data, range(len(data))):
        for j,jj in zip(i, range(len(i))):
            try:
                x = int(float(j))
                data[ii][jj] = x
            except:
                pass
            
    # Remove empty space
    for i,ii in zip(data, range(len(data))):
        x = list(filter(('').__ne__, i))
        data[ii] = x
        
    # Convert lists to np.arrays
    for i,ii in zip(data, range(len(data))):
        gene_names.append(data[ii][1])
        data[ii] = np.array(data[ii][2:])
    
    return data, gene_names

def save_count_positions(transcripts, codon_counts, save_path, save_name):
    """
    This function saves a list of count arrays as a csv file while adding on the gene ID and transcript ID to the file and adding
    a header that shows the position along the transcript for each count.
    """
    #create lists to hold the gene IDs and transcript IDs of the transcripts 
    gene_id = []
    transcript_id = []

    for transcript in transcripts:
        gene_id.append(transcript.attr["gene_name"])
        transcript_id.append(transcript.attr["transcript_id"])
        
    # Insert the gene ids and transcript ids into the codon_count list. 
    for i, j in zip(codon_counts, range(len(gene_id))):
        i.insert(0,gene_id[j])
        i.insert(0,transcript_id[j])
        
    # Calculate the longest cds region in our new list of counts
    l_tr = find_max_list(codon_counts)

    # Define a header that includes labels for the transcript and gene ID as 
    # well as numbers that index the cds region position.
    header=["transcript_id","gene_id"]+list(range(l_tr))

    # insert that header into our counts list. 
    codon_counts.insert(0,header)
    
    # Save the newly altered list as a csv. 
    with open(save_path + save_name, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(codon_counts)

def load_elongation_rates(csv_name, csv_path):
    data = []
    with open(csv_path + csv_name, newline = '') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            data.append(row)
    blank=data.pop(0)
            
    for i,ii in zip(data, range(len(data))):
        for j,jj in zip(i, range(len(i))):
            try:
                x = float(j)
                data[ii][jj] = x
            except:
                pass
            
    # Remove empty space
    for i,ii in zip(data, range(len(data))):
        x = list(filter(('').__ne__, i))
        data[ii] = x
        
    # Convert lists to np.arrays
    for i,ii in zip(data, range(len(data))):
        data[ii] = np.array(data[ii][2:])
    
    return data

# define a function that calculates the smoothed vector of the normalized reads
# using loess and calculates the cumulative sum of said vector.
def get_smoothed_vector(vector, frac = 0.05):
    vector = vector + 0.000000001
    positions = np.array(list(range(len(vector))))
    loess = Loess(positions, vector/sum(vector))
    smoothed_vec = []
    for x in positions:
        y = loess.estimate(x, window=int(len(positions)*frac), use_matrix=False, degree=1)
        smoothed_vec.append(y)
    smoothed_vec = np.array(smoothed_vec)
    cumsum = np.cumsum(smoothed_vec)
    return smoothed_vec, cumsum

# Create a function to obtain a normalized profile (p) of ribosome footprints.
def calculate_p(data):
    p_list=[]
    for i in data:
        i = i+1
        M = sum(i)
        p = i/M
        p_list.append(p)
    return(p_list)

# Calculate the smoothed density vector pbar for xth entry with length n-9
def calculate_pbar(p_list):
    pbar_list=[]
    for p in p_list:
        x=0
        pbar=[]
        for px in p:
            pbar_x = 0.1*sum(p[x:x+10]) #it is x+10 not x+9 because python does not include the final index.
            pbar.append(pbar_x)
            x = x+1
            if x  == len(p)-9:
                break
        pbar_list.append(np.array(pbar))
    return(pbar_list)

# calculate the smoothed, scaled elongation rate lambda bar 
def calculate_lbar(pbar_list):
    lbar_list=[]
    for pbar in pbar_list:
        lbar = []
        for pbarx in pbar:
            if pbarx == 0:
                lbar_x=9999
            else:
                lbar_x = (1-9*pbarx)/(pbarx*(1-pbarx))
            lbar.append(lbar_x)
        lbar_list.append(np.array(lbar))
    return(lbar_list)

def calculate_tau(lbar_list, codon_seq_list):
    tau_list = []
    for lbar, index in tqdm(zip(lbar_list, list(range(len(lbar_list))))):
        A = np.zeros((len(lbar),64))
        for row, i in zip(A, range(len(A))):
            set_of_10 = codon_seq_list[index][i:i+10] # what do I do with the index? 
            for j in set_of_10:
                row[j] = 1
        b = 10*lbar
        ls_result = lsqr(A,b)
        Ci = ls_result[0]
        tau = Ci.mean()
        tau_list.append(tau)
    
    return(tau_list)

def low_density(lamb,a,I):
    '''
    A function that calculates the particle density along a transcript from a set of elongation rates
    inferred from ribosome profiling. This function assumes that elongation is in a low density regime
    (e.g. initiation limiting)
    '''
    Jl = (a*(lamb[0]-a))/(lamb[0] + (I-1)*a)
    pl = 1/(2*I) + (Jl*(I-1))/(2*I*lamb) - np.sqrt((1/(2*I) + (Jl*(I-1))/(2*I*lamb))**2 - Jl/(I*lamb))
    return(pl) 

def high_density(lamb,B,I):
    '''
    A function that calculates the particle density along a transcript from a set of elongation rates
    inferred from ribosome profiling. This function assumes that elongation is in a high density regime
    (e.g. termination limiting)
    '''
    JR = (B*(lamb[-1]-B))/(lamb[-1] + (I-1)*B)
    pR = 1/(2*I) + (JR*(I-1))/(2*I*lamb) + np.sqrt((1/(2*I) + (JR*(I-1))/(2*I*lamb))**2 - JR/(I*lamb))
    return(pR) 

def maximum_current(lamb,a,B,I):
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
    return(p) 

def alter_p(arr_c, arr_m, I = 10):
    '''
    This function is used to create artificiall_gy elongation limited particle density profiles
    to be used as part of the TASEP-KS method. 
    
    The steps taken to create the artificall_gy elongation limited particle densities are as follows:
   
    1. determine the critical initiation and termination rates of the control transcrpt assuming maximum current
    2. arbitrarily set the initiation rate slightly below the critical initation rate and the termination rate 
       slightly above the critical termination rate in order to put the transcript in a low density phase.
    3. Find the maximum pause in the mutant transcript and then locate the pause furthest to the right of the transcript
       that is within 50% of the maximum pause site
    4. recursively lower the elongation rate in the control at the site of elongation limitation in the mutant until 
       elongation limitation occurs in the control
    5. Calculate the particle density of the new set of elongation rates for the control using the maximum current equations
    '''
    lam_c = copy.deepcopy(arr_c)
    lam_m = copy.deepcopy(arr_m)
    Jmax = min(lam_c)/((1+np.sqrt(I))**2)
    crit_a = ((lam_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lam_c[0]*Jmax)/((lam_c[0] - (I - 1)*Jmax)**2)))
    crit_B = ((lam_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lam_c[-1]*Jmax)/((lam_c[-1] - (I - 1)*Jmax)**2)))
    a = crit_a * 0.80
    B = crit_B * 1.2
    mut_min = np.amin(lam_m)
    el_p = np.amax(np.where(lam_m < mut_min*2)[0])
    while True:
        lam_c[el_p] = lam_c[el_p]*0.9 # It keeps doing this every run through. 
        Jmax = min(lam_c)/((1+np.sqrt(I))**2)
        crit_a = ((lam_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lam_c[0]*Jmax)/((lam_c[0] - (I - 1)*Jmax)**2)))
        crit_B = ((lam_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lam_c[-1]*Jmax)/((lam_c[-1] - (I - 1)*Jmax)**2)))
        if crit_a < a and crit_B < B:
            break
    p = maximum_current(lam_c,a=a,B=B,I = 10)
    return p

#Add some label parameters to this and then put it in kat for gods sake. 
def big_dif_smoothed(diff_dist, gene_names, data_mutant, data_control, figsize = (16,50), fontsize = 12, stat_name = "ks_stat ="):
    '''
    A function which creates a large graph showing the smoothed profile arrays for a list of transcripts
    
    returns a matplotlib axis object. 
    '''
    fig,ax = plt.subplots(len(diff_dist), 2, figsize = figsize)
    for axi, stat, gi in zip(ax, diff_dist, diff_dist.index):
        for tr_m, tr_c, name in zip(data_mutant, data_control, gene_names):
            if gi == name:
                my_vec_mutant = tr_m
                my_vec_control = tr_c
        sm_m, cumul_m = get_smoothed_vector(my_vec_mutant)
        sm_c, cumul_c = get_smoothed_vector(my_vec_control)
        maxi = max([max(sm_m), max(sm_c)])*1.1
        axi[0].plot(sm_m)
        axi[0].text(len(sm_m)/2, maxi/1.2, stat_name + str(stat), fontsize = fontsize)
        axi[0].set_ylim([0,maxi])
        axi[0].set_ylabel("Read Counts", fontsize = fontsize)
        axi[0].set_xlabel("Codon Position", fontsize = fontsize)
        axi[0].set_title("mutant " + gi, fontsize = fontsize)
        axi[1].plot(sm_c)
        axi[1].set_ylim([0,maxi])
        axi[1].set_ylabel("Read Counts", fontsize = fontsize)
        axi[1].set_xlabel("Codon Position", fontsize = fontsize)
        axi[1].set_title("control " + gi, fontsize = fontsize)
    fig.tight_layout()
            
    return ax

def big_dif(diff_dist, gene_names, data_mutant, data_control, figsize = (16,50), fontsize = 12, stat_name = "ks_stat ="):
    '''
    A function which creates a large graph showing the profile arrays for a list of transcripts
    
    returns a matplotlib axis object. 
    '''
    fig,ax = plt.subplots(len(diff_dist), 2, figsize = figsize)
    for axi, stat, gi in zip(ax, diff_dist, diff_dist.index):
        for tr_m, tr_c, name in zip(data_mutant, data_control, gene_names):
            if gi == name:
                my_vec_mutant = tr_m
                my_vec_control = tr_c
        maxi = max([max(my_vec_mutant), max(my_vec_control)])*1.1

        axi[0].plot(my_vec_mutant)
        axi[0].text(len(my_vec_mutant)/2, maxi/1.2, stat_name + str(stat), fontsize = fontsize)
        axi[0].set_ylim([0,maxi])
        axi[0].set_ylabel("Read Counts", fontsize = fontsize)
        axi[0].set_xlabel("Codon Position", fontsize = fontsize)
        axi[0].set_title("mutant " + gi, fontsize = fontsize)
        axi[1].plot(my_vec_control)
        axi[1].set_ylim([0,maxi])
        axi[1].set_ylabel("Read Counts", fontsize = fontsize)
        axi[1].set_xlabel("Codon Position", fontsize = fontsize)
        axi[1].set_title("control " + gi, fontsize = fontsize)
    fig.tight_layout()
            
    return ax

def split_equal(value, parts):
    '''
    A simple function that takes a number (value) and then divides that number into a certain number (parts) of equal parts
    '''
    value = float(value)
    return [i*value/parts for i in range(1,parts+1)]

def determine_enrichment(target, all_g, max_stat, N_secs, stat = "ks_stat"): 
    '''
    A function that determine the proportion of target genes from target that are found in
    all_g ks for a specified number of KS fractions. 
    '''
    ratios = []
    sections = split_equal(max_stat, N_secs)
    ratios.append(len(target[stat][target[stat] < sections[0]])/len(all_g[stat][all_g[stat] < sections[0]]))
    for sec, i in zip(sections, list(range(len(sections)))):
        try:
            ratios.append(len(target[stat][(target[stat] > sec) & (target[stat] < sections[i+1])]
                )/len(all_g[stat][(all_g[stat] > sec) & (all_g[stat] < sections[i+1])]))
        except:
            pass
    ratios.append(len(target[stat][target[stat] > sections[-1]])/len(all_g[stat][all_g[stat] > sections[-1]]))
    sections.insert(0,0)
    return ratios, sections

def det_p_values(target, all_g, sections, stat = "ks_stat"):
    '''
    A function that uses the proportion Z test to determine if the enrichment of the target gene is
    significant in any KS fractions. 
    '''
    p_values = []
    for sec, i in zip(sections, list(range(len(sections)))):
        try:
            obs = len(target[stat][(target[stat] > sec) & (target[stat] < sections[i + 1])])
            all_g_p = len(all_g[stat][(all_g[stat] > sec) & (all_g[stat] < sections[i + 1])])
            p_v = proportions_ztest(obs, all_g_p, len(target)/len(all_g))[1]
            p_values.append(p_v)
        except:
            pass
    obs = len(target[stat][target[stat] > sections[-1]])
    all_g_p = len(all_g[stat][all_g[stat] > sections[-1]])
    p_v = proportions_ztest(obs, all_g_p, len(target)/len(all_g))[1]
    p_values.append(p_v)
    return p_values