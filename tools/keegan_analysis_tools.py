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
    arguments to the variable_threeprime_map_function. This wrapper allows me to
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

# Create a function that finds the proteins I need. 
def find_transcripts(gene,transcripts, count_vec_m, count_vec_c):
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
                transcript = i
                index = transcripts.index(i)
                my_vector_m = count_vec_m[transcripts.index(i)]
                my_vector_c = count_vec_c[transcripts.index(i)]
                
    return transcript, my_vector_m, my_vector_c, index

# Create a function that finds the proteins I need. 
def find_trans_by_ID(ID,transcripts, count_vec_m, count_vec_c):
    '''
    A function that takes the name of a gene as input and finds 
    the corresponding transcript from a transcript list. 
    
    returns both the transcript in question and the vector of counts for that transcript.
    
    This function is still a work in progress as for now it simply gives the last 
    transcript in the list that matches the gene ID. 
    '''
    for i in transcripts:
        if i.attr['transcript_id'] == ID:
            transcript = i
            index = transcripts.index(i)
            my_vector_m = count_vec_m[transcripts.index(i)]
            my_vector_c = count_vec_c[transcripts.index(i)]
            break
                
    return transcript, my_vector_m, my_vector_c, index

# Create a function that finds the proteins I need. 
def find_trans_mmus(gene,transcripts, count_vec_m, count_vec_c):
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
                index = transcripts.index(i)
                my_vector_m = count_vec_m[transcripts.index(i)]
                my_vector_c = count_vec_c[transcripts.index(i)]
                
    return my_transcript, my_vector_m, my_vector_c, index

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

    for transcript in protein_coding:
        gene_id.append(transcript.attr["gene_name"])
        transcript_id.append(transcript.attr["transcript_id"])
        
    # Insert the gene ids and transcript ids into the codon_count list. 
    for i, j in zip(codon_counts, range(len(gene_id))):
        i.insert(0,gene_id[j])
        i.insert(0,transcript_id[j])
        
    # Calculate the longest cds region in our new list of counts
    l_tr = kat.find_max_list(codon_counts)

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

def get_smoothed_vector_parallel(vector, frac = 0.05):
    vector = vector + 0.0000000001
    positions = np.array(list(range(len(vector))))
    loess = Loess(positions, vector/sum(vector))
    smoothed_vec = []
    for x in positions:
        y = loess.estimate(x, window=int(len(positions)*frac), use_matrix=False, degree=1)
        smoothed_vec.append(y)
    smoothed_vec = np.array(smoothed_vec)
    cumsum = np.cumsum(smoothed_vec)
    return smoothed_vec, cumsum

#for tr in gtf_reads:
#    if tr.attr["transcript_biotype"] == "protein_coding":
#        if tr.attr["gene_name"] == "sqd":
#            print(gtf_reads.index(tr))
#            my_transcript = tr

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
    This function is used to create artificially elongation limited particle density profiles
    to be used as part of the TASEP-KS method. 
    
    The steps taken to create the artifically elongation limited particle densities are as follows:
   
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

def big_dif_mmus(diff_dist, transcripts, data_mutant, data_control, figsize = (16,50), fontsize = 12, stat_name = "ks stat ="):
    '''
    A function which creates a large graph showing the profile arrays for a list of transcripts
    
    returns a matplotlib axis object. 
    '''
    fig,ax = plt.subplots(len(diff_dist), 2, figsize = figsize)
    for axi, stat, gi in zip(ax, diff_dist, diff_dist.index):
            my_transcript, my_vec_mutant, my_vec_control, index = find_trans_mmus(gi, 
                                           transcripts, data_mutant, data_control)
            maxi = max([max(my_vec_mutant), max(my_vec_control)])*1.1

            axi[0].plot(my_vec_mutant)
            axi[0].text(len(my_vec_mutant)/2, maxi/1.2, stat_name + str(stat), fontsize = fontsize)
            axi[0].set_ylim([0,maxi+5])
            axi[0].set_ylabel("Read Counts", fontsize = fontsize)
            axi[0].set_xlabel("Codon Position", fontsize = fontsize)
            axi[0].set_title("mutant " + gi, fontsize = fontsize)
            axi[1].plot(my_vec_control)
            axi[1].set_ylim([0,maxi+5])
            axi[1].set_ylabel("Read Counts", fontsize = fontsize)
            axi[1].set_xlabel("Codon Position", fontsize = fontsize)
            axi[1].set_title("control " + gi, fontsize = fontsize)
    fig.tight_layout()
    
    return ax

def low_density(lamb,a,I = 10):
    '''
    A function that calculates the particle density along a transcript from a set of elongation rates
    inferred from ribosome profiling. This function assumes that elongation is in a low density regime
    (e.g. initiation limiting)
    '''
    Jl = (a*(lamb[0]-a))/(lamb[0] + (I-1)*a)
    pl = 1/(2*I) + (Jl*(I-1))/(2*I*lamb) - np.sqrt((1/(2*I) + (Jl*(I-1))/(2*I*lamb))**2 - Jl/(I*lamb))
    return(pl) 

def high_density(lamb,B,I = 10):
    '''
    A function that calculates the particle density along a transcript from a set of elongation rates
    inferred from ribosome profiling. This function assumes that elongation is in a high density regime
    (e.g. termination limiting)
    '''
    JR = (B*(lamb[-1]-B))/(lamb[-1] + (I-1)*B)
    pR = 1/(2*I) + (JR*(I-1))/(2*I*lamb) + np.sqrt((1/(2*I) + (JR*(I-1))/(2*I*lamb))**2 - JR/(I*lamb))
    return(pR) 

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
    return(p) 

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
    p = maximum_current(lamb_c,a=a,B=B,I = 10)
    return p

####!!!!#### If I am going to keep this function I better do something about the possibility of a high density regime.
def make_ld(lamb, a, B, I = 10):
    '''
    This function attempts to 
    '''
    lamb_c = copy.deepcopy(lamb)# Create a copy of lamb to work on. 
    Jmax = min(lamb_c)/((1+np.sqrt(I))**2)
    crit_a = ((lamb_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[0]*Jmax)/((lamb_c[0] - (I - 1)*Jmax)**2)))
    crit_B = ((lamb_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[-1]*Jmax)/((lamb_c[-1] - (I - 1)*Jmax)**2)))
    if a < crit_a and B > crit_B:
        p = low_density(lamb_c, a, I)
    else:
        while True:
            min_l = np.where(lamb_c == np.amin(lamb_c))[0][0]
            lamb_c[min_l] = lamb_c[min_l]*1.1 # It keeps doing this every run through. 
            Jmax = min(lamb_c)/((1+np.sqrt(I))**2)
            crit_a = ((lamb_c[0] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[0]*Jmax)/((lamb_c[0] - (I - 1)*Jmax)**2)))
            crit_B = ((lamb_c[-1] - (I-1) * Jmax) / 2)*(1 - np.sqrt(1 - (4*lamb_c[-1]*Jmax)/((lamb_c[-1] - (I - 1)*Jmax)**2)))
            if a < crit_a and B > crit_B:
                break
        p = low_density(lamb_c, a, I) 
    return p

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
        p = low_density(lamb_c, a, I)
        density = "LD"
    elif a > crit_a and B < crit_B:
        p = high_density(lamb_c, B, I)
        density = "HD"
    elif a < crit_a and B < crit_B:
        Jl = (a*(lamb[0]-a))/(lamb[0] + (I-1)*a)
        JR = (B*(lamb[-1]-B))/(lamb[-1] + (I-1)*B)
        sign = Jl - JR 
        if sign > 0:
            p = low_density(lamb_c, a, I)
            density = "LD"
        elif sign < 0:
            p = high_density(lamb_c, B, I)
            density = "HD"
    elif a > crit_a and B > crit_B:
        p = maximum_current(lamb_c, a, B, I)
        density = "MC"
    if intermediates == False:
        return p, density
    elif intermediates == True:
        return p, a, B, crit_a, crit_B, min(lamb_c), lamb_c[0]

# Create a function to automate the simulation process
def simulate_profile(mean_lambda, sd, length, a, B = 2, read_density = 1):
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
    
    '''
    mean_sqr = np.sqrt(mean_lambda)
    sample_size = int(read_density * length)
    lamb_c = np.random.gamma(mean_sqr, sd, length)+0.3
    
    # Create a set of random pause sights in your mutant. 
    pause_sites = np.random.randint(10, length/2, 4)
    lamb_m = copy.deepcopy(lamb_c)
    for ps in pause_sites:
        lamb_m[ps] = lamb_m[ps]* np.random.uniform(0.1, 0.8, 1)[0]
    
    # simulate the reads for the control
    p_c = get_density(lamb_c, a, B)
    prob_c = p_c[0]/sum(p_c[0])
    reads_c = np.zeros(length)
    for i in range(sample_size):
        x = numpy.random.choice(np.arange(0, len(prob_c)), p = prob_c)
        reads_c[x] = reads_c[x]+1

    # simulate the reads for the mutant assuming LD
    mu_a = a* np.random.uniform(0.25, 2, 1)[0]
    p_m = get_density(lamb_m, mu_a, B)
    prob_m = p_m[0]/sum(p_m[0])
    reads_m = np.zeros(length)
    for i in range(sample_size):
        x = numpy.random.choice(np.arange(0, len(prob_m)), p = prob_m)
        reads_m[x] = reads_m[x]+1
    return reads_c, reads_m
