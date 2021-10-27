""" Functions used to regularly in most analyses files 
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
                
    return my_transcript, my_vector

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
        weights = tricubic(distances[min_range] / max_distance)
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

            mean_x = sum_weight_x / sum_weight
            mean_y = sum_weight_y / sum_weight

            b = (sum_weight_xy - mean_x * mean_y * sum_weight) / \
                (sum_weight_x2 - mean_x * mean_x * sum_weight)
            a = mean_y - b * mean_x
            y = a + b * n_x
        return self.denormalize_y(y)
    
def load_count_positions(csv_name, csv_path):
    data = []
    with open(csv_path + csv_name, newline = '') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            data.append(row)
    blank=data.pop(0)
            
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
        data[ii] = np.array(data[ii][2:])
    
    return data

# define a function that calculates the smoothed vector of the normalized reads
# using loess and calculates the cumulative sum of said vector.
def get_smoothed_vector(positions, vector, frac = 0.05):
    loess = Loess(positions, vector/sum(vector))
    smoothed_vec = []
    for x in positions:
        y = loess.estimate(x, window=int(len(positions)*frac), use_matrix=False, degree=1)
        smoothed_vec.append(y)
    smoothed_vec = np.array(smoothed_vec)
    cumsum = np.cumsum(smoothed_vec)
    return smoothed_vec, cumsum


def big_dif(diff_dist, transcripts, data_mutant, data_control, figsize = (16,50), fontsize = 12):
    '''
    A function which creates a large graph showing the profile arrays for a list of transcripts
    
    returns a matplotlib axis object. 
    '''
    fig,ax = plt.subplots(len(diff_dist), 2, figsize = figsize)
    for axi, gi in zip(ax, diff_dist):
            my_transcript, my_vec_mutant, my_vec_control, index = find_transcripts(gi, 
                                           transcripts, data_mutant, data_control)
            maxi = max([max(my_vec_mutant), max(my_vec_control)])

            axi[0].plot(my_vec_mutant)
            axi[0].set_ylim([0,maxi+5])
            axi[0].set_title("mutant " + gi, fontsize = fontsize)
            axi[1].plot(my_vec_control)
            axi[1].set_ylim([0,maxi+5])
            axi[1].set_title("control " + gi, fontsize = fontsize)
            
    return ax

def big_dif_mmus(diff_dist, transcripts, data_mutant, data_control, figsize = (16,50), fontsize = 12):
    '''
    A function which creates a large graph showing the profile arrays for a list of transcripts
    
    returns a matplotlib axis object. 
    '''
    fig,ax = plt.subplots(len(diff_dist), 2, figsize = figsize)
    for axi, gi in zip(ax, diff_dist):
            my_transcript, my_vec_mutant, my_vec_control, index = find_trans_mmus(gi, 
                                           transcripts, data_mutant, data_control)
            maxi = max([max(my_vec_mutant), max(my_vec_control)])

            axi[0].plot(my_vec_mutant)
            axi[0].set_ylim([0,maxi+5])
            axi[0].set_title("mutant " + gi, fontsize = fontsize)
            axi[1].plot(my_vec_control)
            axi[1].set_ylim([0,maxi+5])
            axi[1].set_title("control " + gi, fontsize = fontsize)
            
    return ax

#for tr in gtf_reads:
#    if tr.attr["transcript_biotype"] == "protein_coding":
#        if tr.attr["gene_name"] == "sqd":
#            print(gtf_reads.index(tr))
#            my_transcript = tr