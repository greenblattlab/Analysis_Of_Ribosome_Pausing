from plastid import BAMGenomeArray, Transcript
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import csv
import math

def variable_threeprime_map_function(alignments,segment,p_offsets):
        '''
        This function is used to map read alignments to the location of the ribosomal p-site 
        from their 3' end. The offsets to use for each read length are specified by a file
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
        count_array = np.zeros(len(segment))
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
    arguments to the variable_threeprime_map_function. This wrapper allows one to
    specify the offset that needs to be passed to the function. 
    '''
    def new_func(alignments,segment):
        return variable_threeprime_map_function(alignments,segment,p_offsets=p_offsets)

    return new_func

def find_max_list(list):
    ''' 
    A function that finds the longest list/array in a list of lists. 
    '''
    list_len = [len(i) for i in list]
    return(max(list_len))

def add_gene_ids(transcripts, count_vectors):
    #create lists to hold the gene IDs and transcript IDs of the transcripts 
    gene_id = []
    transcript_id = []

    for transcript in transcripts:
        gene_id.append(transcript.attr["gene_name"])
        transcript_id.append(transcript.attr["transcript_id"])
        
    # Insert the gene ids and transcript ids into the codon_count list. 
    for i, j in zip(count_vectors, range(len(gene_id))):
        i.insert(0,gene_id[j])
        i.insert(0,transcript_id[j])

def save_count_positions(count_vectors, path_to_save):
    # Calculate the longest cds region in our new list of counts
    l_tr = find_max_list(count_vectors)

    # Define a header that includes labels for the transcript and gene ID as 
    # well as numbers that index the cds region position.
    header=["transcript_id","gene_id"]+list(range(l_tr))

    # insert that header into our counts list. 
    count_vectors.insert(0,header)
    
    # Save the newly altered list as a csv. 
    with open(path_to_save, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(count_vectors)
        
def load_count_positions(path_to_csv):
    # Create a list to hold the data and then fill it
    data = []
    gene_names = []
    with open(path_to_csv, newline = '') as csvfile:
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

def get_smoothed_vector(vector, frac = 0.05):
    '''
    a function that calculates the smoothed vector of the normalized reads
    using loess and calculates the cumulative sum of said vector.
    '''
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

def split_equal(value, parts):
    '''
    A simple function that takes a number (value) and then divides that number into a certain number (parts) of equal parts
    '''
    value = float(value)
    return [i*value/parts for i in range(1,parts+1)]

def determine_enrichment(target, non_target, max_stat, N_secs, stat = "ks_stat"): 
    '''
    A function that determine the proportion of target genes from a target table that are found in
    non_target ks for a specified number of KS fractions. 
    '''
    frac_t = []
    frac_n = []
    enrich = []
    sections = split_equal(max_stat, N_secs)
    frac_t.append(len(target[stat][target[stat] < sections[0]])/len(target[stat]))
    frac_n.append(len(non_target[stat][non_target[stat] < sections[0]])/len(non_target[stat]))
    enrich.append(frac_t[0]/frac_n[0])
    for sec, i in zip(sections, list(range(len(sections)))):
        try:
            frac_t.append(len(target[stat][(target[stat] > sec) & (target[stat] < sections[i+1])]
                )/len(target[stat]))
            frac_n.append(len(non_target[stat][(non_target[stat] > sec) & (non_target[stat] < sections[i+1])]
                )/len(non_target[stat]))
            enrich.append(frac_t[i+1]/frac_n[i+1])
        except:
            pass
    frac_t.append(len(target[stat][target[stat] > sections[-1]])/len(target[stat]))
    frac_n.append(len(non_target[stat][non_target[stat] > sections[-1]])/len(non_target[stat]))
    enrich.append(frac_t[-1]/frac_n[-1])
    sections.insert(0,0)
    return enrich, sections

def Fisher_exact_p_values(target, non_target, sections, stat = "ks_stat"):
    '''
    A function that uses Fisher's exact test to determine if the enrichment of the target gene is
    significant in any KS fractions. 
    '''
    p_values = []
    for sec, i in zip(sections, list(range(len(sections)))):
        try:
            obs = len(target[stat][(target[stat] > sec) & (target[stat] < sections[i + 1])])
            non_target_p = len(non_target[stat][(non_target[stat] > sec) & (non_target[stat] < sections[i + 1])])
            table = [[obs,len(target)-obs],[non_target_p,len(non_target)-non_target_p]]
            p_v = stats.fisher_exact(table)[1]
            p_values.append(p_v)
        except:
            pass
    obs = len(target[stat][target[stat] > sections[-1]])
    non_target_p = len(non_target[stat][non_target[stat] > sections[-1]])
    table = [[obs,len(target)-obs],[non_target_p,len(non_target)-non_target_p]]
    p_v = stats.fisher_exact(table)[1]
    p_values.append(p_v)
    return p_values