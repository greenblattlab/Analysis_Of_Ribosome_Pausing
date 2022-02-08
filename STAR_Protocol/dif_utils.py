import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import csv
import math

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

def determine_enrichment(target, non_target, max_stat, N_secs, stat = "ks_stat"): 
    '''
    A function that determine the proportion of target genes from target that are found in
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