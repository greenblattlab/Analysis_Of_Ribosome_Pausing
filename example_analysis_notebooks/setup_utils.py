from plastid import BAMGenomeArray, Transcript
import numpy as np
import csv

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
    '''
    
    '''
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
    '''
    
    '''
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