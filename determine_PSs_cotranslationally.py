"""
Co-transcriptionally folds the genome. Tested on HPeV genome sets provided.

Dependencies:   Bio     (python-biopython)      availble via aptitude
                forgi   (python-forgi)          available at https://pypi.python.org/pypi/forgi/0.20

Author:         German Leonov
                University of York
                german.leonov@york.ac.uk
                
Created:        13/03/2017
Last modified:  04/07/2017
"""

import multiprocessing
from operator import itemgetter
import os
import re
import subprocess
import sys
import time, datetime

from Bio import SeqIO
from Bio.Seq import Seq

import forgi.graph.bulge_graph as cgb


class Strain:
    def __init__(self, name, sequence):
        self.name = name                            # strain identifier, typically the GenBank accession number
        self.sequence = sequence                    # the genomic sequence of the RNA
        self.SLs = []                               # contains the onjects of the SL class, i.e. stem loops


class SL:
    def __init__(self, structure, start, end):
        self.structure = structure                  # symmetric structure
        self.start = start                          # start position of the first basepair
        self.end = end                              # end position of the last basepair
        self.motif_position = None                  # start position of the recognition motif
        self.apical_loop_position = None            # start position of the apical loop
        self.basepairs = set()                      # set of basepair positions
        
        self.min_free_energy = None                 # Gibbs free energy of structure formation
        self.fold_count = 1                         # individual stability based on sampling (DO NOT ADJUST)
        self.fold_stability = 1                     # shared stability due to basepairing similarity to other structures (using other.fold_count, adjusted for recognition motif matching)
        self.structures = []                        # shared structures
        self.SLs_in_common = []                     # shared SLs (pointers)
        self.structures_common_basepairs = []       # shared structures common basepairs
    
    def check_if_same_SL(self, other):
        """
        Determines if one SL is the same as the other SL.
        """
        if self.start == other.start and self.structure == other.structure:
            return True
        else:
            return False
    
    def get_apical_loop_info(self, genome):
        """
        Returns the apical loop sequence.
        """
        match = re.finditer('\(\.+?\)', self.structure).next()
        self.apical_loop_position = self.start+match.start()+1
        return genome[self.start+match.start()+1:self.start+match.end()-1], self.apical_loop_position


class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


def time_stamp():
    """
    Returns the date and time.
    """
    return '%s%s%s' % (color.RED, datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), color.END)


def run_jobs(jobs_to_do, max_CPUs):
    """
    Runs jobs that need to be submitted for processing.
    """
    running_jobs = {}
    total_jobs = len(jobs_to_do)

    print time_stamp(), '%sTotal jobs to do: %s%d%s' % (color.BOLD, color.BLUE,  total_jobs, color.END)
    
    while len(jobs_to_do) > 0:
        finished_jobs = []

        # identify finished jobs
        for job in running_jobs:
            alive = job.is_alive()
            if not alive:
                finished_jobs.append(job)
                

        # remove finished jobs from running jobs list
        for job in finished_jobs:
            running_jobs.pop(job)

        # try start more jobs
        jobs_to_run = max_CPUs - len(running_jobs)
        if jobs_to_run > 0:
            for job in range(jobs_to_run):
                if len(jobs_to_do) > 0:
                    new_job = jobs_to_do.pop(0)
                    new_job.start()
                    running_jobs[new_job] = 0              
    
    while len(running_jobs) > 0:
        finished_jobs = []
        for job in running_jobs:
            alive = job.is_alive()
            if not alive:
                finished_jobs.append(job)
        for job in finished_jobs:
            running_jobs.pop(job)
    
    print time_stamp(), '%s%s%d%%%s%s of jobs completed.%s' % (color.BOLD, color.BLUE, int((total_jobs - len(jobs_to_do)) / float(total_jobs) * 100), color.END, color.BOLD, color.END)


def get_genomes(genomes_file):
    """
    Returns a dictionary of strain and its [sequence, start_codon_position].
    """
    print time_stamp(), 'Reading in genomes ...'
    
    strains = []
    genomes_fasta = SeqIO.parse(open(os.path.realpath(genomes_file)), 'fasta')
    
    nucleotides = set('ACGTU')
    
    for genome in genomes_fasta:
        
        # skip genomes with ambiguous bases
        complete = set(str(genome.seq)).difference(nucleotides)
        if len(complete) > 0:
            print '\nSkipping FASTA entry', genome.description, 'as it contains inappropriate nucleotides: %s' % str(', '.join(list(complete))), ''
            continue
        
        # store genomes for processing
        strain_name = genome.description.split('|')[3].split('.')[0]
        strains.append(Strain(strain_name, str(genome.seq)))
    return strains


def run_fortran_tfold(strain, position, sequence, directory, tfold_sample_size, structure=False):
    """
    Run fortran Tfold. It's in the name.
    """
    file_name = '%s%s_%d_%d' % (directory, strain, len(sequence), position)
    
    temp_file = open(file_name, 'w')
    if structure:
        temp_file.write('%s\n%s' % (sequence, structure))
        temp_file.close()
        bashCommand = './Tfold.x -i %s -o %s.folds -T %f -E' % (file_name, file_name, 273.15+36.7)
    else:
        temp_file.write('>1\n%s' % (sequence))
        temp_file.close()
        bashCommand = './Tfold.x -p %d -i %s -o %s.folds -s %d -T %f' % (tfold_sample_size, file_name, file_name, 123456789, 273.15+36.7)
    
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    process.communicate()
    
    os.remove(file_name)
    
    return file_name + '.folds'


def __get_symmetric_structure(start, end, structure):
    """
    Finds the symmetric structure and positions based on the position of the apical loop.
    """
    structure = list(structure)
    left_bp_count = structure.count('(')
    right_bp_count = structure.count(')')
    
    while left_bp_count > right_bp_count or structure[0] == '.':
        structure.pop(0)
        left_bp_count = structure.count('(')
        start += 1
            
    while right_bp_count > left_bp_count or structure[-1] == '.':
        structure.pop()
        right_bp_count = structure.count(')')
        end -= 1
    
    return start, end, ''.join(structure)
    

def __get_symmetric_structure_energy(sl, genome, directory, strain):
    """
    Returns the energy of the symmetric structure.
    """
    folds_file = run_fortran_tfold(strain, sl.start, genome[sl.start:sl.end], directory, 1, sl.structure)
    energy = float(open(folds_file).read().split('\n')[1].split()[1])
    os.remove(folds_file)
    return energy


def _get_SLs(SLs, fold, genome, position, motif_pattern, directory, strain, tfold_sample_size):
    """
    Returns the SLs with the update SL info.
    """
    for match in re.finditer('[\(\.]+\(\.+?\)[\)\.]+', fold):
        start = match.start()
        end = match.end()
        start, end, structure = __get_symmetric_structure(start, end, match.group())
        
        # get initial stem loop object
        sl = SL(structure, start + position, end + position)
        motif, apical_start = sl.get_apical_loop_info(genome)
        
        # get basepairs
        bulge_graph = cgb.BulgeGraph()
        bulge_graph.from_dotbracket(structure)
        for pos in range(1, len(structure)+1):
            pos_prime = bulge_graph.pairing_partner(pos)
            if pos < pos_prime and pos and pos_prime:
                sl.basepairs.add((pos+position+start, pos_prime+position+start))
        
        # check if the apical loop has the motif
        match_motif = [m for m in re.finditer(motif_pattern, motif)]
        
        if match_motif:
            sl.motif_position = match_motif[0].start() + apical_start
        
        matches = filter(sl.check_if_same_SL, SLs)
        
        # update SLs list with SL entry
        if not matches:
            # calculate the folding energy of the symmetric structure
            sl.min_free_energy = __get_symmetric_structure_energy(sl, genome, directory, strain)
            SLs.append(sl)
        else:
            matches[0].fold_count += sl.fold_count
            matches[0].fold_stability += sl.fold_stability
    
    return SLs


def _update_shared_SL_structures_and_stability(SLs, min_PS_stability_preference):
    """
    Determines and records the strctures and their stability that is shared by multiple SLs, i.e. SLs that have the same minimum structure with some variation in base pairs.
    """
    # reset previously computed variables 
    for parent_sl in SLs:
        parent_sl.fold_stability = parent_sl.fold_count
        parent_sl.structures = []
        parent_sl.structures_common_basepairs = []
        parent_sl.SLs_in_common = []
    
    # update SL stability and neighbour map
    for parent_sl in SLs:
        # update parent fold stability by the number of its basepairs shared
        if parent_sl.motif_position:
            parent_sl.fold_stability *= min_PS_stability_preference
        
        for child_sl in SLs:
            if parent_sl.check_if_same_SL(child_sl):
                continue
            
            common_basepairs = parent_sl.basepairs.intersection(child_sl.basepairs)
            if common_basepairs:
                # update parent fold stability with related structure fold stabilities as described above
                child_score = (len(common_basepairs) / len(parent_sl.basepairs)) * child_sl.fold_count
                
                if child_sl.motif_position:
                    child_score *= min_PS_stability_preference
                
                parent_sl.fold_stability += child_score
                parent_sl.structures.append(child_sl.structure)
                parent_sl.structures_common_basepairs.append(common_basepairs)
                parent_sl.SLs_in_common.append(child_sl)


def _write_to_file(f, sl, strain, window_position, window_size, sequence_fragment):
    """
    Writes data to file.
    """
    # write sequence information
    matched_RM = str(len(set([s.motif_position for s in sl.SLs_in_common if s.motif_position])))
    f.write('%s,%d,%s,%s,%d,%d,%s,%s,%s,%s,%s,%s\n' % (strain.name, len(strain.SLs), matched_RM, 'N/A', window_position+1, window_position+window_size-1, 'N/A', sequence_fragment, 'Sequence', 'N/A', 'N/A', 'N/A'))
    
    # write parent SL first
    motif_position = sl.motif_position
    if motif_position:
        motif_position = str(motif_position+1)
    else:
        motif_position = 'None'
    
    f.write('%s,%d,%s,%s,%d,%d,%d,%s,%s,%d,%d,%3.2f\n' % (strain.name, len(strain.SLs), motif_position, sl.get_apical_loop_info(strain.sequence)[0], sl.start+1, sl.end+1, sl.get_apical_loop_info(strain.sequence)[1]+1, ' '*(sl.start-window_position) + sl.structure, 'Top Hit', sl.fold_count, sl.fold_stability, sl.min_free_energy))
    
    # write grouped SLs sorted by a criteria (see: lambda function)
    for grouped_sl in sorted(sl.SLs_in_common, key=lambda x: x.fold_stability, reverse=True):
        motif_position = grouped_sl.motif_position
        if motif_position:
            motif_position = str(motif_position+1)
        else:
            motif_position = 'None'
        
        f.write('%s,%d,%s,%s,%d,%d,%d,%s,%s,%d,%d,%3.2f\n' % (strain.name, len(strain.SLs), motif_position, grouped_sl.get_apical_loop_info(strain.sequence)[0], grouped_sl.start+1, grouped_sl.end+1, grouped_sl.get_apical_loop_info(strain.sequence)[1]+1, ' '*(grouped_sl.start-window_position) + grouped_sl.structure, 'Structure', grouped_sl.fold_count, grouped_sl.fold_stability, grouped_sl.min_free_energy))


def _write_genome_fold(strain, f_name):
    """
    Writes out full genome with the corresponding top fold.
    """
    # all SLs
    f = open(f_name.replace('.csv', '.fasta'), 'w')
    fold = ['.' for i in range(len(strain.sequence))]
    for sl in strain.SLs:
        for basepair in sl.basepairs:
            fold[basepair[0]] = '('
            fold[basepair[1]] = ')'
    f.write('%s\n' % strain.sequence)
    f.write('%s\n' % ''.join(fold[1:]))
    f.close()
    
    # only SLs with a recognition motif
    f = open(f_name.replace('.csv', '_SLs_with_RM.fasta'), 'w')
    fold = ['.' for i in range(len(strain.sequence))]
    for sl in strain.SLs:
        matched_RM = len(set([s.motif_position for s in sl.SLs_in_common if s.motif_position]))
        if matched_RM:
            for basepair in sl.basepairs:
                fold[basepair[0]] = '('
                fold[basepair[1]] = ')'
    f.write('%s\n' % strain.sequence)
    f.write('%s\n' % ''.join(fold[1:]))
    f.close()
    
    # only SLs with a high score
    f = open(f_name.replace('.csv', '_high_score.fasta'), 'w')
    fold = ['.' for i in range(len(strain.sequence))]
    for sl in strain.SLs:
        if sl.fold_stability > 100:
            for basepair in sl.basepairs:
                fold[basepair[0]] = '('
                fold[basepair[1]] = ')'
    f.write('%s\n' % strain.sequence)
    f.write('%s\n' % ''.join(fold[1:]))
    f.close()


def _fold(strain, motif_pattern, max_window_size, min_PS_stability_preference, tfold_sample_size, back_step_size, temp_directory, output_directory, min_distance_between_SLs):
    """
    Subprocess of a threaded folding method.
    """
    print time_stamp(), 'Started folding strain %s%s%s' % (color.BOLD, strain.name, color.END)
    
    f_name = '%sCTFolds_%s_maxWindow_%d_backStepSize_%d_StabPrefOfPS_a-%d-b-%d_MinSLdistance_%d.csv' % (output_directory, strain.name, max_window_size, back_step_size, min_PS_stability_preference[0], min_PS_stability_preference[1], min_distance_between_SLs)
    f = open(f_name, 'w')
    f.write('Strain,SL Number,RM Position (RM Count),Motif,Start,End,Loop Position,Structure,Type,Fold Count,Calculated Stability,Min. Free Energy\n')
    
    window_size = 7
    window_position = 0
    genome_size = len(strain.sequence)
    SLs = []            # holds SLs of the current folding window; reset only when a new SL has been found (appended)
    
    PS_counter = 0
    while window_position + window_size < genome_size:
    
        sequence_fragment = strain.sequence[window_position:window_position+window_size]
        
        folds_file = run_fortran_tfold(strain.name, window_position, sequence_fragment, temp_directory, tfold_sample_size)
        folds = open(folds_file).read().split('\n')[1:-1]
        os.remove(folds_file)
        
        for fold in folds:
            SLs = _get_SLs(SLs, fold.split()[0], strain.sequence, window_position, motif_pattern, temp_directory, strain.name, tfold_sample_size)
        
        # increase the window size for the next iteration
        window_size += 1
        
        # decide which SL to fix permanently, then reset working objects
        if window_size > max_window_size and SLs:
            
            # update based on a concentration dependent event and followed by 4 independent of concentration binding events
            if PS_counter % 5 == 0:
                _update_shared_SL_structures_and_stability(SLs, min_PS_stability_preference[0])
            else:
                _update_shared_SL_structures_and_stability(SLs, min_PS_stability_preference[1])
            
            # determine if a compatible SL can be appended within the search window
            appended = False
            for sl in sorted(SLs, key=lambda x: x.fold_stability, reverse=True):
                if len(strain.SLs) > 0:
                    # check if the free energy of the new SL is above the set criteria
                    if sl.min_free_energy >= 0.:
                        continue
                    
                    # avoid very low scoring SLs to allow for larger folds to form
                    if sl.fold_stability < 10:
                        continue
                    
                    # check if new SL structure does not overlap with the previous fixed SL sructure; also check that another SL can be added
                    if sl.start <= strain.SLs[-1].end or appended:
                        continue
                    
                    # check if previously fixed SL is the same as the new SL
                    if sl.check_if_same_SL(strain.SLs[-1]):
                        continue
                    
                    # check that there is enough space between two adjacent PSs
                    if strain.SLs[-1].motif_position and sl.motif_position:
                        if strain.SLs[-1].motif_position + min_distance_between_SLs > sl.motif_position:
                            continue
                
                strain.SLs.append(sl)
                appended = True
                _write_to_file(f, sl, strain, window_position, window_size, sequence_fragment)
            
            # if a compatible SL is found, reset the settings
            if appended:
                window_position = strain.SLs[-1].end - back_step_size
                if window_position < 0:
                    window_position = 1
                window_size = 7
                SLs = []
                
                # update PS counter
                if strain.SLs[-1].motif_position:
                    PS_counter += 1
    
    _write_genome_fold(strain, f_name)
    
    print time_stamp(), 'Folding completed for strain %s%s%s' % (color.BOLD, strain.name, color.END)


def multithread_fold(strains, motif_pattern, max_window_sizes, min_PS_stability_preferences, tfold_sample_sizes, back_step_sizes, max_CPUs, min_distance_between_SLs):
    """
    Initiate co-folding sequence - multithreaded.
    """
    jobs_to_do = []
    temp_directories = []
    for max_window_size in max_window_sizes:
        for min_distance_between_SL in min_distance_between_SLs:
            for back_step_size in back_step_sizes:
                for min_PS_stability_preference in min_PS_stability_preferences:
                    temp_directory = os.getcwd() + os.sep + 'temp_folding_jobs_%d_%d_%d_%d_%d' % (max_window_size, back_step_size, min_PS_stability_preference[0], min_PS_stability_preference[1], min_distance_between_SL) + os.sep
                    temp_directories.append(temp_directory)
                    try:
                        os.makedirs(temp_directory)
                    except OSError:
                        pass
                    
                    output_directory = os.getcwd() + os.sep + 'ParamSampling_results' + os.sep
                    try:
                        os.makedirs(output_directory)
                    except OSError:
                        pass
                    
                    for strain in strains:
                        job = multiprocessing.Process(target=_fold, args=(strain, motif_pattern, max_window_size, min_PS_stability_preference, tfold_sample_size, back_step_size, temp_directory, output_directory, min_distance_between_SL))
                        jobs_to_do.append(job)
    
    print time_stamp(), '%s%sSubmitting genomes for co-transcriptional folding ...%s' % (color.BLUE, color.BOLD, color.END)
    
    if len(jobs_to_do) > 0:
        run_jobs(jobs_to_do, max_CPUs)
    else:
        print '\nDid not find any jobs to do! Stopping the application.'
    print time_stamp(), '%s%sRemoving temporary files.%s\n' % (color.BLUE, color.BLUE, color.END)
    
    for directory in temp_directories:
        os.rmdir(directory)


if __name__ == '__main__':
    # Parameters
    max_CPUs = 70                                   # max CPUs to utilise
    motif_pattern = 'G.T..'                         # RegEx to identify the SLs with a motif for binding
    back_step_size = [0]                            # distance the window is shifted back from the end of the last fixed SL
    max_window_size = [120]                         # maximum distance the folding window can grow before fixing structures in place; window size resets from last fixed SL
    min_PS_stability_preference = [(1,1), (4,4)]    # the motif binding factor score contribution 1...10 where 1 is no contribution and 10 is 10x contribution to the score - provided as a tuple of the same value (variation can be used for probing other settings, not yet described for the reader)
    min_distance_between_SLs = [26]                 # the minimum nucleotide distance between the start positions of the recognition motifs 
    tfold_sample_size = 10                          # Tfold sample size setting
    
    # Run
    strains = get_genomes(sys.argv[1])              # genomes FASTA file
    multithread_fold(strains, motif_pattern, max_window_size, min_PS_stability_preference, tfold_sample_size, back_step_size, max_CPUs, min_distance_between_SLs)