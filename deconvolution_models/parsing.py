

from epiread_tools.epireadToBedgraph import EpiRunner
from epiread_tools.em_utils import Mapper
from epiread_tools.naming_conventions import *
import numpy as np
import pandas as pd
import pickle

#%%
'''
Read mixture file to list of methylation matrices (one np array per interval)
save to .npy file


- read bed intervals (TIMs) for GenomicInterval list
- split_intervals_to_chromosomes
- parse epireads
- save


mixture should be in epiread format 
atlas for celfie/celfie+ should be bedgraph per cpg, regions file  

for each chromosome:
    parse intervals to matrix
    split large matrix to list per region
    save stats on each mat in list
    
epireadToBedgraph has only one mat
    
'''
class MixtureReader(EpiRunner):

    def __init__(self, regions, mixture, CpG_file, outfile, epiformat="old_epiread", header=False):
        super().__init__(regions, CpG_file, [mixture], outfile, epiformat, header=header, bedfile=True)

    def sum_regions(self): #celfie output, one value per TIM
        '''
        parse mixture and return sum per region
        :return:
        '''
        self.parse_multiple_chromosomes()
        self.calc_coverage()
        self.calc_methylation()
        return [np.sum(x) for x in self.methylation], [np.sum(x) for x in self.coverage]

    def get_matrices(self): #for celfie-plus
        '''
        parse mixture, retrieve raw reads
        :return:
        '''
        self.parse_multiple_chromosomes()
        return self.matrices



class AtlasReader(EpiRunner):
    '''
    a note on strandedness:
    All positions must match the cpg file. To process
    strand separately, make sure the CoG file includes separate coordinates
    e.g. chrN:200-201, chrN:201-202. This has to be the same file the epireads
    are aligned to
    '''
    def __init__(self, tims, atlas_file, cpg_file, outfile):
        '''

        :param tims: Tissue informative markers, bed format
        :param atlas_file: methylation values per cpg
        :param cpg_file: coordinates o cpgs in genome assembly
        :param outfile: path to putput file
        '''
        self.atlas = atlas_file
        super().__init__(tims, cpg_file, [atlas_file], outfile)
        self.init_beta()

    def init_beta(self):
        df = pd.read_csv(self.atlas, sep="\t")
        vals = df.iloc[:,3:].values
        beta = vals[:,::2]/vals[:,1::2]
        self.beta = beta.T
        self.atlas_chrom = df["CHROM"].values
        self.atlas_start = df["START"].values  #remove if 0 based

    def align_beta(self, chrom, mapper):

        mat = np.ones(shape=(self.beta.shape[0], mapper.max_cpgs)) #cell types by cpgs
        mat.fill(np.nan)
        chrom_filter = (self.atlas_chrom==chrom)
        cpgs = self.atlas_start[chrom_filter]
        for i, cpg in enumerate(cpgs):
            if cpg in mapper.abs_to_rel and mapper.abs_to_rel[cpg] in mapper.all_rel:
                mat[:,mapper.abs_to_ind(cpg)] = self.beta[:,chrom_filter][:,i]
            elif (cpg in mapper.abs_to_rel and mapper.abs_to_rel[cpg] not in mapper.all_rel):
                #not in intervals, skip
                print("not in intervals", cpg)
                pass
            else:
                raise KeyError(cpg, "atlas start doesn't match cpg file")
        return mat

    def parse_one_chromosome(self,chrom, intervals):
        intervals = sorted(intervals, key=lambda x: x.start)
        self.interval_order.extend([str(x) for x in intervals])
        mapper = Mapper(chrom, intervals, [], self.cpg_locations, False)  # init mapping
        window_list = mapper.get_ind_intervals(intervals)
        mat = self.align_beta(chrom, mapper)
        for start, end in window_list:
            slice = mat[:,start:end]
            self.matrices.append(slice)

    def write_plus_output(self):
        data = dict(zip(self.interval_order, self.matrices))
        sorted_data = [data[str(x)] for x in self.genomic_intervals]  # sort by input order
        with open(self.outfile, 'wb') as outfile:
            pickle.dump(sorted_data, outfile, protocol=pickle.HIGHEST_PROTOCOL)


