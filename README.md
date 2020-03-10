# Trans-chromosomal frequency score adjusted by copy number (cnadjTIF)


Description:  An R script that takes genome-wide matrix (cis can be NA) and copy number prediction. For each bin in the matrix, 
the average trans-chromosomal interaction frequency (TIF) is adjusted by copy number. 


Main code: cnadj.R  


Inputs: 
1) HF2354_i50000_rawPET0.mat (N x N matrix). 
2) hg19_i50000.ind (N rows of bin info).
3) HF2354.subtractBlacklist.CN.bdg (copy number prediction from readDepth https://github.com/chrisamiller/readDepth).
4) ecDNA_regions.txt (Bed format ecDNA regions).
5) z3_encode_gap.blacklist.bed (blacklist, a concatenation of UCSC hg19 gaps, ENCODE blacklist, and ChIA-PET greylist).


Additional code: clusters2matrix.hg19.py that prints out aggregate interaction matrix from pair tag reads with interaction counts.

