# Trans-chromosomal frequency score adjusted by copy number (cnadjTIF)


Description:  An R script that takes genome-wide matrix (cis can be NA) and copy number prediction. For each bin in the matrix, 
the average trans-chromosomal interaction frequency (TIF) is adjusted by copy number. 


Main code: cnadj.R  


Inputs included here: 
1) HF2354_i50000_rawPET0.mat (N x N matrix). 
2) hg19_i50000.ind (N rows of bin info).
3) HF2354.subtractBlacklist.CN.bdg 
4) z3_encode_gap.blacklist.bed (blacklist, a concatenation of UCSC hg19 gaps, ENCODE blacklist, and ChIA-PET greylist).
5) ecDNA_regions.txt (Bed format ecDNA regions).


Additional code: 
clusters2matrix.hg19.py, prints out aggregate interaction matrix from pair tag reads with interaction counts; works directly with ChIA-PET Utilities ([CPU](https://github.com/cheehongsg/CPU/wiki)).

Reference:   
Zhu Y, Gujar AD, Wong CH, Tjong H, Ngan CY, Gong L, Chen YA, Kim H, Liu J, Li M, Mil-Homens A, Maurya R, Kuhlberg C, Sun F, Yi E, deCarvalho AC, Ruan Y, Verhaak RGW, Wei CL. "Oncogenic extrachromosomal DNA functions as mobile enhancers to globally amplify chromosomal transcription." Cancer Cell. 2021 May 10;39(5):694-707.e7.   
doi: 10.1016/j.ccell.2021.03.006. Epub 2021 Apr 8. PMID: 33836152; PMCID: PMC8119378.
