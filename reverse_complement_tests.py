#!/usr/bin/env python
import random
import timeit
import numpy as np
import os
from Bio.Seq import Seq
import string
import sys
import subprocess


#we need to build the seqpy extension if we have not yet built it
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
print(subprocess.check_output(['python','setup.py', 'build_ext', '-i']))
os.chdir(cwd)
import seqpy

cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
print(subprocess.check_output(['python','cython_setup.py', 'build_ext', '--inplace']))
os.chdir(cwd)
from revcomp_c import reverse_complement_c_v1

cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
print(subprocess.check_output(['python','cython_setup2.py', 'build_ext', '--inplace']))
os.chdir(cwd)
from revcomp_c2 import reverse_complement_c_v2


global complement
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

DNAlength = 17

#randomly generate 100k bases
int_to_basemap = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}
num_strings = 10000
num_trials = 100
random.seed(90210)
DNAstrings = ["".join([int_to_basemap[random.randint(1,4)] for i in range(DNAlength)])
              for j in range(num_strings)]
#get an idea of what the DNAstrings look like
print(DNAstrings[0:5])

def reverse_complement_naive(seq):
    this_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(this_complement.get(base, base) for base in reversed(seq))

def reverse_complement(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))

##https://stackoverflow.com/questions/19570800
def complement_fromSO(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)
def revcom_fromSO(s):
    return complement_fromSO(s[::-1])

#https://stackoverflow.com/questions/19570800
revcomplSO = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

#https://stackoverflow.com/questions/19570800
def revcomp_translateSO(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

#https://codereview.stackexchange.com/questions/151329/
def string_replace(seq):
    return(seq.upper().replace('A', 'temp').replace('T', 'A').replace('temp', 'T').replace('G', 'temp').replace('C','G').replace('temp','C')[::-1])

#https://bioinformatics.stackexchange.com/questions/3583/
global tab
if sys.version_info[0] < 3:
    tab = string.maketrans("WSATUGCYRKMBDHVNwsatugcyrkmbdhvn", "WSTAACGRYMKVHDBNwstaacgrymkvhdbn")
else:
    tab = "WSATUGCYRKMBDHVNwsatugcyrkmbdhvn".maketrans("WSATUGCYRKMBDHVNwsatugcyrkmbdhvn", "WSTAACGRYMKVHDBNwstaacgrymkvhdbn")
def reverse_complement_table_devonryan(seq):
    return seq.translate(tab)[::-1]

global tab_b
tab_b = bytes.maketrans(b"ACTG", b"TGAC")

def reverse_complement_bytes_jackaidley(seq):
    return seq.translate(tab_b)[::-1]

def reverse_complement_bytesthenstring_jackaidley(seq):
    return str(seq.translate(tab_b)[::-1])


bl= []
for n in range(num_trials):
    tic=timeit.default_timer()
    rcs = [reverse_complement_naive(seq) for seq in DNAstrings]
    toc=timeit.default_timer()
    bl.append(toc - tic)
baseline = np.mean(bl)


namefunc = {"naive implementation": reverse_complement_naive,
            "global dict implementation": reverse_complement,
            "biopython seq then rc": Seq,
            "revcom from SO": revcom_fromSO,
            "lambda from SO": revcomplSO,
            "revcomp_translateSO": revcomp_translateSO,
            "string_replace": string_replace,
            "devonryan string": reverse_complement_table_devonryan,
            "jackaidley bytes": reverse_complement_bytes_jackaidley,
            "jackaidley bytesstring": reverse_complement_bytesthenstring_jackaidley,
            "user172818 seqpy.c": seqpy.revcomp,
            "alexreynolds Cython implementation (v1)": reverse_complement_c_v1,
            "alexreynolds Cython implementation (v2)": reverse_complement_c_v2}

for function_name in sorted(namefunc):
    times_list = []
    for n in range(num_trials):
        func = namefunc[function_name]
        if function_name == "biopython seq then rc":
            tic=timeit.default_timer()
            rcs = [func(seq).reverse_complement() for seq in DNAstrings]
        else:
            tic=timeit.default_timer()
            rcs = [func(seq) for seq in DNAstrings]
        toc=timeit.default_timer()
        times_list.append(toc-tic)
    walltime = np.mean(times_list)
    print("""{}
       {:.5f}s total,
       {:.1f} strings per second
       {:.1f}% increase over baseline""".format(
           function_name,
           walltime,
           num_strings/walltime,
           100- ((walltime/baseline)*100)   ))
