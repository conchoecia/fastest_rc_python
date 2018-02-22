from libc.stdlib cimport malloc

cdef int seq_len = 17
cdef char *seq_dest = <char *>malloc(seq_len + 1)
seq_dest[seq_len] = '\0'

cdef char *basemap = [ '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',   '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a' ]

def reverse_complement_c_v2(str seq):
    cdef bytes py_bytes = seq.encode('UTF-8')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    for i in range(seq_len):
        seq_dest[seq_len - i - 1] = basemap[<int>seq_src[i]]
    return seq_dest[:seq_len].decode('UTF-8')