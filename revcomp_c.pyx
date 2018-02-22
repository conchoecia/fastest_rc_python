from libc.stdlib cimport malloc

cdef int seq_len = 17
cdef char *seq_dest = <char *>malloc(seq_len + 1)
seq_dest[seq_len] = '\0'

def reverse_complement_c_v1(str seq):
    cdef bytes py_bytes = seq.encode('ascii')
    cdef char *seq_src = py_bytes
    cdef int i = 0
    cdef int d = 0
    for i in range(seq_len):
        d = seq_len - i - 1
        if seq_src[i] == 'A':
            seq_dest[d] = 'T'
        elif seq_src[i] == 'G':
            seq_dest[d] = 'C'
        elif seq_src[i] == 'C':
            seq_dest[d] = 'G'
        elif seq_src[i] == 'T':
            seq_dest[d] = 'A'
    return seq_dest[:seq_len].decode('ascii')