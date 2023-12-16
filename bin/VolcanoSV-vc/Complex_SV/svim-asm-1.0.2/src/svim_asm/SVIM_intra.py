from __future__ import print_function

import sys

from svim_asm.SVCandidate import CandidateDeletion, CandidateInsertion


def analyze_cigar_indel(tuples, min_length):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos_ref = 0
    pos_read = 0
    indels = []
    for operation, length in tuples:
        if operation == 0:                     # alignment match
            pos_ref += length
            pos_read += length
        elif operation == 1:                   # insertion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "INS"))
            pos_read += length
        elif operation == 2:                   # deletion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "DEL"))
            pos_ref += length
        elif operation == 4:                   # soft clip
            pos_read += length
        elif operation == 7 or operation == 8:        # match or mismatch
            pos_ref += length
            pos_read += length
    return indels


def analyze_alignment_indel(alignment, bam, query_name, options):
    sv_candidates = []
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start
    indels = analyze_cigar_indel(alignment.cigartuples, options.min_sv_size)
    for pos_ref, pos_read, length, typ in indels:
        if typ == "DEL":
            sv_candidates.append(CandidateDeletion(ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, [query_name], bam))
        elif typ == "INS":
            insertion_seq = alignment.query_sequence[pos_read:pos_read+length]
            sv_candidates.append(CandidateInsertion(ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, [query_name], insertion_seq, bam))
    return sv_candidates


