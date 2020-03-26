"""
File containing utils for dealing with FASTA files
"""

import subprocess


def parse_fasta(lines):
    """
    Return dict of {label:seq} from FASTA file

    Arguments
    ---------
    lines: list of lines, open file object

    Returns
    -------
    res: dict {label:seq}
    """
    res = {}  # result
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            label = line[1:]
            res[label] = []
        else:
            res[label].append(line)
    for k, v in res.items():
        res[k] = ''.join(v)
    return res


def alt_parse_fasta(lines):
    """
    Alternative function to parse FASTAs.
    Useful when a header appears multiple times (e.g. some alignments
    for correlated mutations)
    """
    # fasta = lines.readlines()
    # headers = fasta[::2]
    # seqs = fasta[1::2]

    # headers = [x.strip().lstrip('>') for x in headers]
    # seqs = [x.strip() for x in seqs]

    # return headers, seqs

    headers = []
    seqs = []
    seq = ''
    did_seq_finish = True

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            header = line[1:]
            headers.append(header)
            if seq:
                seqs.append(seq)
                did_seq_finish = True
        else:
            if did_seq_finish:
                seq = line
                did_seq_finish = False
            else:
                seq = ''.join([seq, line])
    seqs.append(seq)
    return headers, seqs


def write_fasta(fasta_dict, out_path,flag='w'):
    with open(out_path, flag) as f:
        for k, v in fasta_dict.items():
            f.write(''.join(['>', k, '\n']))
            f.write(''.join([v, '\n']))


def alt_write_fasta(fasta_headers, fasta_seqs, out_path, flag='w'):
    with open(out_path, flag) as f:
        for idx, val in enumerate(fasta_headers):
            f.write(''.join(['>', val, '\n']))
            f.write(''.join([fasta_seqs[idx], '\n']))


def write_paired_fastas(col1_ids, col2_ids, seqs, interm, out_path1, out_path2):
    """
    Write two FASTA files containing pairs of interacting sequences.
    Only writes heterodimers.
    """
    fasta1 = open(out_path1, 'w')
    fasta2 = open(out_path2, 'w')

    for idx, val in enumerate(col1_ids):
        if val != col2_ids[idx]:  # Ignore homodimers
            header1 = interm[val]
            seq1 = seqs[header1]

            header2 = interm[col2_ids[idx]]
            seq2 = seqs[header2]

            fasta1.write(''.join(['>', header1, '\n']))
            fasta2.write(''.join(['>', header2, '\n']))

            fasta1.write(''.join([seq1, '\n']))
            fasta2.write(''.join([seq2, '\n']))

    fasta1.close()
    fasta2.close()


def call_mafft(infile, outfile, threads=2):
    """
    Call MAFFT with default options.
    """
    try:
        with open(outfile, "w") as msa:
            subprocess.call(
                ["mafft", "--auto",  "--thread", str(threads), infile], stdout=msa)
    except:
        raise ValueError("MAFFT failed!")
