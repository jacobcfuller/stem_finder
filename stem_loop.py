# Rule 1: substring length >= 12 and <=100
# Rule 2: A ‘GNRA’ motif is at the center of the subsequence
# Rule 3: The subsequence contains at least 4 pairs of complementary
#         nucleotides.
# Rule 4: Each half of each complementary pair must be located equidistant from
#         the center of the ‘GNRA’ motif as shown in the diagram of the hairpin
#         structure above.
# Rule 5: Both nucleotides of one of the complementary pairs must be a distance
#         of <= 2 from either end of the ‘GNRA’ motif
# Rule 6: At least 70% of nucleotides in a complementarity region must be
#         participating in a pairing interaction.


def import_seq(input_txt_file):
    """Return the input sequence with only dna sequence from question3_input

    Parameters
    ----------
    input_txt_file : .txt file

    Returns
    -------
    str
        Sequence from question3_input.txt with non-nucleotides removed

    Notes
    -----
    Only works with the provided question3_input.txt
    """
    with open(input_txt_file, 'r') as input_seq:
        data = input_seq.read().replace('\n', '')
    return data.split(":")[1]


def complement(s):
    """Return the complement of the given DNA or RNA string

    Parameters
    ----------
    s : str
        DNA or RNA string

    Returns
    -------
    str
        Complemented DNA or RNA string

    Notes
    -----
    Degenerate bases (W, Y, N, etc) are allowed

    Raises
    ------
    ValueError
        Unknown character(s) passed
    """
    known_bases = set('ATUCGNatucgnKRSBDMYWVHkm')
    unknown = set(s) - known_bases
    if len(unknown) > 0:
        raise ValueError('Unknown bases passed: %s' % ', '.join(unknown))

    complement = str.maketrans('ATUCGNatucgnKRSBDMYWVHkm',
                               'TAAGCNtaagcnMYSVHKRWBDmk')
    return s.translate(complement)


def compl_ratio(seq):
    """Given a DNA sequence, remove the 'GNRA motif, then find the
    complementarity ratio in a potential stem

    Parameters
    ----------
    seq : str
        DNA string with 'GNRA' motif in the center

    Returns
    -------
    float
        Ratio of complementary/total length

    Notes
    -----
    seq must have an even # of nucleotides

    Raises
    ------
    ValueError
        Seq length must be even
        Unknown character(s) passed
    """
    seq_len = len(seq)
    if seq_len % 2 != 0:
        raise ValueError('Seq length must be even')
    front = seq[:int(seq_len/2 - 2)]
    end = seq[int(seq_len/2 + 2):]
    revc = complement(front[::-1])
    same = "".join(n for idx, n in enumerate(revc) if n == end[idx])

    ratio = len(same)/len(front)

    return ratio


def find_motif(seq_str):
    """Finds all "GNRA" motifs where at least one complementary pair exists
    within 2 bp (rule 5). Returns list of location and motif tuple

    Parameters
    ----------
    seq_str : str
        input sequence

    Returns
    -------
    motif : list
        [(locus int, motif str)]
        locus is location of 1st 'G' in GNRA motif
        motif is seq with GNRA in center, and 2 bp on each side of GNRA

    Raises
    ------
    ValueError
        Unknown character(s) passed
    """
    seq = seq_str

    motif = []
    # motif can't be first 4 in seq (or last 4)
    i = 3
    while i <= (len(seq) - 6):
        # get potential motif & surrounding 2 bp on each side
        win = seq[(i - 2):(i + 6)]
        # see if there is a complementary pair within 2 bp of motif spot
        if complement(win[0]) == win[7] or complement(win[1]) == win[6]:
            # now see if motif is there
            if win[2] == 'G' and win[5] == 'A' and win[4] in ('G', 'A'):
                motif.append((i, win))
        i += 1
    return motif


def search(motif_list, seq):
    """Return list of stem loop sequences according to rules. Uses motif list
    location, then starts at potential biggest stem loop and narrows down
    smaller until >= 70 percent complementarity.

    Parameters
    ----------
    list - motif_list
        List of (locus, motif) tuples from find_motif
    str - seq
        sequence from question3_input.txt

    Returns
    -------
    list - loop_seqs
        list stem loop sequences

    """
    loop_seqs = []
    for motif in motif_list:
        pos = motif[0]
        # if far enough into sequence to test for max size (100)
        if 48 <= pos < (len(seq)-52):
            i = 48
            while i > 3:
                seq_chunk = seq[pos-i:pos+i+4]
                if seq_chunk[0] == complement(seq_chunk[-1]):
                    ratio = compl_ratio(seq_chunk)
                    if ratio >= .70:
                        loop_seqs.append(seq_chunk)
                        i = 0
                i -= 1
        # if at start or end of input seq (stem loop smaller than 100)
        else:
            i = pos
            while i > 3:
                seq_chunk = seq[pos-i:pos+i+4]
                if seq_chunk[0] == complement(seq_chunk[-1]):
                    ratio = compl_ratio(seq_chunk)
                    if ratio >= .70:
                        loop_seqs.append(seq_chunk)
                        i = 0
                i -= 1

    return loop_seqs


# To test. Will print loop sequences to console
if __name__ == '__main__':
    input_path = input("Please enter path to question3_input.txt: ")
    seq = import_seq(input_path)
    motif = find_motif(seq)
    loops = search(motif, seq)

    for loop in loops:
        print(loop + '\n')
