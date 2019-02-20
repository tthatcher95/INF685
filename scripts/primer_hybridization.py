from itertools import groupby, compress


def dimerization(primer_1, primer_2):
    """
    Returns a dictionary with the highest total dimerization score and
    the score of the longest run (weighted by C-G vs. A-T pairings).
    The values are each the maximum (worst-case) possiblities
    among all alignments of the primers.
    primer_1 (str, Seq): First primer in 5' to 3' direction
    primer_2 (str, Seq): Second primer reverse_complemented
                         in 3' to 5' direction.
    """
    p1 = primer_1 + ("*" * (len(primer_2) - 1))
    p2 = "*" * (len(primer_1) - 1) + primer_2[::-1]
    max_alignment = []
    max_run = []
    while len(p1):
        filter_ = [0 if "*" in primers else 1 for primers in zip(p1, p2)]
        total, run, = dimerization_worker(
            compress(p1, filter_),
            compress(p2, filter_))
        max_run.append(run)
        max_alignment.append(total)
        p1 = p1[:-1]
        p2 = p2[1:]    
    return {
        'total': max(max_alignment),
        'run': max(max_run)}


HYBRID_SCORES = {'A': 2, 'T': 2, 'C':4, 'G': 4}


def dimerization_worker(primer_1, primer_2):
    """
    Returns the total score of complementary bases and the score of
    the longest run of complementary bases.
    """
    complementary = [
        HYBRID_SCORES[s1] if s1 == s2 else 0 for s1, s2 in zip(
            primer_1, primer_2)]
    total = sum(complementary)
    complementary_runs = [list(comp) for run, comp in groupby(
        complementary, key=lambda x: x > 0)]
    max_run = sum(max(complementary_runs, key=sum))
    return (total, max_run)


