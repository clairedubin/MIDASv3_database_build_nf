
def has_ambiguous_bases(sequence):
    # Check if sequence contains lower-case letters, which usually indicate soft-masked bases
    ambiguous_bases = ['N', 'X', 'n', 'x']
    return any(base in ambiguous_bases for base in sequence)

