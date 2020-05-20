from collections import defaultdict

def create_profile_hmm(alignment, threshold, alphabet):
    star_alignment = set()
    num_patterns = len(alignment)
    for idx, _ in enumerate(alignment[0]):
        deletions = sum(1 for pattern in alignment
                            if pattern[idx] == '-')
        if deletions <= threshold * num_patterns:
            star_alignment.add(idx)

    paths = []
    for pattern in alignment:
        path = ['S']
        curr = 0
        for idx, ltr in enumerate(pattern):
            if idx in star_alignment:
                curr += 1
                if ltr != '-':
                    path.append('M'+str(curr)) # Letter emitted
                else:
                    path.append('D'+str(curr)) # Deletion

            else:
                if ltr != '-':
                    path.append('I'+str(curr)) # Insertion
        path.append('E')
        paths.append(path)




if __name__ == '__main__':
    threshold = 0.289
    alphabet = 'ABCDE'
    multiple_alignment = ['EBA', 'E-D', 'EB-', 'EED', 'EBD', 'EBE', 'E-D', 'E-D']
    create_profile_hmm(multiple_alignment, threshold, multiple_alignment)