def suffix_array(text):
    idx_to_suffix = {idx: text[idx:] for idx, _ in enumerate(text)}
    arr = list(sorted(idx_to_suffix, key=idx_to_suffix.get))
    return arr

def burrows_wheeler_transform(text):
    cyclic_shifts = [text[shift:] + text[:shift] for shift in range(len(text))]
    lexic_shifts = sorted(cyclic_shifts)
    bwt = ''.join([pattern[-1] for pattern in lexic_shifts])
    return bwt

def invert_bwt(bwt):
    pass