def suffix_array(text):
    idx_to_suffix = {idx: text[idx:] for idx, _ in enumerate(text)}
    arr = list(sorted(idx_to_suffix, key=idx_to_suffix.get))
    return arr

def burrows_wheeler_transform(text):
    pass
