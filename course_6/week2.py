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
    last_col = bwt
    first_col = ''.join(sorted(bwt))

    first_row = '$'
    last_ltr_count = 0
    for row_idx in range(1, len(bwt)-1):
        last_ltr = first_row[-1]

        # Find col idx of #last_ltr_count occurence of last_ltr in last_col
        count = 0
        for col_idx, ch in enumerate(last_col):
            if ch == last_ltr:
                if count == last_ltr_count:
                    break
                else:
                    count += 1
        nxt_ltr = first_col[col_idx]
        first_row = first_row + nxt_ltr
        #Find count of nxt_ltr
        last_ltr_count = len([ch for ch in first_col[:col_idx] if ch==nxt_ltr])

    first_row = first_row + last_col[0] # Add the last character of first row
    text = first_row[1:] + first_row[0] # Move the $ to the end
    return text


def bw_match(bw_test, pattern):
    pass

def match_patterns(bw_text, patterns):
    num_matches = [bw_match(bw_text, pattern) for pattern in patterns]
    return num_matches