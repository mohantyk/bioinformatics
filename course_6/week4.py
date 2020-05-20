import numpy as np
import pandas as pd

def viterbi(emitted, alphabet, states, transitions, emissions):
    transition_prob = pd.DataFrame( data=transitions, columns=list(states), index=list(states) )
    emission_prob = pd.DataFrame( data=emissions, columns=list(alphabet), index=list(states))

    graph = np.zeros((len(states), len(emitted)), dtype=float)
    backtrack = np.empty((len(states), len(emitted)), dtype=str)
    for idx, ltr in enumerate(emitted):
        if idx == 0:
            graph[:, idx] = np.log(emission_prob.loc[:, ltr]*1/len(states))
            continue
        for s_idx, state in enumerate(states):
            all_links = graph[:,idx-1] + np.log(transition_prob.loc[:,state].values) + np.log(emission_prob.loc[state, ltr])
            graph[s_idx, idx] = np.max(all_links)
            backtrack[s_idx, idx] = states[np.argmax(all_links)]

    path = []
    s_idx = np.argmax(graph[:, -1])
    state = states[s_idx]
    path.append(state)
    for idx in reversed(range(len(emitted))):
        s_idx = states.index(state)
        state = backtrack[s_idx, idx]
        path.append(state)

    final = ''.join(path[::-1])
    return final


def outcome_likelihood(emitted, alphabet, states, transitions, emissions):
    transition_prob = pd.DataFrame( data=transitions, columns=list(states), index=list(states) )
    emission_prob = pd.DataFrame( data=emissions, columns=list(alphabet), index=list(states))

    graph = np.zeros((len(states), len(emitted)), dtype=float)
    for idx, ltr in enumerate(emitted):
        if idx == 0:
            graph[:, idx] = emission_prob.loc[:, ltr]*(1/len(states))
            continue
        for s_idx, state in enumerate(states):
            graph[s_idx, idx] = np.dot(graph[:,idx-1], transition_prob.loc[:,state].values) * emission_prob.loc[state, ltr]
    final = np.sum(graph[:,-1])
    return final
