import networkx as nx
import numpy as np
from collections import Counter

FOLDER = "Np_216_random_p_0.05-1/p_0.1/1/"

P_data = np.loadtxt(FOLDER + 'positionP.csv', delimiter=',')
P_connections = P_data[:, -3:].astype(np.int)  # P, Pt, Pt

G = {f'P{i}': [f'P{j}', f'Pt{k}', f'Pt{l}']
     for i, (j, k, l) in enumerate(P_connections, 1)}

graph = nx.Graph(G)
digraph = nx.DiGraph(graph)
cycles = nx.simple_cycles(digraph)
# print(cycles)
# for cycle in cycles:
#     if len(cycle) == 5:
#         print(cycle)

cycle_lens = Counter(map(len, cycles))
for l in [4, 5, 6]:
    print(l, ':', cycle_lens[l] / 2)
