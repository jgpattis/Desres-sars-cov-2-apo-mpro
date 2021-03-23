import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

mi_name = 'sasa_mi.npy'
exp_name = 'exposons.npy'
out_name = 'plots/exposons_try1'
o_type = 'pdf'

mi = np.load(open(mi_name,"rb"))
exp = np.load(open(exp_name,"rb"))


non_pocket = 1
pockets = 6

unique, un_inds, un_counts = np.unique(exp, return_index=True, return_counts=True)
order = np.argsort(un_counts)

mi -= mi.diagonal() * np.eye(*mi.shape)

indsn = np.where(exp == unique[order][-1])[0]
#for i in indsn:
#    mi[i] = np.nan
#    mi[:,i] = np.nan

#for j in range(10, len(order)):
#    inds = np.where(exp == unique[order][-j])[0]
#    for k in inds:
#        mi[k] = np.nan
#        mi[:,k] = np. nan

important = []
for j in range(1, 10):
    inds = np.where(exp == unique[order][-j])[0]
    important.append(inds)

imp1 = np.concatenate(important)
imp2 = set(imp1)
not_imp = []
for i in range(mi.shape[0]):
    if i not in imp2:
        not_imp.append(i)

for i in not_imp:
    mi[i] = np.nan
    mi[:,i] = np.nan

result = np.transpose(np.where(mi > 0.001))

graph = nx.DiGraph()

for i, j in result:
    graph.add_edge(i, j, weight=mi[i, j])


## layouts are more of an art than a science 

#pos = nx.spiral_layout(graph)
pos = nx.nx_agraph.graphviz_layout(graph, prog="neato")
#pos = nx.shell_layout(graph)
#pos = nx.kamada_kawai_layout(graph, scale=4)
#pos = nx.circular_layout(graph)
#pos = nx.spring_layout(graph) #, scale=4, iterations=100)
#pos = nx.spectral_layout(graph)
#pos = nx.drawing.nx_agraph.graphviz_layout(graph)

node_size = 500

node_color = [exp[i] for i in graph.nodes()]
#node_color='pomegranate'

fig, ax = plt.subplots()

nx.draw(graph, pos=pos, node_color=node_color, arrows=False,
                     node_size=node_size, ax=ax, with_labels=True)

fig.savefig('exposon_network_try12.pdf')
