__author__ = "Keenan Manpearl"
__date__ = "2023/03/06"

def create_graph(kmers):
    graph = []
    for k1 in range(len(kmers)):
        for k2 in range(len(kmers)-1):
            if kmers[k1].sufix == kmers[k2].prefix:
                graph.append([kmers[k1].sequence, kmers[k2].sequence])
    return(graph)