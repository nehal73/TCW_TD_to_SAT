# coding=utf-8
# coding=utf-8
import argparse
import os,sys
import networkx as nx
from random import choice
import subprocess
import time
from networkx.drawing.nx_agraph import *
import networkx as nx
# import matplotlib.pyplot as plt
import numpy as np
A = nx.Graph()


def counter(s,n_var):
    lines = s.split('\n')
    n_clauses=0
    encoding = ''
    for s in lines:
        s=s.split()
        if s==[]:
            continue
        count = int(s.pop())
        operator = s.pop()
        if count==0:
            for u in s:
                encoding+='%i 0\n'%(-1*int(u))
                n_clauses+=1
            continue
        counter_var = range(n_var + 1, n_var + count * len(s) + 1)
        n_var += count * len(s)
        counter_var = np.reshape(counter_var, (len(s), count))
        for u in xrange(len(s)):
            encoding += '%i %i 0\n' % (-1 * int(s[u]), counter_var[u][0])
            n_clauses += 1
        for u in xrange(len(s) - 1):
            for i in xrange(count):
                encoding += '-%i %i 0\n' % (counter_var[u][i], counter_var[u + 1][i])
                n_clauses += 1
        for u in xrange(count - 1, len(s)):
            encoding += '%i -%i 0\n' % (-1 * int(s[u]), counter_var[u - 1][count - 1])
            n_clauses += 1
        for u in xrange(len(s) - 1):
            for i in xrange(count - 1):
                encoding += '%i -%i %i 0\n' % (-1 * int(s[u + 1]), counter_var[u][i], counter_var[u + 1][i + 1])
                n_clauses += 1
    return encoding,n_var,n_clauses


def construction(g, s, N):
    if N == {s}:
        return
    t = choice(list(N - {s}))
    x, S = nx.minimum_cut(g, s, t)
    A.add_edge(s, t, capacity=x)
    construction(g, s, N.intersection(S[0]))
    construction(g, t, N.intersection(S[1]))


def reduce_graph(g):
    # show_graph(g,3)
    nx.set_edge_attributes(g, 'capacity', 1)
    g = nx.convert_node_labels_to_integers(g, 0)
    g1 = nx.Graph()
    for edge in g.edges():
        if edge in g1.edges():
            g1.edge[edge[0]][edge[1]]['capacity'] += 1
        else:
            g1.add_edge(edge[0], edge[1], capacity=1)
    construction(g1, 0, set(g1.nodes()))
    # show_graph(A, 3)
    capacities=nx.get_edge_attributes(A,'capacity')
    if max(capacities.values())<=1:
        return g,None,1
    if max(capacities.values())==2:
        return g,None,2
    remove_edges = list()
    biconnected = list()
    for e in A.edges(data=True):
        if e[2]['capacity'] <= 2:
            remove_edges.append(e)
        if e[2]['capacity'] == 2:
            biconnected.append(e)
    A.remove_edges_from(remove_edges)
    for e in biconnected:
        cut_edges = list(nx.minimum_edge_cut(g1, e[0], e[1]))
        if cut_edges[1][0] in nx.node_connected_component(A, cut_edges[0][0]):
            g.add_edge(cut_edges[0][0], cut_edges[1][0])
            if cut_edges[1][1] in nx.node_connected_component(A, cut_edges[0][1]):
                g.add_edge(cut_edges[0][1], cut_edges[1][1])
            else:
                raise Exception("Biconnected edge does not match")
        elif cut_edges[1][0] in nx.node_connected_component(A, cut_edges[0][1]):
            g.add_edge(cut_edges[0][0], cut_edges[1][1])
            if cut_edges[1][1] in nx.node_connected_component(A, cut_edges[0][0]):
                g.add_edge(cut_edges[0][1], cut_edges[1][0])
            else:
                raise Exception("Biconnected edge does not match")
        else:
            raise Exception("Biconnected edge does not match")
        # print e[0],g1[e[0]]
        # print e[1],g1[e[1]]
    # show_graph(A, 1)
    # show_graph(g,1)
    G = list(nx.connected_components(A))
    connected_components=list(nx.connected_components(A))
    for i in connected_components:
        if len(i)<3:
            G.remove(i)
    if G==[]:
        return g,None,2
    return g,G,3


def read_edge(filename):
    with open(filename, 'r') as in_file:
        edge = in_file.read()
    edge = edge.replace('e ', '')
    edge = edge.split('\n')
    while edge[0][0] != 'p':
        edge.pop(0)
    attr = edge.pop(0)
    attr = attr.split()
    attr = int(attr.pop())
    while len(edge) > attr:
        edge.pop()
    int_edge = list()
    for e in edge:
        e = e.split()
        int_edge.append(map(int, e))
    if int_edge[len(int_edge)-1]==[]:
        int_edge.pop()
    return int_edge


def make_vars(g, d, width):
    nv = g.number_of_nodes()
    ne = g.number_of_edges()
    s = [[[0 for i1 in xrange(d)] for j in xrange(nv)] for k in xrange(nv)]
    l = [[0 for i1 in xrange(d)] for j in xrange(nv)]
    ad = [[[0 for i1 in xrange(d)] for j in xrange(ne)] for k in xrange(nv)]
    tor = [[[0 for i1 in xrange(d)] for j in xrange(nv)] for k in xrange(nv)]
    ctr_A = [[[[0 for i1 in xrange(width)] for j in xrange(d)] for k in xrange(ne)] for l1 in xrange(nv)]
    ctr_t = [[[[0 for i1 in xrange(width)] for j in xrange(d)] for k in xrange(nv)] for l1 in xrange(nv)]
    nvar = 1
    for k in xrange(d):
        for i in xrange(nv):
            l[i][k] = nvar
            nvar += 1
            for j in xrange(i, nv):
                s[i][j][k] = nvar
                nvar += 1
            for j in xrange(nv):
                tor[i][j][k] = nvar
                nvar += 1
                for l1 in xrange(width):
                    ctr_t[i][j][k][l1] = nvar
                    nvar += 1
            for j in xrange(ne):
                ad[i][j][k] = nvar
                nvar += 1
                for l1 in xrange(width):
                    ctr_A[i][j][k][l1] = nvar
                    nvar += 1
    return s, l, ad, tor, ctr_A, ctr_t, nvar-1


def generate_encoding(g, d, width):
    nv = g.number_of_nodes()
    ne = g.number_of_edges()
    s, l, ad, tor, ctr_a, ctr_t, n_var = make_vars(g, d, width)
    encoding = ''
    n_clauses = 0
    for u in xrange(nv):
        for v in xrange(u, nv):
            encoding += '%i 0\n' % s[u][v][d - 1]
            n_clauses += 1
    for u in xrange(nv):
        encoding += '-%i 0\n' % s[u][u][0]
        n_clauses += 1
    for u in xrange(nv):
        for v in xrange(u, nv):
            for i in xrange(d - 1):
                encoding += '-%i %i 0\n' % (s[u][v][i], s[u][v][i + 1])
                n_clauses += 1
    for u in xrange(nv):
        for v in xrange(u + 1, nv):
            for w in xrange(v + 1, nv):
                for i in xrange(d):
                    encoding += '-%i -%i %i 0\n' % (s[u][v][i], s[u][w][i], s[v][w][i])
                    n_clauses += 1
                    encoding += '-%i -%i %i 0\n' % (s[u][v][i], s[v][w][i], s[u][w][i])
                    n_clauses += 1
                    encoding += '-%i -%i %i 0\n' % (s[u][w][i], s[v][w][i], s[u][v][i])
                    n_clauses += 1
    for u in xrange(nv):
        for v in xrange(u + 1, nv):
            for i in xrange(1,d):
                encoding += '-%i %i 0\n' % (s[u][v][i], s[u][u][i])
                encoding += '-%i %i 0\n' % (s[u][v][i], s[v][v][i])
                n_clauses += 2
    for i in xrange(d):
        for u in xrange(nv):
            encoding += '-%i %i 0\n' % (l[u][i], s[u][u][i])
            encoding += '%i -%i ' % (l[u][i], s[u][u][i])
            for v in xrange(u):
                encoding += '%i ' % (s[v][u][i])
            encoding += '0\n'
            n_clauses += 1
            for v in xrange(u):
                encoding += '-%i -%i 0\n' % (l[u][i], s[v][u][i])
                n_clauses += 1
    # for i in xrange(1, d - 1):
    #     for u in xrange(nv):
    #         encoding += '%i -%i 0\n' % (l[u][i], l[u][i + 1])
    #         n_clauses += 1

    for e in xrange(ne):
        # print e, g.edges()[e]
        for u in xrange(nv):
            v = min(g.edges()[e][0], g.edges()[e][1])
            w = max(g.edges()[e][0], g.edges()[e][1])
            if u <= v:
                # print u,
                for i in xrange(1,d-1):
                    encoding += '-%i -%i %i %i 0\n' % (l[u][i], s[u][v][i], s[u][w][i], ad[u][e][i])
                    encoding += '-%i -%i %i %i 0\n' % (l[u][i], s[u][w][i], s[u][v][i], ad[u][e][i])
                    n_clauses += 2
            if v < u <= w:
                # print u,
                for i in xrange(1,d-1):
                    encoding += '-%i -%i %i 0\n' % (l[u][i], s[u][w][i], ad[u][e][i])
                    n_clauses += 1
        # print "\n","*"*10
    for u in xrange(nv):
        for v in xrange(u, nv):
            for i in xrange(1,d - 1):
                encoding += '-%i -%i -%i %i 0\n' % (l[u][i + 1], l[v][i], s[u][v][i + 1], tor[u][v][i + 1])
                n_clauses += 1
    for u in xrange(nv):
        for v in xrange(u, nv):
            for i in xrange(d - 1):
                encoding += '-%i -%i %i %i 0\n' % (l[u][i + 1], s[u][v][i + 1], s[v][v][i], tor[u][v][i + 1])
                n_clauses += 1
    # for i in xrange(d):
    #     for u in xrange(nv):
    #         for e in xrange(ne):
    #             encoding += '-%i %i 0\n' % (ad[u][e][i], ctr_a[u][e][i][0])
    #             n_clauses += 1
    # for i in xrange(d):
    #     for u in xrange(nv):
    #         for e in xrange(ne - 1):
    #             for j in xrange(width):
    #                 encoding += '-%i %i 0\n' % (ctr_a[u][e][i][j], ctr_a[u][e + 1][i][j])
    #                 n_clauses += 1
    # for i in xrange(d):
    #     for u in xrange(nv):
    #         for e in xrange(ne - 1):
    #             for j in xrange(width - 1):
    #                 encoding += '-%i -%i %i 0\n' % (ad[u][e + 1][i], ctr_a[u][e][i][j], ctr_a[u][e + 1][i][j + 1])
    #                 n_clauses += 1
    # for i in xrange(d):
    #     for u in xrange(nv):
    #         for e in xrange(ne - 1):
    #             encoding += '-%i -%i 0\n' % (ad[u][e + 1][i], ctr_a[u][e][i][width - 1])
    #             n_clauses += 1
    # for i in xrange(d-1):
    #     for u in xrange(nv):
    #         for v in xrange(u, nv):
    #             encoding += '-%i %i 0\n' % (tor[u][v][i], ctr_t[u][v][i][0])
    #             n_clauses += 1
    # for i in xrange(d-1):
    #     for u in xrange(nv):
    #         for v in xrange(u,nv - 1):
    #             for j in xrange(width-1):
    #                 encoding += '-%i %i 0\n' % (ctr_t[u][v][i][j], ctr_t[u][v - 1][i][j])
    #                 n_clauses += 1
    # for i in xrange(d-1):
    #     for u in xrange(nv):
    #         for v in xrange(u,nv - 1):
    #             for j in xrange(width - 2):
    #                 encoding += '-%i -%i %i 0\n' % (tor[u][v + 1][i], ctr_t[u][v][i][j], ctr_t[u][v + 1][i][j + 1])
    #                 n_clauses += 1
    # for i in xrange(d-1):
    #     for u in xrange(nv):
    #         for v in xrange(u,nv - 2):
    #             encoding += '-%i -%i 0\n' % (tor[u][v][i], ctr_t[u][v - 1][i][width - 1])
    #             n_clauses += 1
    # i=d-1
    # for u in xrange(nv):
    #     for v in xrange(u, nv):
    #         encoding += '-%i %i 0\n' % (tor[u][v][i], ctr_t[u][v][i][0])
    #         n_clauses += 1
    # for u in xrange(nv):
    #     for v in xrange(u,nv - 1):
    #         for j in xrange(width):
    #             encoding += '-%i %i 0\n' % (ctr_t[u][v][i][j], ctr_t[u][v - 1][i][j])
    #             n_clauses += 1
    # for u in xrange(nv):
    #     for v in xrange(u,nv - 1):
    #         for j in xrange(width - 1):
    #             encoding += '-%i -%i %i 0\n' % (tor[u][v + 1][i], ctr_t[u][v][i][j], ctr_t[u][v + 1][i][j + 1])
    #             n_clauses += 1
    # for u in xrange(nv):
    #     for v in xrange(u,nv - 1):
    #         encoding += '-%i -%i 0\n' % (tor[u][v][i], ctr_t[u][v - 1][i][width - 1])
    #         n_clauses += 1
    a_count = ''
    for i in xrange(1,d-1):
        for u in xrange(nv):
            # counter += 'card '
            for e in xrange(ne):
                a_count += '%i ' % ad[u][e][i]
            a_count += '<= %i\n' % width
    a_clauses=0
    a_count,n_var,a_clauses=counter(a_count,n_var)
    t_count=''
    for i in xrange(d - 1):
        for u in xrange(nv):
            # counter += 'card '
            for v in xrange(u, nv):
                t_count += '%i ' % tor[u][v][i]
            t_count += '<= %i\n' % (width - 1)
    t_clauses=0
    t_count,n_var,t_clauses=counter(t_count,n_var)
    clauses=0
    count=''
    for u in xrange(nv):
        # counter += 'card '
        for v in xrange(u, nv):
            count += '%i ' % tor[u][v][d - 1]
        count += '<= %i\n' % width
    count,n_var,clauses=counter(count,n_var)
    preamble = 'p cnf %i %i\n' % (n_var, n_clauses+a_clauses+t_clauses+clauses)
    final = preamble +t_count+a_count+ count + encoding
    # print nvar-1,n_clauses
    # preamble = 'p cnf %i %i\n' % (nvar - 1, n_clauses)
    # final = preamble + encoding
    if ' 0 ' in final or ' -0 ' in final:
        raise Exception("problem!")
    return final


def verify_decomp(g, decomp, width, length):
    sys.stderr.write("\nChecking decomposition\n")
    sys.stderr.flush()
    if not nx.is_tree(decomp):
        raise Exception("Decomposition is not tree")
    if set(decomp.nodes(data=True)[0][1]['component']) != set(g.nodes()):
        s = str(set(g.nodes()) - set(decomp.nodes(data=True)[0][1]['component']))
        raise Exception("All nodes are not covered" + s)
    for curr in decomp.nodes():
        tor = len(decomp[curr])
        vertices = set()
        children = len(decomp[curr])
        component = set(decomp.nodes(data=True)[curr][1]['component'])
        for node in decomp[curr]:
            # print decomp.nodes(data=True)[node][1]
            # exit(0)
            nodecomp = set(decomp.nodes(data=True)[node][1]['component'])
            if not nodecomp <= component:
                s = "Parent node missing vertices " + str(nodecomp) + str(component)
                raise Exception(s)
            vertices = vertices.union(nodecomp)
            for node1 in decomp[curr]:
                if node1 == node:
                    continue
                node1comp = set(decomp.nodes(data=True)[node1][1]['component'])
                if node1comp & nodecomp:
                    s = 'Appearing in multiple branches ' + str(node1comp) + str(nodecomp)
                    raise Exception(s)

        if decomp.nodes(data=True)[curr][1]['level'] < length-1:
            if len(set(component) - vertices) + len(decomp[curr]) + 1 > width:
                raise Exception("Higher torso size for bag %i than width %i" % (curr, width))
        elif len(set(component) - vertices) + len(decomp[curr]) > width:
            raise Exception("Higher torso size for bag %i than width %i" % (curr, width))
        incident_edges = list()
        for e in g.edges(list(vertices)):
            if not (e[0] in component and e[1] in component):
                incident_edges.append(e)
                # try:
                #     incident_edges.remove(e)
                #     sys.stderr.write("\ndeleted edge %s\n"%str(e))
                # except:
                #     try:
                #         incident_edges.remove((e[1],e[0]))
                #         sys.stderr.write("\ndeleted edge %s\n" % str(e[1],e[0]))
                #     except:
                #         sys.stderr.write("\nskipped edge %s\n"%str(e))
                #         pass
        if len(incident_edges) > width:
            raise Exception("Higher adhesion than width")
    sys.stderr.write("Correct decomposition\n")
    sys.stderr.flush()
    # for node in decomp.nodes(data=True):
    #     print node, decomp[node[0]]


def decode_output(sol, g, d, width):
    with open(sol, 'r') as out_file:
        out = out_file.read()
    out = out.split('\n')
    # print out
    out.pop()
    out = out[0]
    out = out.split(' ')
    out = map(int, out)
    out.pop()
    nv = g.number_of_nodes()
    ne = g.number_of_edges()
    # print g.edges()
    s, l, ad, tor, ctr_a, ctr_t, nvar = make_vars(g, d, width)
    components = list()
    adh = list()
    for i in xrange(d):
        temp = list()
        ad_temp = list()
        for u in xrange(nv):
            if out[l[u][i] - 1] > 0:
                curr_component = list()
                curr_ad = list()
                for v in xrange(u, nv):
                    if out[s[u][v][i] - 1] > 0:
                        curr_component.append(v)
                for e in xrange(ne):
                    if out[ad[u][e][i] - 1] > 0:
                        curr_ad.append((g.edges()[e]))
                temp.append(curr_component)
                ad_temp.append(curr_ad)
        if temp not in components:
            components.append(temp)
            adh.append(ad_temp)
    for i in components:
        sys.stderr.write(str(i)+'\n')
    decomp = nx.DiGraph()
    node = dict()
    level = dict()
    for i in xrange(len(components) - 1, -1, -1):
        for j1 in xrange(len(components[i])):
            j = components[i][j1]
            node_len = len(node)
            if set(j) in node.values():
                level[frozenset(j)] = i
                continue
            level[frozenset(j)] = i
            decomp.add_node(node_len, component=j, level=i)
            node[node_len] = frozenset(j)
            for n in decomp.nodes():
                if node[n] > set(j) and level[node[n]] == i + 1:
                    decomp.add_edge(n, node_len, label=adh[i][j1])
    # show_graph(decomp,1)
    # show_graph(g,1)
    # for n in decomp.nodes(data=True):
    #     sys.stderr.write(str(n)+'\n')
    # for e in decomp.edges(data=True):
    #     sys.stderr.write(str(e)+'\n')
    verify_decomp(g, decomp, width, len(components))


def parse_args():
    parser = argparse.ArgumentParser(description='%(prog)s -f instance')
    parser.add_argument('-f', '--file', dest='instance', action='store', type=lambda x: os.path.realpath(x),
                        default=None, help='instance')
    parser.add_argument('-d', '--depth', dest='d', action='store', type=int, default=-1, help='depth')
    parser.add_argument('-w', '--width', dest='width', action='store', type=int, default=-1, help='width')
    parser.add_argument('-o', '--timeout', dest='timeout', action='store', type=int, default=900, help='timeout for each SAT call')
    parser.add_argument('-t', '--temp', dest='temp', action='store', type=str, default='/home/neha/temp/',
                        help='temporary folder')
    parser.add_argument('-s', '--solver', dest='solver', action='store', type=str, default='glucose',
                        help='SAT solver')
    args = parser.parse_args()
    return args


def show_graph(graph, layout, nolabel=0):
    """ show graph
    layout 1:graphviz,
    2:circular,
    3:spring,
    4:spectral,
    5: random,
    6: shell
    """

    global pos
    m = graph.copy()
    if layout == 1:
        pos = graphviz_layout(m)
    elif layout == 2:
        pos = nx.circular_layout(m)
    elif layout == 3:
        pos = nx.spring_layout(m)
    elif layout == 4:
        pos = nx.spectral_layout(m)
    elif layout == 5:
        pos = nx.random_layout(m)
    elif layout == 6:
        pos = nx.shell_layout(m)
    if not nolabel:
        nx.draw_networkx_edge_labels(m, pos)
        node_labels = nx.get_node_attributes(m, 'component')
        nx.draw_networkx_labels(m, pos, labels=node_labels)
    nx.draw_networkx_nodes(m, pos)
    nx.draw_networkx_labels(m, pos)
    # write_dot(m, "m1.dot")
    # os.system("dot -Tps m1.dot -o m1.ps")
    nx.draw(m, pos)
    plt.show()


def main():
    args = parse_args()
    instance = args.instance
    d = args.d
    width = args.width
    temp = args.temp
    timeout=args.timeout
    cpu_time=time.time()
    solver=args.solver
    if instance is not None:
        edge = read_edge(instance)
        g = nx.MultiGraph()
        g.add_edges_from(edge)
        instance = os.path.basename(instance)
        instance = instance.split('.')
        instance = instance[0]
    else:
        g = nx.complete_bipartite_graph(4, 4)
        # g = nx.gnp_random_graph(30, 0.1)
        # g = nx.complete_graph(6)
        # g=nx.barbell_graph(5,1)
        # g=nx.cycle_graph(4)
        # g.add_edge(0,2)
        # g=nx.grid_2d_graph(4,4)
        # g=nx.configuration_model([3,3],)
        # g1=nx.complete_graph(6)
        # g=nx.complete_graph(7)
        # g=nx.disjoint_union(g1,g)
        # g.add_edge(g1.number_of_nodes()-1,g1.number_of_nodes())
        # g.add_edge(g1.number_of_nodes()-2,g1.number_of_nodes()+1)
        # g1=nx.complete_graph(5)
        # g=nx.disjoint_union(g1,g)
        # g.add_edge(g1.number_of_nodes()-1,g1.number_of_nodes())
        instance = 'random'
        # show_graph(g, 1)
    g=nx.convert_node_labels_to_integers(g,first_label=0)
    max_degree=max(g.degree().values())
    n=g.number_of_nodes()
    m=g.number_of_edges()
    lb=0
    ub=0
    to=False
    prep_time=time.time()
    g,G,biconn = reduce_graph(g)
    prep_time=time.time()-prep_time
    print 'treecut2sat', instance, n, m, max_degree,
    if biconn<=2:
        print biconn,lb,ub,to, time.time()-cpu_time,prep_time
        exit(0)
    comp_ver=list()
    solving_time=list()
    encoding_time=list()
    for ver in G:
        comp_ver.append(len(ver))
        if len(ver) == 1:
            continue
        if d == -1:
            d = len(ver)
        if width == -1:
            u = 2
            width=len(ver)
        else:
            u = width
        for i in xrange(len(ver)-6, u-1,-1):
            g_cur = g.subgraph(ver)
            g_cur = nx.convert_node_labels_to_integers(g_cur)
            # print g_cur.edges()
            encode_time=time.time()
            encoding = generate_encoding(g_cur, d, i)
            cnf = temp + instance + '_' + str(i) + ".cnf"
            with open(cnf, 'w') as ofile:
                ofile.write(encoding)
            encode_time=time.time()-encode_time
            encoding_time.append(encode_time)
            sol = temp + instance + '_' + str(i) + '.sol'
            solving=time.time()
            cmd = [solver, '-cpu-lim=%i'%timeout, cnf, sol]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, err = p.communicate()
            solving_time.append(time.time()-solving)
            rc = p.returncode
            sys.stderr.write('*' * 10+'\n')
            # print err
            sys.stderr.write("%i %i\n"%(i, rc))
            if rc==0:
                to=True
                if lb<i:
                    ub=i
            if rc == 20:
                if to:
                    if lb>i:
                        lb=i-1
                # decode_output(sol, g_cur, d, i)
                if width > i:
                    width = i
                d = -1
                break
    cpu_time=time.time()-cpu_time
    sys.stderr.write("*" * 10+'\n')
    print width,lb,ub,to,cpu_time,prep_time,sum(comp_ver),len(comp_ver),
    print sum(encoding_time),
    print sum(solving_time),
    for i in solving_time:
        print i,
    for i in comp_ver:
        print i,


if __name__ == "__main__":
    main()
