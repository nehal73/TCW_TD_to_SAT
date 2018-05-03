# coding=utf-8
import time
import argparse
import os, sys
import networkx as nx
import numpy as np
import subprocess
from networkx.drawing.nx_agraph import *


import matplotlib.pyplot as plt


def counter(s, n_var):
    lines = s.split('\n')
    n_clauses = 0
    encoding = ''
    for s in lines:
        # print s
        s = s.split()
        if s == []:
            continue
        count = int(s.pop())
        operator = s.pop()
        if int(count)==0:
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
    return encoding, n_var, n_clauses


def apex_vertices(g):
    buff = 0
    delete_vertices = list()
    for u, degree in g.degree().iteritems():
        if degree == g.number_of_nodes() - 1:
            delete_vertices.append(u)
    g.remove_nodes_from(delete_vertices)
    buff += len(delete_vertices)
    nx.convert_node_labels_to_integers(g, first_label=0)
    return g, buff


def degree_one_reduction(g):
    """
    Removes all but one degree one neighbours of one vertex 
    :returns g: reduced graph
    :type g: networkx graph
    """
    nodes = set()
    for u in g.nodes():
        deg = 0
        for v in g.neighbors(u):
            if g.degree(v) == 1:
                if deg == 0:
                    deg = 1
                else:
                    nodes = nodes.union({v})
    g.remove_nodes_from(list(nodes))
    g = nx.convert_node_labels_to_integers(g, first_label=0)
    return g


def read_edge(filename):
    """
    read dimacs graph file
    :returns int_edge: a list of edges
    :type filename: filename of the graph file
    """
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
    if int_edge[len(int_edge) - 1] == []:
        int_edge.pop()
    return int_edge


def make_vars(g, width):
    """
    Makes variables for treedepth
    :type g: networkx graph (input graph)
    :type width: int (width to encode)
    """
    nv = g.number_of_nodes()
    # p = range(1, nv * nv + 1)
    # nvar = nv * nv + 1
    # p = np.reshape(p, (nv, nv))
    nvar = 1
    p = [[0 for u in xrange(nv)] for v in xrange(nv)]
    for u in xrange(nv):
        for v in xrange(nv):
            if u != v:
                p[u][v] = nvar
                nvar += 1
    a = range(nvar, nvar + nv * nv)
    nvar += nv * nv
    a = np.reshape(a, (nv, nv))
    r = range(nvar, nvar + nv)
    nvar += nv
    return (p, a, r, nvar)


def generate_encoding(g, width):
    """
    :returns
    :type g: networkx graph
    :param g: input graph
    :type width: int (width)
    """
    (p, a, r, nvar) = make_vars(g=g, width=width)
    s = ''
    nclauses = 0
    nv = g.number_of_nodes()
    g = nx.convert_node_labels_to_integers(g, first_label=0)
    # print g.nodes()
    for e in g.edges():
        s += '%i %i 0\n' % (a[e[0]][e[1]], a[e[1]][e[0]])
        nclauses += 1
    for u in xrange(nv):
        for v in xrange(nv):
            if u == v:
                # s+='-%i 0\n'%p[u][u]
                # nclauses+=1
                s += '-%i 0\n' % (a[u][u])
                nclauses += 1
                continue
            s += '-%i -%i 0\n' % (a[u][v], a[v][u])
            s += '-%i %i 0\n' % (p[v][u], a[v][u])
            nclauses += 1
            for w in xrange(nv):
                if u == w or v == w:
                    continue
                s += '-%i -%i %i 0\n' % (p[w][u], a[v][w], a[v][u])
                nclauses += 1
                s += '-%i -%i %i 0\n' % (p[w][u], a[v][u], a[v][w])
                nclauses += 1
            for w in xrange(v + 1, nv):
                if w == u:
                    continue
                s += '-%i -%i 0\n' % (p[v][u], p[w][u])
                nclauses += 1
    for u in xrange(nv):
        for v in xrange(nv):
            if u == v:
                continue
            s += '-%i -%i 0\n' % (r[u], a[v][u])
            nclauses += 1
    for u in xrange(nv):
        s += '%i ' % r[u]
    s += '0\n'
    for u in xrange(nv):
        for v in xrange(u+1,nv):
            s+='-%i -%i 0\n'%(r[u],r[v])
            nclauses+=1
    # s+='%i 0\n'%r[1]
    for u in xrange(nv):
        s += '%i ' % r[u]
        for v in xrange(nv):
            if u == v:
                continue
            s += '%i ' % p[v][u]
        s += '0\n'
        nclauses += 1
    for u, v in g.edges():
        if set(g.neighbors(u)) - {v} <= set(g.neighbors(v)) - {u}:
            s += '%i 0\n' % a[v][u]
            nclauses += 1
        elif set(g.neighbors(v)) - {u} < set(g.neighbors(u)) - {v}:
            s += '%i 0\n' % a[u][v]
            nclauses += 1
    # s+='%i 0\n'%p[1][2]
    # nclauses+=1
    # s+='%i 0\n'%p[4][1]
    # nclauses+=1
    count = ''
    for u in xrange(nv):
        for v in xrange(nv):
            count += '%i ' % a[v][u]
        count += ' <= %i\n' % width
        nclauses += 1
    clauses = 0
    count, nvar, clauses = counter(count, nvar)
    if ' 0 ' in s or '-0' in s:
        raise Exception("0 in formula")
    preamble = 'p cnf %i %i\n' % (nvar, nclauses + clauses)
    return preamble + count + s


def verify_decomp(g, s, anc, width, root):
    sys.stderr.write("Validating tree depth decomposition\n")
    sys.stderr.flush()
    # print g.edges()
    for i in anc:
        if len(i) > width:
            raise ValueError("more ancestors than width\n")
    for e in g.edges():
        try:
            nx.shortest_path(s, e[0], e[1])
        except:
            try:
                nx.shortest_path(s, e[1], e[0])
            except:
                raise Exception("Edge %i %i not covered\n" % (e[0], e[1]))
    for v, d in g.degree().iteritems():
        count = 0
        if d != 1:
            continue
        for i in root:
            try:
                if len(nx.shortest_path(s, i, v)) - 1 > width:
                    raise ValueError("depth of tree more than width\n")
            except:
                count += 1
                if count == len(root):
                    raise Exception("No root found\n")
                continue
    sys.stderr.write("Valid treedepth decomp\n")
    sys.stderr.flush()


def decode_output(sol, g, width):
    with open(sol, 'r') as out_file:
        out = out_file.read()
    out = out.split('\n')
    out.pop()
    out = out.pop()
    out = out.split(' ')
    out = map(int, out)
    out.pop()
    # print out
    nv = g.number_of_nodes()
    ne = g.number_of_edges()
    (p, a, r, nvar) = make_vars(g=g, width=width)
    s = nx.DiGraph()
    for u in xrange(nv):
        s.add_node(u)
    anc = list()
    root = list()
    for u in xrange(nv):
        if out[r[u] - 1] > 0:
            root.append(u)
        ancu = list()
        for v in xrange(nv):
            if u == v:
                continue
            if out[a[v][u] - 1] > 0:
                ancu.append(v)
            if out[p[v][u] - 1] > 0:
                s.add_edge(v, u)
        anc.append(ancu)
    show_graph(s, 1)
    verify_decomp(g=g, s=s, anc=anc, width=width, root=root)


def show_graph(graph, layout, nolabel=0):
    """ show graph
    layout 1:graphviz,
    2:circular,
    3:spring,
    4:spectral,
    5: random,
    6: shell
    """

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
    nx.draw_networkx_labels(m, pos)
    nx.draw_networkx_nodes(m, pos)
    # write_dot(m, "m1.dot")
    # os.system("dot -Tps m1.dot -o m1.ps")
    nx.draw(m, pos)
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description='%(prog)s -f instance')
    parser.add_argument('-f', '--file', dest='instance', action='store', type=lambda x: os.path.realpath(x),
                        default=None, help='instance')
    parser.add_argument('-o', '--timeout', dest='timeout', action='store', type=int, default=900,
                        help='timeout for each SAT call')
    parser.add_argument('-d', '--depth', dest='d', action='store', type=int, default=-1, help='depth')
    parser.add_argument('-w', '--width', dest='width', action='store', type=int, default=-1, help='width')
    parser.add_argument('-t', '--temp', dest='temp', action='store', type=str, default='/home/neha/temp/',
                        help='temporary folder')
    parser.add_argument('-s', '--solver', dest='solver', action='store', type=str, default='glucose',
                        help='SAT solver')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    instance = args.instance
    d = args.d  # 'glucose'
    # 'minicard_encodings_static'
    width = args.width
    temp = args.temp
    solver = args.solver
    # solver = 'minicard_encodings_static'
    timeout = args.timeout
    cpu_time = time.time()
    if instance != None:
        edge = read_edge(filename=instance)
        g = nx.MultiGraph()
        g.add_edges_from(ebunch=edge)
        instance = os.path.basename(instance)
        instance = instance.split('.')
        instance = instance[0]
    else:
        # g = nx.complete_bipartite_graph(6,9)
        # g=nx.complete_graph(7)
        # g=nx.cycle_graph(7)
        # g.add_edge(6,7)
        # g.add_edge(0,7)
        # g.add_edge(1,7)
        # g=nx.balanced_tree(2,3)
        instance = 'random'
        g = nx.path_graph(7)
        # show_graph(g,1)
    n = g.number_of_nodes()
    m = g.number_of_edges()
    prep_time = time.time()
    g=degree_one_reduction(g=g)
    buff = 0
    lb = 0
    ub = 0
    to = False
    g,buff=apex_vertices(g=g)
    prep_time = time.time() - prep_time
    print 'treedepthtcb2sat', instance, n, m,g.number_of_nodes()
    if g.number_of_nodes() <= 1:
        print buff,lb,ub,to, time.time() - cpu_time, prep_time
        exit(0)
    encoding_time = list()
    solving_time = list()
    g = nx.convert_node_labels_to_integers(g, first_label=0)
    if width == -1:
        for i in xrange(g.number_of_nodes() + 3,1,-1):
            encode_time = time.time()
            encoding = generate_encoding(g=g, width=i - 1)
            cnf = temp + instance + '_' + str(i) + ".cnf"
            with open(cnf, 'w') as ofile:
                ofile.write(encoding)
            # with open(cnf,'r') as ifile:
            #     s=ifile.read()
            #     print s
            encode_time = time.time() - encode_time
            encoding_time.append(encode_time)
            sol = temp + instance + '_' + str(i) + '.sol'
            cmd = [solver, '-cpu-lim=%i' % timeout, cnf, sol]
            # print cmd
            solving = time.time()
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, err = p.communicate()
            rc = p.returncode
            solving = time.time() - solving
            solving_time.append(solving)

            sys.stderr.write('*' * 10 + '\n')
            # print err
            sys.stderr.write("%i %i\n" % (i, rc))
            if rc == 0:
                to = True
                if lb == 0:
                    lb = i - 1
            if rc == 10:
                if to:
                    ub = i
                decode_output(sol=sol, g=g, width=i)
                print buff + i, buff + lb, buff + ub, to, time.time() - cpu_time, prep_time, sum(
                    encoding_time), sum(solving_time),
                for j in solving_time:
                    print j,
                exit(0)


if __name__ == "__main__":
    main()

#
# def generate_encoding(g,width):
#     p,a,nvar=make_vars(g,width)
#     s=''
#     nclauses=0
#     nv=g.number_of_nodes()
#     for e in g.edges():
#         s+='%i %i 0\n'%(p[e[0]][e[1]],p[e[1]][e[0]])
#         nclauses+=1
#     for u in xrange(nv):
#         for v in xrange(nv):
#             if u==v:
#                 continue
#             s+='-%i -%i 0\n'%(p[u][v],p[v][u])
#             s+='-%i %i 0\n'%(p[u][v],a[u][v][0])
#             for i in xrange(width):
#                 s+='-%i %i 0\n'%(a[u][v][i],p[u][v])
#                 nclauses += 1
#             for w in xrange(nv):
#                 if u==w or v==w:
#                     continue
#                 for i in xrange(width-1):
#                     s+='-%i -%i %i 0\n'%(p[u][w],a[w][v][i],a[u][v][i+1])
#                     nclauses += 1
#                 s+='-%i -%i -%i 0\n'%(p[u][v],p[w][v],p[u][w])
#                 nclauses += 1
#                 s += '-%i -%i -%i 0\n' % (p[u][v], p[w][v], p[w][u])
#                 nclauses += 1
#                 s+='-%i -%i %i 0\n'%(p[u][v],p[w][u],p[w][v])
#                 nclauses+=1
#                 for i in xrange(width):
#                     s+='-%i -%i -%i 0\n'%(a[u][v][i],a[w][v][i],p[u][w])
#                     nclauses += 1
#                     s+='-%i -%i -%i 0\n'%(a[u][v][i],a[w][v][i],p[w][u])
#                     nclauses += 1
#                 s+='-%i -%i 0\n'%(p[u][w],a[w][v][width-1])
#                 nclauses += 1
#     preamble='p cnf %i %i\n'%(nvar,nclauses)
#     return preamble+s
