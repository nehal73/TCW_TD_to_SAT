# coding=utf-8
# coding=utf-8
import argparse
import os
import sys
import time
import networkx as nx
import numpy as np
import subprocess

from networkx.drawing.nx_agraph import *
import networkx as nx
# import matplotlib.pyplot as plt
import signal


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


def make_vars(g,width):
    nv=g.number_of_nodes()
    p=[[[0 for i in xrange(width)]for j in xrange(nv)] for k in xrange(nv)]
    nvar=1
    for u in xrange(nv):
        for v in xrange(u,nv):
            for i in xrange(width):
                p[u][v][i]=nvar
                nvar+=1
    return p,nvar-1


def generate_encoding(g,width):
    s,nvar=make_vars(g,width)
    encoding=''
    nclauses=0
    nv=g.number_of_nodes()
    for u in xrange(nv):
        for v in xrange(u,nv):
            encoding += '%i 0\n' % (s[u][v][width-1])
            nclauses += 1
            encoding += '-%i 0\n' % (s[u][v][0])
            nclauses += 1
    for u in xrange(nv):
        for v in xrange(u, nv):
            for i in xrange(1,width):
                encoding += '-%i %i 0\n' % (s[u][v][i-1], s[u][v][i])
                nclauses += 1
    for u in xrange(nv):
        for v in xrange(u + 1, nv):
            for w in xrange(v + 1, nv):
                for i in xrange(width):
                    encoding += '-%i -%i %i 0\n' % (s[u][v][i], s[u][w][i], s[v][w][i])
                    nclauses += 1
                    encoding += '-%i -%i %i 0\n' % (s[u][v][i], s[v][w][i], s[u][w][i])
                    nclauses += 1
                    encoding += '-%i -%i %i 0\n' % (s[u][w][i], s[v][w][i], s[u][v][i])
                    nclauses += 1
    for u in xrange(nv):
        for v in xrange(u + 1, nv):
            for i in xrange(width):
                    encoding += '-%i %i 0\n' % (s[u][v][i], s[u][u][i])
                    encoding += '-%i %i 0\n' % (s[u][v][i], s[v][v][i])
                    nclauses += 2
    for u in xrange(nv):
        for v in xrange(u+1,nv):
            for i in xrange(1,width):
                encoding+='-%i %i %i 0\n'%(s[u][v][i],s[u][u][i-1],s[v][v][i-1])
                nclauses+=1
    for e in g.edges():
        u=min(e)
        v=max(e)
        for i in xrange(1,width):
            encoding+='-%i %i -%i %i 0\n'%(s[u][u][i],s[u][u][i-1],s[v][v][i],s[u][v][i])
            nclauses+=1
            encoding+='-%i %i -%i %i 0\n'%(s[u][u][i],s[v][v][i-1],s[v][v][i],s[u][v][i])
            nclauses+=1
    preamble='p cnf %i %i\n'%(nvar,nclauses)
    return preamble+encoding


def decode_output(sol,g,width):
    with open(sol, 'r') as out_file:
        out = out_file.read()
    out = out.split('\n')
    out = out[0]
    out = out.split(' ')
    out = map(int, out)
    out.pop()
    nv = g.number_of_nodes()
    ne = g.number_of_edges()
    s,nvar = make_vars(g, width)
    components=list()
    for i in xrange(width-1,0,-1):
        level=list()
        for u in xrange(nv):
            ver=list()
            for v in xrange(u,nv):
                if out[s[u][v][i]-1]>0:
                    ver.append(v)
            do_not_add=0
            for v in level:
                if set(ver).issubset(set(v)):
                    do_not_add=1
            if do_not_add==0:
                level.append(ver)
        components.append(level)
    for i in components:
        sys.stderr.write(str(i)+'\n')
    sys.stderr.write('\n'+"*"*10+'\n')
    decomp=nx.DiGraph()
    root=list()
    level_i=list()
    for i in xrange(width-1,0,-1):
        level = list()
        # sys.stderr.write('\n'+'*'*10+'\n')
        for u in xrange(nv):
            if out[s[u][u][i]-1]>0 and out[s[u][u][i-1]-1]<0:
                edge_add=False
                if i==width-1:
                    root.append(u)
                decomp.add_node(u,level=i)
                # sys.stderr.write("%i "%u)
                level.append(u)
                if level_i!=[]:
                    for v in level_i[len(level_i)-1]:
                        if out[s[min(u,v)][max(u,v)][i+1]-1]>0:
                            decomp.add_edge(v,u)
                            edge_add=True
                if not edge_add:
                    if level_i!=[]:
                        level_u=g.number_of_nodes()
                        v_u=-1
                        for v in decomp.nodes(data=True):
                            if v[0]==u:
                                continue
                            if out[s[min(u,v[0])][max(u,v[0])][v[1]['level']]-1]>0:
                                if level_u>v[1]['level']:
                                    level_u=v[1]['level']
                                    v_u=v[0]
                        decomp.add_edge(v_u,u)
        level_i.append(level)
        # print level
    # show_graph(decomp,1)
    verify_decomp(g=g,s=decomp,width=width,root=root)


def verify_decomp(g, s, width, root):
    sys.stderr.write("\nValidating tree depth decomposition\n")
    sys.stderr.flush()
    # print g.edges()
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
                    raise Exception("No root found for %i\n"%v)
                continue
    sys.stderr.write("Valid treedepth decomp\n")
    sys.stderr.flush()


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
    #plt.show()


class Timer(object):
    def __init__(self, time_list=None):
        self.time_list = time_list

    def __enter__(self):
        self.start = time.time()
        self.end = self.duration = None
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()
        self.duration = self.end - self.start
        if self.time_list is not None:
            self.time_list.append(self.duration)
        if exc_val is None:
            self.err = None
        else:
            self.err = (exc_type, exc_val, exc_tb)
            print >> sys.stderr, "\ntimed block terminated abruptly after", self.duration, "seconds"
            print >> sys.stderr, self.err


def solve_component(g, cli_args):
    lb = 0
    ub = 0
    to = False
    encoding_times = list()
    solving_times = list()
    n = g.number_of_nodes()
    if n <= 1:
        return n, n, n, to, encoding_times, solving_times
    temp = os.path.abspath(cli_args.temp)
    instance = cli_args.instance
    for i in xrange(g.number_of_nodes() + 2, 1, -1):
        with Timer(time_list=encoding_times):
            encoding = generate_encoding(g, i)
            cnf = os.path.join(temp, instance + '_' + str(i) + ".cnf")
            with open(cnf, 'w') as ofile:
                ofile.write(encoding)
        sol = os.path.join(temp, instance + '_' + str(i) + ".sol")
        cmd = [cli_args.solver, '-cpu-lim={}'.format(cli_args.timeout), cnf, sol]
        # print cmd
        with Timer(time_list=solving_times):
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, err = p.communicate()
            rc = p.returncode
        sys.stderr.write('*' * 10 + '\n')
        # print output, err
        sys.stderr.write("%i %i\n" % (i - 1, rc))
        if rc == 0:
            to = True
            if lb == ub == 0:  # first timeout, record ub
                ub = i
        if rc == 20:
            if to:
                #lb = i - 2
                lb = i
            if lb == ub == 0:  # never timed out
                    lb = ub = i
            return i, lb, ub, to, encoding_times, solving_times


def signal_handler(signum, frame):
    print "aborting due to signal", signum
    print "* final treedepth ?"
    sys.exit(0)


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
                        help='SAT solver')  # or 'minicard_encodings_static'
    parser.add_argument('-n', '--no-preprocess', dest="preprocess", action='store_false', help="Turn off preprocessing")
    args = parser.parse_args()
    return args


def main():
    cpu_time = time.time()
    args = parse_args()
    instance = args.instance
    if instance != None:
        edge = read_edge(instance)
        g = nx.MultiGraph()
        g.add_edges_from(edge)
        instance = os.path.basename(instance)
        instance = instance.split('.')
        instance = instance[0]
    else:
        # g = nx.complete_bipartite_graph(2,2)
        # g=nx.complete_graph(7)
        g=nx.balanced_tree(2,2)
        # g=nx.cycle_graph(15)
        # g=nx.path_graph(70)
        instance = 'random'
        # show_graph(g,6)
    args.instance = instance
    n = g.number_of_nodes()
    m = g.number_of_edges()
    buff = 0
    with Timer() as prep_timer:
        if args.preprocess:
            print >> sys.stderr, "preprocessing..."
            g=degree_one_reduction(g=g)
            g,buff=apex_vertices(g=g)
    print "* buffer verts:", buff
    print 'treedepthp2sat', instance, n, m,g.number_of_nodes(),buff
    if args.width != -1:
        return
    ncomps = nx.number_connected_components(g)
    icomp = 1
    global_lb, global_ub = 1e9, -1
    if ncomps == 0:  # only empty graph remains after preprocessing, for loop won't be triggered
        global_lb = global_ub = 0
    for subgraph in nx.connected_component_subgraphs(g):
        print >> sys.stderr, "\ncomponent", icomp, "of", ncomps
        component = nx.convert_node_labels_to_integers(subgraph,first_label=0)
        i, lb, ub, to, encoding_time, solving_time = solve_component(component, args)
        # decode_output(sol=os.path.join(args.temp, instance + '_' + str(i+1) + ".sol"), g=component, width=i+1)
        print i-2, lb, ub, to, time.time() - cpu_time, prep_timer.duration, \
            sum(encoding_time), sum(solving_time),
        for j in solving_time:
            print j,
        print
        print >> sys.stderr, "* component treedepth range: [{}-{}]".format(lb, ub)
        global_lb = min(global_lb, lb)
        global_ub = max(global_ub, ub)
        icomp += 1
    print "\n* final treedepth:",
    if global_ub == global_lb:
        print buff + global_ub,
    else:
        print "[{}-{}]".format(buff + global_lb, buff + global_ub),
    print "\ttotal-time: {:.2f}s".format(time.time()-cpu_time)


if __name__ == "__main__":
    signal.signal(signal.SIGHUP, signal_handler)
    main()
