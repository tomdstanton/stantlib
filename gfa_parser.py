#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    GFA Parser in Python using modified code from here:
    https://github.com/songbowang125/SVision/blob/d4daaa6ab5a4c374b2edce48d3a4be0d77d083e5/src/collection/graph.py
"""

from textwrap import wrap

__version__ = '0.0.1'
__author__ = "Tom Stanton"
__maintainer__ = "Tom Stanton"
__email__ = "tomdstanton@gmail.com"
__status__ = "Development"


# ............... GFA Parsing functions and classes ...............
class Node:
    """
    """

    def __init__(self, chr, ref_start, ref_end, read_start, read_end, seq, is_reverse, id, host):

        self.chr = chr
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_start = read_start
        self.read_end = read_end
        self.is_reverse = is_reverse
        self.id = id
        self.seq = seq
        self.host = host
        self.depth = 0
        self.node_is_dup = False
        self.dup_from = -1

    def add_depth(self, num):
        self.depth += num

    def set_dup_node(self, dup_from, dup_from_cord):
        if dup_from != -1:
            self.node_is_dup = True
            self.dup_from = dup_from
            self.dup_from_cord = dup_from_cord

    def to_string(self):
        if self.node_is_dup:
            return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}'.format(self.id, self.chr, self.ref_start, self.ref_end,
                                                                   self.read_start, self.read_end,
                                                                   '-' if self.is_reverse else '+', self.dup_from)
        else:
            return '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(self.id, self.chr, self.ref_start, self.ref_end,
                                                              self.read_start, self.read_end,
                                                              '-' if self.is_reverse else '+')


class Edge:
    """
    Used code from here:
    https://github.com/songbowang125/SVision/blob/d4daaa6ab5a4c374b2edce48d3a4be0d77d083e5/src/collection/graph.py
    """

    def __init__(self, node1, node1_is_reverse, node2, node2_is_reverse, id):
        self.node1 = node1
        self.node1_is_reverse = node1_is_reverse
        self.node2 = node2
        self.node2_is_reverse = node2_is_reverse
        self.id = id
        self.edge_is_dup = False

    def set_dup_edge(self):
        self.edge_is_dup = True

    def to_string(self):
        return '{0}\t{1}\t{2}\t{3}\t{4}'.format(self.id, self.node1, '-' if self.node1_is_reverse else '+', self.node2,
                                                '-' if self.node2_is_reverse else '+')


class Graph:
    """
    """

    def __init__(self, nodes, edges, qname=""):
        self.nodes = nodes
        self.edges = edges
        self.appear_time = 1
        self.qname = qname


def parse_gfa_file(gfa_path):
    """
    """
    nodes = []
    edges = []
    with open(gfa_path) as fin:
        for line in fin.readlines():
            line_split = line.strip().split('\t')
            if line_split[0] == 'S':
                node_id = line_split[1]
                node_seq = line_split[2]
                node_host = line_split[3].split(':')[-1]
                node_start = line_split[4].split(':')[-1]

                if len(line_split) == 8:
                    node_is_dup = True
                    dup_from = line_split[7].split(':')[2]
                    dup_from_cord = int(line_split[7].split(':')[3])

                else:
                    node_is_dup = False
                    dup_from = -1
                    dup_from_cord = -1

                node = Node(-1, node_start, -1, node_start, -1, node_seq, False, node_id, node_host)
                node.set_dup_node(dup_from, dup_from_cord)
                nodes.append(node)
            elif line_split[0] == 'L':
                edge_node1 = line_split[1]

                edge_node1_is_reverse = True if line_split[2] == '-' else False
                edge_node2 = line_split[3]
                edge_node2_is_reverse = True if line_split[4] == '-' else False
                edges.append(Edge(edge_node1, edge_node1_is_reverse, edge_node2, edge_node2_is_reverse, 0))
            else:
                pass
    return Graph(nodes, edges)


def fasta_wrap(fasta_string):
    return '\n'.join(wrap(str(fasta_string), 60))


