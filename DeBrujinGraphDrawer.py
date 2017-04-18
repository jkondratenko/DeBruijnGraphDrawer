#!/usr/bin/env python3

import argparse
from tqdm import tqdm

# A small program which parses provided fastq or fasta file with genome
# reads, builds DeBrujin Graph, compress non-branching edges and outputs
# dot.file. Dot file can be used for visualization.


class Edge:
    show_sequences = False

    def __init__(self, v1, v2, seq="", coverage=0):
        # Edge instance should contains
        # starting and finishing vertices,
        # coverage and edge sequence
        self.v1 = v1
        self.v2 = v2
        if not coverage:
            self.coverage = 1
        else:
            self.coverage = coverage
        if not seq:
            self.sequence = v1.sequence + v2.sequence[-1]
        else:
            self.sequence = seq

    def inc_coverage(self):
        self.coverage += 1

    def __len__(self):
        return len(self.sequence)

    def merge(self, following_edge):
        raise NotImplementedError

    def __str__(self):
        return "sequence: %s, coverage: %s" % (str(self.sequence), str(self.coverage))


class Vertex:
    show_sequences = False

    def __init__(self, seq):
        # Vertex instance contains
        # sequence, output and input edges
        self.sequence = seq
        self.inn = set()
        self.out = set()

    def add_edge(self, other):
        # Increases coverage if the edge already exists
        assert (self.sequence[1:] == other.sequence[:-1])
        for out_edge in self.out:
            if out_edge.v2 == other:
                out_edge.inc_coverage()
                break
        else:
            new_edge = Edge(self, other)
            self.out.add(new_edge)
            other.inn.add(new_edge)

    def __str__(self):
        return "sequence: %s, in_edges: %s, out_edges: %s" % (
            self.sequence, str([str(x) for x in self.inn]), str([str(y) for y in self.out]))

    def can_compress(self):
        return len(self.inn) == 1 and len(self.out) == 1

    def compress(self):
        # Returns False, if cannot be compressed
        # Otherwise compresses this vertex and returns true
        in_edge, out_edge = self.inn.pop(), self.out.pop()
        from_v, to_v = in_edge.v1, out_edge.v2
        for out_edge in from_v.out:
            if out_edge.v2 == self:
                from_v.out.discard(out_edge)
                break
        for in_edge in to_v.inn:
            if in_edge.v1 == self:
                to_v.inn.discard(in_edge)
                break
        new_seq = in_edge.sequence + out_edge.sequence[len(self.sequence):]
        new_cov = (in_edge.coverage * len(in_edge.sequence) + out_edge.coverage * len(out_edge.sequence)) / len(
                     in_edge.sequence) + len(out_edge.sequence)
        compressed_edge = Edge(from_v, to_v, new_seq, new_cov)
        from_v.out.add(compressed_edge)
        to_v.inn.add(compressed_edge)


class Graph:
    k = None

    def __init__(self):
        # Contains all vertices
        self.d = {}

    def add_edge(self, seq1, seq2):
        # Increases coverage if the edge already exists
        from_v = None
        to_v = None
        if seq1 in self.d:
            from_v = self.d[seq1]
        if seq2 in self.d:
            to_v = self.d[seq2]
        if not from_v:
            from_v = Vertex(seq1)
            self.d[seq1] = from_v
        if not to_v:
            to_v = Vertex(seq2)
            self.d[seq2] = to_v
        from_v.add_edge(to_v)

    def add_seq(self, seq):
        # Adds edges between all k-mers in the sequence
        for i in range(len(seq) - Graph.k + 1):
            k_plus = seq[i:i + Graph.k]
            first, second = k_plus[:-1], k_plus[1:]
            self.add_edge(first, second)

    def compress(self):
        while True:
            for seq, v in self.d.items():
                if v.can_compress():
                    v.compress()
                    del self.d[seq]
                    break
            else:
                break

    def save_dot(self, outp):
        printed_edges = set()
        with open(outp, "w") as outp:
            outp.write("digraph {\n")
            for v in self.d.values():
                for out_edge in v.out:
                    if out_edge not in printed_edges:
                        outp.write("%s -> %s\n" % (out_edge.v1.sequence, out_edge.v2.sequence))
                        printed_edges.add(out_edge)
            outp.write("}")


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def reverse_complement(seq):
    return ''.join(complement[nt] for nt in seq[::-1])


def read_fastq(f):
    for line in f:
        name = line.strip()
        seq = next(f).strip()
        next(f)
        next(f)
        yield name, seq


def read_fasta(f):
    name = None
    seq = None
    for line in f:
        if line.startswith('>'):
            if name:
                yield name, seq
            name = line.lstrip('>').strip()
            seq = ''
        else:
            seq += line.strip()
    yield name, seq


def read(f):
    if f.name.endswith('a'):
        return read_fasta(f)
    else:
        return read_fastq(f)


def main():
    parser = argparse.ArgumentParser(description='De Bruijn graph')
    parser.add_argument('-i', '--input', help='Input fastq', metavar='File',
                        type=argparse.FileType(), required=True)
    parser.add_argument('-k', help='k-mer size (default: 55)', metavar='Int',
                        type=int, default=55)
    parser.add_argument('-o', '--output', help='Output dot', metavar='File', required=True)
    parser.add_argument('-c', '--compress', help='Shrink graph', default=True)
    parser.add_argument('--vertex', help='Show vertex sequences', action='store_true')
    parser.add_argument('--edge', help='Show edge sequences', action='store_true')
    args = parser.parse_args()

    Graph.k = args.k
    Vertex.show_sequences = args.vertex
    Edge.show_sequences = args.edge

    graph = Graph()
    for name, seq in tqdm(read(args.input)):
        graph.add_seq(seq)
        graph.add_seq(reverse_complement(seq))

    if args.compress:
        graph.compress()
    graph.save_dot(args.output)


if __name__ == '__main__':
    main()
