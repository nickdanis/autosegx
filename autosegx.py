import pydot
import re
import graphviz
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
from networkx.utils.misc import graphs_equal
import networkx.algorithms.isomorphism as iso
from itertools import chain, combinations
from collections import defaultdict
from functools import cached_property, cache


def powerset(iterable):
    # from https://stackoverflow.com/questions/1482308/how-to-get-all-subsets-of-a-set-powerset
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def part_of_class(seg, natclass):
    '''takes a segment and a factor
    returns true if the induced subgraph of the segment with the nodes of that factor
    is isomorphic to the natural class factor'''
    nm = iso.categorical_node_match("label",'')
    for fac in seg.factors:
        if nx.is_isomorphic(fac,natclass,node_match=nm):
            return True
    return False



class Geometry(nx.DiGraph):
    def __init__(self):
        nx.DiGraph.__init__(self)
        self.raw_nodes = {}
        self.labeled_nodes = {}
        self._factors = []
    
    def __str__(self):
        graph_string = ""
        labels = dict(self.nodes(data="label"))
        if list(self.edges()) != []:
            nice_edges = [f"{a} -> {labels[b]}" for a, b in self.edges()]
            graph_string += ', '.join(nice_edges)
        else:
            nice_nodes = [f"{labels[n]}" for n in self.nodes()]
            graph_string += ', '.join(nice_nodes)
        return graph_string

    def label_nodes(self, raw_edges):
        '''converts node to explicit label, removes any digits
        so n1 and n2 will be unique nodes, but will share the label 'n' '''
        self.raw_nodes = list(set(chain(*raw_edges)))
        labeled_nodes = [(node, {'label' : re.sub(r'\d+','',node)}) for node in self.raw_nodes]
        self.labeled_nodes = labeled_nodes
        return labeled_nodes

    def create_geometry(self, raw_edges):
        self.label_nodes(raw_edges)
        self.add_nodes_from(self.labeled_nodes)
        self.add_edges_from(list(raw_edges)) 

    @cached_property  
    def factors(self):
        '''calculates all connected factors of the representation
        results are cached, and accessed through calling factors again'''
        factor_nodes = [sub for sub in powerset(self.nodes) if len(sub)>0]
        factors = []
        for nodes in factor_nodes:
            factor = nx.subgraph(self,nodes)
            if nx.is_connected(factor.to_undirected()):
                factors.append(factor)
        self._factors = factors
        return self._factors

    def draw_phono(self):
        '''draws graph as a tree'''
        nx.draw_networkx(self,graphviz_layout(self,prog='dot'),node_color='white',with_labels=True)

    def gv(self):
        title = f"[{self.ipa}]" if self.ipa != '' else ''
        graph_attr = {'splines':'false', 'label':title, 'labelloc':'t'}
        node_attr={'shape': 'plaintext'}
        edge_attr={'arrowhead' : 'none', 'tailport':'s', 'headport':'n'}
        labels = dict(self.nodes(data="label"))
        dot = graphviz.Digraph(graph_attr = graph_attr, node_attr=node_attr, edge_attr=edge_attr)
        for node in self.nodes():
            dot.node(node,labels[node])
        for edge in self.edges():
            dot.edge(edge[0],edge[1])
        return dot
        

class Segment(Geometry):
    def __init__(self, ipa=''):
        Geometry.__init__(self)
        self.ipa = ipa

class Theory():
    def __init__(self, theory_dict = dict(), name = ''):
        self.name = name
        self.segments = dict()
        self._factors = []
        self._natural_classes = defaultdict(list)
        self._nce = defaultdict(list)
        self.theory_dict = theory_dict
        for seg, rep in self.theory_dict.items():
            G = Segment()
            G.ipa = seg
            G.create_geometry(rep)
            self.segments[seg] = (G)

    @cached_property  
    def factors(self):
        '''todo: remove duplicate factors based on label name
        '''
        nm = iso.categorical_node_match("label",'')
        all_factors = [fac for seg in self.segments.values() for fac in seg.factors]
        unique_factors = []
        for fac in all_factors:
            if len(unique_factors) == 0:
                unique_factors.append(fac)
            else:
                equal_to = []
                for f in unique_factors:
                    # if graphs_equal(f,fac):
                    if nx.is_isomorphic(f,fac,node_match=nm):
                        equal_to.append(f)
                if len(equal_to) == 0:
                    unique_factors.append(fac)
        self._factors = unique_factors
        return self._factors

    @cached_property
    def natural_classes(self):
        for fac in self.factors:
            for seg in self.segments.values():
                if part_of_class(seg,fac):
                    self._natural_classes[fac].append(seg.ipa)
        return self._natural_classes

    @cached_property
    def nce(self):
        for key, value in self.natural_classes.items():
            self._nce[tuple(sorted(value))].append(key)
        return self._nce


class Comparison():
    def __init__(self, theories):
        self.theories = theories
        self.t1 = theories[0]
        self.t2 = theories[1]
        self._unique = dict()
        self._shared = set()
        self._preserving = bool
        self._results = dict()

    def unique_to(self, t1, t2):
        # return list(t1.nce.keys() - t2.nce.keys())
        return { k : t1.nce[k] for k in set(t1.nce) - set(t2.nce) }
    
    @cached_property
    def unique(self):
        self._unique[self.t1] = self.unique_to(self.t1, self.t2)
        self._unique[self.t2] = self.unique_to(self.t2, self.t1)
        return self._unique
    
    @cached_property
    def shared(self):
        self._shared = { k : (self.t1.nce[k],self.t2.nce[k]) for k in set(self.t1.nce).intersection(self.t2.nce) }
        return self._shared
    
    @cached_property
    def preserving(self):
        if self.unique[self.t1] == self.unique[self.t2] == {}:
            self._preserving = True
        else:
            self._preserving = False
        return self._preserving
    
    def print_classes(self, classes, verbose):
        for idx, nce in enumerate(classes):
            print(f"\t{idx+1}. [ {' '.join(nce)} ]")
            if verbose >= 2:
                print(f"\t\tDefining factors:")
                for idx, factor in enumerate(classes[nce]):
                    if type(factor) == list:
                        for subidx, seg in enumerate(factor):
                            print(f"\t\t({idx+1}-{subidx+1}) {seg}")
                    else:
                        print(f"\t\t({idx+1}) {factor}")

    def results(self, verbose = 0, shared=True):
        print(f"The theories are {'' if self.preserving else 'NOT '}natural class preserving.")
        if not self.preserving:
            for theory in self.unique.keys():
                print(f"\n{len(self.unique[theory])} natural class(es) unique to {theory.name}.")
                if verbose >=1:
                    self.print_classes(self.unique[theory],verbose)
        print(f"\n{len(self.shared)} natural class(es) shared.")
        if verbose >=1 and shared:
            self.print_classes(self.shared,verbose)




