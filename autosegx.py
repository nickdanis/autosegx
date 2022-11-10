import pydot
import re
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
from networkx.utils.misc import graphs_equal
import networkx.algorithms.isomorphism as iso
from itertools import chain, combinations
from collections import defaultdict
from functools import cached_property


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
        if list(self.edges()) == []:
            return str(list(self.nodes()))
        else:
            return str(list(self.edges()))

    def label_nodes(self, raw_edges):
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
        nx.draw_networkx(self,graphviz_layout(self,prog='dot'),node_color='white')
        

class Segment(Geometry):
    def __init__(self, ipa=''):
        Geometry.__init__(self)
        self.ipa = ipa

class Theory():
    def __init__(self, theory_dict = dict()):
        self.name = ''
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
        all_factors = [fac for seg in self.segments.values() for fac in seg.factors]
        unique_factors = []
        for fac in all_factors:
            if len(unique_factors) == 0:
                unique_factors.append(fac)
            else:
                equal_to = []
                for f in unique_factors:
                    if graphs_equal(f,fac):
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


def compare_theories(t1,t2,verbose=0):
    unique_t1 = list(t1.nce.keys() - t2.nce.keys())
    unique_t2 = list(t2.nce.keys() - t1.nce.keys())
    shared = set(t1.nce.keys()).intersection(t2.nce.keys())
    if unique_t1 == unique_t2 == []:
        nc_preserving = True
    else:
        nc_preserving = False
    if verbose >= 1:
        c1 = 12
        c2 = 10
        c3 = 10
        if nc_preserving:
            print("The theories are natural class preserving")
        else:
            print("The theories are NOT natural class preserving")
        print(f"{'nat. classes':>{c1}}{'T1':^{c2}}{'T2':^{c3}}")
        print(f"{'unique':>{c1}}{len(unique_t1):^{c2}}{len(unique_t2):^{c3}}")
        print(f"{'shared':>{c1}}{len(shared):^{c2}}{len(shared):^{c3}}")
        print(f"{'total':>{c1}}{len(t1.nce.keys()):^{c2}}{len(t2.nce.keys()):^{c3}}")
    if verbose >= 2:
        if unique_t1 != []:
            print(f"Natural classes unique to T1:")
            for nce in unique_t1:
                print(nce)
                for fac in t1.nce[nce]:
                    print(fac)
        if unique_t2 != []:
            print(f"Natural classes unique to T2:")
            for nce in unique_t2:
                print("\t",nce)
                for fac in t2.nce[nce]:
                    print("\t\t",fac)
    return nc_preserving



