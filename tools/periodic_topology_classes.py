"""Exact isomorphism of periodic cubical flag graphs.

Lattice offsets identify flags across seams, then are omitted from labels.
Flags preserve attaching maps when quotient cells repeat vertices or edges. Geometry, vertex numbering, deck gauge and the chosen tile
cut are irrelevant. Bottom and top boundary roles are preserved. Matching is
combinatorial equivalence, not general homeomorphism of the underlying slab.
"""
from collections import Counter
import numpy as np
from mesh_tools import EDGES, FACES
from periodic_transition import labels, key, incidence


def quotient_graph(p, h, o):
    # Flags (cell, face, edge, corner) retain attachment maps even where a
    # quotient cell repeats a vertex or edge. A plain Hasse diagram would lose
    # that information. Four edge colors encode changing one flag component.
    kinds=[];adj=[]
    from collections import defaultdict
    inside=[defaultdict(list) for _ in range(3)];across=defaultdict(list)
    physical=incidence(h,o)
    for ci,cell in enumerate(labels(h,o)):
        for fi,indices in enumerate(FACES):
            face=cell[indices];fk=key(face);anchor=min(map(tuple,face))
            boundary=len(physical[fk])==1
            kind=('bottom' if abs(float(np.mean(p[face[:,0],2]))-float(p[:,2].min()))<1e-9 else 'top') if boundary else 'interior'
            for j,corner in enumerate(indices):
                for edge in [(indices[j-1],corner),(corner,indices[(j+1)%4])]:
                    edge=tuple(sorted(map(int,edge)));v=len(kinds);kinds.append(kind);adj.append({})
                    inside[0][(ci,fi,edge)].append(v)
                    inside[1][(ci,fi,int(corner))].append(v)
                    inside[2][(ci,edge,int(corner))].append(v)
                    def shifted(index):
                        q=cell[index];return int(q[0]),int(q[1]-anchor[1]),int(q[2]-anchor[2])
                    across[(fk,shifted(corner),tuple(sorted(shifted(i) for i in edge)))].append(v)
    for color,groups in zip((1,2,4,8),[*inside,across]):
        for members in groups.values():
            if len(members)==1:
                if color!=8:raise ValueError('Broken cell flag')
                continue
            if len(members)!=2:raise ValueError('Nonmanifold flag adjacency')
            a,b=members;adj[a][b]=color;adj[b][a]=color
    return kinds,adj


def isomorphism(a,b):
    """Return a verified full node bijection, or None after exact rejection."""
    ak,aa=a;bk,ba=b
    if len(ak)!=len(bk) or Counter(ak)!=Counter(bk):return None
    palette={k:i for i,k in enumerate(sorted(set(ak)))}
    ac=[palette[k] for k in ak];bc=[palette[k] for k in bk]
    def refine(ac,bc):
        while True:
            signatures=[(colors[v],tuple(sorted((colors[w],m) for w,m in adjacency[v].items())))
                        for colors,adjacency in [(ac,aa),(bc,ba)] for v in range(len(colors))]
            names={s:i for i,s in enumerate(sorted(set(signatures)))}
            new=[names[s] for s in signatures];an=new[:len(ac)];bn=new[len(ac):]
            if Counter(an)!=Counter(bn):return None
            if len(set(an))==len(set(ac)):return an,bn
            ac,bc=an,bn
    def solve(ac,bc):
        refined=refine(ac,bc)
        if refined is None:return None
        ac,bc=refined;counts=Counter(ac)
        if len(counts)==len(ac):
            inverse={color:v for v,color in enumerate(bc)};mapping=[inverse[color] for color in ac]
            if all({mapping[w]:m for w,m in aa[v].items()}==ba[mapping[v]] for v in range(len(ac))):return mapping
            return None
        color=min((c for c in counts if counts[c]>1),key=lambda c:counts[c])
        v=ac.index(color);unique=max(ac)+1
        for w,c in enumerate(bc):
            if c!=color:continue
            an=ac.copy();bn=bc.copy();an[v]=bn[w]=unique
            mapping=solve(an,bn)
            if mapping is not None:return mapping
        return None
    return solve(ac,bc)


def classify(named_meshes):
    groups=[];graphs={};witnesses={}
    for name,mesh in named_meshes:
        graph=quotient_graph(*mesh[:3]);graphs[name]=graph
        for group in groups:
            mapping=isomorphism(graph,graphs[group[0]])
            if mapping is not None:
                group.append(name);witnesses[name]=dict(reference=group[0],node_mapping=mapping);break
        else:
            groups.append([name]);witnesses[name]=dict(reference=name,node_mapping=list(range(len(graph[0]))))
    return groups,witnesses
