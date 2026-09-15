#!/usr/bin/env python3
"""Export full polar stabilizer actions for finite parent-coset traversal."""
import argparse, json, time
from pathlib import Path
from canonical_class2 import graph, stabilizer_actions

def rref(B, n):
    B = list(B); row = 0
    for p in range(n):
        k = next((k for k in range(row, len(B)) if B[k] >> p & 1), None)
        if k is None: continue
        B[row], B[k] = B[k], B[row]
        for i in range(len(B)):
            if i != row and B[i] >> p & 1: B[i] ^= B[row]
        row += 1
    return tuple(B[:row])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path)
    ap.add_argument('m',type=int);ap.add_argument('output',type=Path);a=ap.parse_args()
    parents={};start=time.monotonic()
    for line in a.input.read_text().splitlines():
        ident,m,r,*B=map(int,line.split())
        if m!=a.m:continue
        parents.setdefault(rref(B,m*(m-1)//2),[]).append(ident)
    metadata=[]
    with a.output.open('w') as out:
        for idx,(B,idents) in enumerate(sorted(parents.items()),1):
            m=a.m;r=len(B)
            acts,order=stabilizer_actions(m,r,B,graph(m,r,B))
            out.write(' '.join(map(str,[m,r,len(acts),*B]))+'\n')
            for T,invS in acts:out.write(' '.join(map(str,T))+'\n')
            metadata.append(dict(parent_number=idx,m=m,r=r,B=B,source_ids=idents,
                                 graph_aut_order_scientific=order,
                                 T_generators=[T for T,invS in acts]))
            out.flush()
            print('PARENT',idx,'B',B,'generators',len(acts),'sec',round(time.monotonic()-start,3),flush=True)
    a.output.with_suffix('.json').write_text(json.dumps(metadata,indent=2))
    print('COMPLETE',len(parents),'parents',flush=True)

if __name__=='__main__':main()
