#!/usr/bin/env python3
"""Prepare the portable GAP verification/export for all new stem strata."""
import argparse, json
from pathlib import Path

def matrices(m,B):
    mats=[]
    for mask in B:
        M=[[0]*m for _ in range(m)];k=0
        for i in range(m):
            for j in range(i+1,m):
                M[i][j]=M[j][i]=(mask>>k)&1;k+=1
        mats.append(M)
    return mats

def main():
    ap=argparse.ArgumentParser();ap.add_argument('json',nargs='+',type=Path)
    ap.add_argument('--output',type=Path,default=Path('verify_special_representatives.g'));a=ap.parse_args()
    lines=['Read("sr_class2.g");;','Read("sr1024_class2_m8r2_representatives.g");;',
      'SR2_special_candidates:=List(SR2_m8r2,R->rec(m:=8,r:=2,B:=R.B,Q:=R.Q,family:=R.family));;']
    for path in a.json:
        data=json.loads(path.read_text())
        assert data['complete_for_supplied_input']
        for polar in data['polar_orbits']:
            m,r=polar['m'],polar['r'];B=matrices(m,polar['B'])
            for orbno,orb in enumerate(polar['quadratic_orbits'],1):
                Q=[[(orb['Q']>>(i*r+j))&1 for j in range(r)] for i in range(m)]
                family=f"m{m}r{r}_polar{polar['polar_orbit']}_q{orbno}"
                lines.append(f'Add(SR2_special_candidates,rec(m:={m},r:={r},B:={B},Q:={Q},family:="{family}",orbit_size:={orb["size"]}));')
    lines.append('''
SR2_VerifySpecial:=function()
local num,R,G,H,Xdata,n,roots,classes,s2,s3,c3,inv,k,os;
num:=0;
SR2_PrintTo("sr1024_special_representatives.g","# Complete stem list; coverage proved in accompanying manuscript.\\nSR2_special_representatives:=[\\n");
SR2_PrintTo("sr1024_special_representatives.csv","No,m,r,Family,Order,Class,CenterSize,DerivedSize,Involutions,k,WignerS2,WignerS3,PcCode\\n");
for R in SR2_special_candidates do
 num:=num+1;G:=SR2_Group(R.B,R.Q);n:=Size(G);
 if n<>1024 or NilpotencyClassOfGroup(G)<>2 or Size(Centre(G))<>2^R.r or Size(DerivedSubgroup(G))<>2^R.r then
  Error("special group metadata mismatch",num);fi;
 if not SR2_IsSRData(R.B,R.Q) then Error("form criterion failed",num);fi;
 roots:=Collected(List(Elements(G),x->x^2));classes:=ConjugacyClasses(G);k:=Length(classes);
 s2:=Sum(roots,x->x[2]^2);s3:=Sum(roots,x->x[2]^3);c3:=Sum(classes,C->n^2/Size(C));
 if s2<>n*k or s3<>c3 then Error("Wigner failed",num);fi;
 inv:=Number(Elements(G),x->Order(x)=2);R.code:=CodePcGroup(G);R.number:=num;
 H:=PcGroupCode(R.code,1024);Xdata:=SR2_Extract(H);
 if Size(H)<>1024 or Size(Centre(H))<>2^R.r or Xdata=fail or not SR2_IsSRData(Xdata.B,Xdata.Q) then
  Error("saved pc code roundtrip failed",num);fi;
 SR2_AppendTo("sr1024_special_representatives.g",R,",\\n");
 SR2_AppendTo("sr1024_special_representatives.csv",num,",",R.m,",",R.r,",",R.family,",1024,2,",2^R.r,",",2^R.r,",",inv,",",k,",",s2,",",s3,",",R.code,"\\n");
 Print("PASS ",num," ",R.family," inv=",inv," k=",k,"\\n");
od;
SR2_AppendTo("sr1024_special_representatives.g","];;\\n");
Print("COMPLETE ",num," special representatives; all metadata, form criterion, independent Wigner S2/S3 and saved pc-code roundtrips passed; GAP ",GAPInfo.Version,".\\n");
end;
SR2_VerifySpecial();
QUIT;
''')
    a.output.write_text('\n'.join(lines))

if __name__=='__main__':main()
