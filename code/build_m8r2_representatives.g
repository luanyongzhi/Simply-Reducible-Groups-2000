# Build 18 representatives across the four admissible dimension-eight pencils.
# Requires classify_last_pencil.g to have produced its five representatives.
Read("sr_class2.g");;
Read("last_pencil_representatives.g");;
SR2_WignerMoments:=function(G)
 local roots,classes,n,s2,s3,c3;
 n:=Size(G);roots:=Collected(List(Elements(G),x->x^2));classes:=ConjugacyClasses(G);
 s2:=Sum(roots,x->x[2]^2);s3:=Sum(roots,x->x[2]^3);
 c3:=Sum(classes,C->n^2/Size(C));
 return rec(ok:=s2=n*Length(classes) and s3=c3,s2:=s2,s3:=s3,
     k:=Length(classes),centralizer_square_sum:=c3);
end;
SR2_BuildM8R2:=function()
 local reps,spec,counts,types,B,Q,F,colors,start,color,plane,positions,i,a,
       t,H,G,W,P,recdata,num,inv;
 F:=GF(2);colors:=[[1,0],[0,1],[1,1]]*One(F);reps:=[];
 for spec in [
  rec(name:="R310",counts:=[3,1,0],types:=[[0,0,0],[0,1,0],[1,0,0],[1,1,0]]),
  rec(name:="R220",counts:=[2,2,0],types:=[[0,0,0],[0,1,0],[1,1,0]]),
  rec(name:="R211",counts:=[2,1,1],types:=[[0,0,0],[0,0,1],[0,1,1],[1,0,0],[1,0,1],[1,1,1]])
 ] do
  B:=List([1,2],a->NullMat(8,8,F));positions:=[];start:=1;
  for color in [1..3] do
   Add(positions,start);
   for plane in [1..spec.counts[color]] do
    for a in [1,2] do B[a][start][start+1]:=colors[color][a];B[a][start+1][start]:=colors[color][a];od;
    start:=start+2;
   od;
  od;
  for types in spec.types do
   Q:=NullMat(8,2,F);
   for color in [1..3] do
    if types[color]=1 then
     Q[positions[color]]:=ShallowCopy(colors[color]);Q[positions[color]+1]:=ShallowCopy(colors[color]);
    fi;
   od;
   Add(reps,rec(family:=spec.name,arf_types:=types,B:=B,Q:=Q));
  od;
 od;
 for recdata in SR2_last_pencil_representatives do
  Add(reps,rec(family:="S1S1R1",arf_types:=[],B:=recdata.B,Q:=recdata.Q));
 od;
 SR2_PrintTo("sr1024_class2_m8r2.csv","No,Family,Order,Class,DerivedSize,CenterSize,Involutions,k,RankProfile,WignerS2,WignerS3,PcCode\n");
 SR2_PrintTo("sr1024_class2_m8r2_representatives.g","# GAP ",GAPInfo.Version,". Exactly 18 (m,r)=(8,2) stem SR representatives.\nSR2_m8r2:=[\n");
 num:=0;
 for recdata in reps do
  num:=num+1;G:=SR2_Group(recdata.B,recdata.Q);P:=SR2_PolarData(recdata.B);
  if Size(G)<>1024 or NilpotencyClassOfGroup(G)<>2 or Size(DerivedSubgroup(G))<>4 or Size(Centre(G))<>4 then
   Error("incorrect group invariants");fi;
  if not SR2_IsSRData(recdata.B,recdata.Q) then Error("form SR test failed");fi;
  W:=SR2_WignerMoments(G);
  if not W.ok then Error("independent Wigner moment check failed");fi;
  inv:=Number(Elements(G),x->Order(x)=2);
  recdata.code:=CodePcGroup(G);recdata.number:=num;
  SR2_AppendTo("sr1024_class2_m8r2.csv",num,",",recdata.family,",1024,2,4,4,",inv,",",W.k,",\"",SortedList(P.ranks{[2..4]}),"\",",W.s2,",",W.s3,",",recdata.code,"\n");
  SR2_AppendTo("sr1024_class2_m8r2_representatives.g",recdata,",\n");
  Print("PASS ",num," ",recdata.family," inv=",inv," k=",W.k,"\n");
 od;
 SR2_AppendTo("sr1024_class2_m8r2_representatives.g","];;\n");
 Print("COMPLETE: all ",num," representatives independently pass both Wigner moment identities in GAP ",GAPInfo.Version,".\n");
end;
SR2_BuildM8R2();
QUIT;
