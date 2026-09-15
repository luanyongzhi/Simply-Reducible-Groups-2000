# Exact orbit certificate for the two 3-dimensional singular blocks plus
# one 2-dimensional regular block. This proves the five Q-orbits for THIS B.
# Completeness among all pencils uses the separate Kronecker-form theorem.
Read("sr_class2.g");;
SR2_ClassifyLastPencil:=function()
local B,F,P,Qs,ids,gens,T,C,s,i,j,a,b,add,changes,gen,M,seen,orbit,queue,
      pos,n,k,Q,Qt,id,reps,inv,tzero,key,invariants,G,code,total;
F:=GF(2); B:=List([1,2],i->NullMat(8,8,F));
for s in [0,3] do
 B[1][s+1][s+2]:=One(F);B[1][s+2][s+1]:=One(F);
 B[2][s+1][s+3]:=One(F);B[2][s+3][s+1]:=One(F);
od;
B[1][7][8]:=One(F);B[1][8][7]:=One(F);
P:=SR2_PolarData(B);Qs:=[];
SR2_ForEachRefinement(B,function(B,Q) Add(Qs,Q);end);
code:=Q->Sum([1..16],i->IntFFE(Flat(Q)[i])*2^(i-1));
ids:=List(Qs,code);gens:=[];
add:=function(changes)
 local T,p;
 T:=IdentityMat(8,F);
 for p in changes do T[p[1]][p[2]]:=T[p[1]][p[2]]+One(F);od;
 Add(gens,rec(T:=T,C:=IdentityMat(2,F)));
end;
# GL(2,2) on the two singular blocks, with dual action on b_i,c_i.
T:=NullMat(8,8,F);
for i in [1..8] do T[i][[4,5,6,1,2,3,7,8][i]]:=One(F);od;
Add(gens,rec(T:=T,C:=IdentityMat(2,F)));
add([[1,4],[5,2],[6,3]]);
# Symmetric shears on the a_i.
for s in [0,3] do add([[s+1,s+2]]);add([[s+1,s+3]]);od;
add([[1,5],[4,2]]);add([[1,6],[4,3]]);
# The symplectic group on the regular block.
T:=IdentityMat(8,F);T[7][7]:=Zero(F);T[8][8]:=Zero(F);
T[7][8]:=One(F);T[8][7]:=One(F);
Add(gens,rec(T:=T,C:=IdentityMat(2,F)));add([[8,7]]);
# These mix regular and singular blocks while preserving both forms.
for s in [0,3] do
 add([[s+1,7],[8,s+2]]);add([[s+1,8],[7,s+2]]);
od;
# Interchange the two rank-six forms; the rank-four form is fixed.
T:=IdentityMat(8,F);T[3][2]:=One(F);T[6][5]:=One(F);
C:=[[1,0],[1,1]]*One(F);Add(gens,rec(T:=T,C:=C));
# Check each claimed generator against the full tensor, not a fingerprint.
for gen in gens do
 if RankMat(gen.T)<>8 or RankMat(gen.C)<>2 then Error("singular generator");fi;
 for a in [1,2] do
  M:=NullMat(8,8,F);
  for b in [1,2] do M:=M+gen.C[b][a]*(gen.T*B[b]*TransposedMat(gen.T));od;
  if M<>B[a] then Error("generator does not preserve B");fi;
 od;
od;
Print("GAP ",GAPInfo.Version,"; checked polar generators=",Length(gens),"; Q candidates=",Length(Qs),"\n");
seen:=List(Qs,Q->false);reps:=[];invariants:=[];total:=0;
SR2_PrintTo("last_pencil_orbits.csv","representative,orbit_size,t_zero,involutions,pc_code\n");
SR2_PrintTo("last_pencil_representatives.g","# Complete five-orbit certificate for the fixed pencil in classify_last_pencil.g\nSR2_last_pencil_representatives:=[\n");
for n in [1..Length(Qs)] do
 if seen[n] then continue;fi;
 queue:=[n];seen[n]:=true;pos:=1;
 while pos<=Length(queue) do
  Q:=Qs[queue[pos]];pos:=pos+1;
  for gen in gens do
   Qt:=List(gen.T,v->SR2_QuadraticValue(B,Q,v)*gen.C);
   id:=Position(ids,code(Qt));
   if id=fail then Error("transformation leaves the reality solution space");fi;
   if not seen[id] then seen[id]:=true;Add(queue,id);fi;
  od;
 od;
 Q:=Qs[n];G:=SR2_Group(B,Q);
 inv:=4*Number([0..255],k->IsZero(SR2_QuadraticValue(B,Q,SR2_Bits(k,8))))-1;
 tzero:=IsZero(Q[2]) and IsZero(Q[5]);key:=[tzero,inv];
 if key in invariants then Error("the orbit-separating invariants collide");fi;
 Add(invariants,key);Add(reps,Q);total:=total+Length(queue);
 SR2_AppendTo("last_pencil_orbits.csv",Length(reps),",",Length(queue),",",tzero,",",inv,",",CodePcGroup(G),"\n");
 SR2_AppendTo("last_pencil_representatives.g","rec(B:=",B,",Q:=",Q,",code:=",CodePcGroup(G),"),\n");
 Print("orbit=",Length(reps)," size=",Length(queue)," t_zero=",tzero," involutions=",inv,"\n");
od;
SR2_AppendTo("last_pencil_representatives.g","];;\n");
if total<>256 or Length(reps)<>5 then Error("unexpected orbit count");fi;
Print("COMPLETE: 256 refinements partition into five orbits under verified tensor automorphisms, separated by characteristic invariants.\n");
end;
SR2_ClassifyLastPencil();
QUIT;
