Read("sr_class2.g");;
LoadPackage("smallgrp");;
SR2_ExportTwins:=function()
 local id,G,D,B,Q,mask,i,j,k,a;
 SR2_PrintTo("twin_tripwire_forms.txt","");
 for id in [6249567,6249622,6249623] do
  G:=SmallGroup(512,id);D:=SR2_Extract(G);B:=[];
  for a in [1..D.r] do
   mask:=0;k:=0;
   for i in [1..D.m] do for j in [i+1..D.m] do
    mask:=mask+IntFFE(D.B[a][i][j])*2^k;k:=k+1;
   od;od;
   Add(B,mask);
  od;
  Q:=Sum([1..D.m*D.r],i->IntFFE(Flat(D.Q)[i])*2^(i-1));
  SR2_AppendTo("twin_tripwire_forms.txt",id," ",D.m," ",D.r," ",Q);
  for a in B do SR2_AppendTo("twin_tripwire_forms.txt"," ",a);od;
  SR2_AppendTo("twin_tripwire_forms.txt","\n");
 od;
end;
SR2_ExportTwins();
QUIT;
