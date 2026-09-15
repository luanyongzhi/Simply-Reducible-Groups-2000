if IsReadableFile("sr_class2.g") then
  Read("sr_class2.g");
else
  Read("../output/sr_review/sr_class2.g");
fi;;
LoadPackage("smallgrp");;
Print("GAP ",GAPInfo.Version,"\n");
ids:=[6249567, 6249622, 6249623, 7535895, 7535943, 7536102, 7536155, 7536284, 7536911, 7536940, 7536941, 7537438, 7537457, 7537458, 7538135, 7538394, 7538415, 7538416, 7538445, 7540510, 7540544, 7540545, 7577399, 7579142, 7579170, 7579560, 7579947, 7579948, 7579952, 7579957, 7606661, 7606662, 7606711, 7607280, 7640345, 7640799, 7640967, 7644607, 7644708, 7645035, 7645152, 7645402, 7645737, 7646520, 7646542, 7646585, 8840079, 8840080, 8840081, 8841562, 8841573, 8841813, 8841816, 8841960, 10493082, 10493091, 10493092, 10493114, 10493115, 10493131, 10493137, 10493156, 10493169, 10493170, 10493176, 10493220, 10493221, 10493237, 10493238, 10493262, 10493279, 10493293, 10493294, 10493299, 10493360, 10493377, 10493378, 10493526, 10493531, 10493555, 10493560, 10493668, 10493669, 10493714, 10493750, 10493751, 10493754, 10494202, 10494203, 10494205, 10494206, 10494208, 10494209, 10494211, 10494212];;
PrintTo("parents512_forms.txt","");;
RunExtract:=function()
local id,G,X,codes,M,code,bit,i,j;
for id in ids do
  G:=SmallGroup(512,id); X:=SR2_Extract(G); codes:=[];
  for M in X.B do
    code:=0; bit:=0;
    for i in [1..X.m] do for j in [i+1..X.m] do
      if not IsZero(M[i][j]) then code:=code+2^bit; fi;
      bit:=bit+1;
    od; od;
    Add(codes,code);
  od;
  AppendTo("parents512_forms.txt",id," ",X.m," ",X.r);
  for code in codes do AppendTo("parents512_forms.txt"," ",code); od;
  AppendTo("parents512_forms.txt","\n");
od;
end;;
RunExtract();
Print("EXPORTED ",Length(ids)," parents\n");
QUIT;
