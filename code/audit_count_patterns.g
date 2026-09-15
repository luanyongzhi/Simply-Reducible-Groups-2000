# Structural audit of every supplied record at orders congruent to 4 modulo 8.
# The proof of exhaustive coverage is the SmallGroups census; this checks records.
LoadPackage("smallgrp");;
Read("count_v2_2_ids.g");;
Read("sr_class2.g");;
AuditCountPatterns := function()
local pair,G,M,H,exceptions,cnt;
exceptions:=[]; cnt:=0;
for pair in patternids do
 G:=SmallGroup(pair);;
 M:=HallSubgroup(G,Filtered(Set(FactorsInt(Size(G))),p->p<>2));;
 if Size(M)<>Size(G)/4 or not IsNormal(G,M) then Error("normal Hall failure",pair);fi;
 if not IsElementaryAbelian(FactorGroup(G,M)) then Error("V4 quotient failure",pair);fi;
 if not IsAbelian(M) then
  Add(exceptions,pair);
  Print("EXCEPTION ",pair," HallOrder=",Size(M)," HallStructure=",StructureDescription(M)," SR=",SR2_CharacterCheck(G),"\n");
 fi;
 cnt:=cnt+1;
od;
if exceptions<>[[300,25],[1500,37],[1620,422]] then Error("wrong exceptions");fi;
Print("COMPLETE ",cnt," records: all normal odd Hall and V4 quotient; ",cnt-Length(exceptions)," abelian Hall; exceptions=",exceptions,"\n");
G:=AlternatingGroup(5);;
Print("A5 SR=",SR2_CharacterCheck(G),"\n");
G:=DihedralGroup(16);;
Print("D16 exponent=",Exponent(G)," SR=",SR2_CharacterCheck(G),"\n");
end;;
AuditCountPatterns();
QUIT;
