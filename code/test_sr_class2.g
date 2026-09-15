# Run from the directory containing sr_class2.g:
# gap -q -b test_sr_class2.g
Read("sr_class2.g");;
if LoadPackage("smallgrp")=fail then Error("SmallGrp is required for this test"); fi;
Print("GAP ",GAPInfo.Version,"; SmallGrp ",PackageInfo("smallgrp")[1].Version,"\n");
SR2_tested:=0;; SR2_roundtrips:=0;;
for SR2_n in [8,16,32,64] do
    SR2_count:=0; SR2_stems:=0;
    for SR2_id in [1..NrSmallGroups(SR2_n)] do
        SR2_G:=SmallGroup(SR2_n,SR2_id);
        SR2_X:=SR2_Extract(SR2_G);
        if SR2_X<>fail then
            SR2_got:=SR2_IsSRData(SR2_X.B,SR2_X.Q);
            SR2_want:=SR2_CharacterCheck(SR2_G);
            if SR2_got<>SR2_want then Error("SR criterion mismatch",[SR2_n,SR2_id]); fi;
            SR2_H:=SR2_Group(SR2_X.B,SR2_X.Q);
            if Size(SR2_H)<>SR2_n or IdGroup(SR2_H)<>[SR2_n,SR2_id] then
                Error("extract/reconstruct mismatch",[SR2_n,SR2_id]);
            fi;
            SR2_tested:=SR2_tested+1; SR2_roundtrips:=SR2_roundtrips+1;
            if SR2_got then
                SR2_count:=SR2_count+1;
                if SR2_X.stem then SR2_stems:=SR2_stems+1; fi;
                SR2_P:=SR2_PolarData(SR2_X.B);
                SR2_A:=SR2_RealityAffine(SR2_P);
                if SR2_A=fail then Error("missing known reality solution"); fi;
                if SR2_A.equations<>[] and
                   SR2_A.equations*Flat(SR2_X.Q)<>SR2_A.rhs then
                    Error("known refinement fails the linear system");
                fi;
                if SR2_X.stem and SR2_X.r=2 and SR2_A.number<>2^SR2_X.m then
                    Error("r=2 refinement-dimension theorem failed");
                fi;
            fi;
        elif NilpotencyClassOfGroup(SR2_G)=2 and SR2_CharacterCheck(SR2_G) then
            Error("SR class-2 group fails extraction hypotheses");
        fi;
    od;
    Print("order=",SR2_n," class2_SR=",SR2_count," stem_SR=",SR2_stems,"\n");
od;
# Every refinement of the standard 2-dimensional symplectic form is SR.
SR2_B:=[[[0,1],[1,0]]];;
SR2_refcount:=0;;
SR2_stream:=SR2_ForEachRefinement(SR2_B,function(B,Q)
    local G;
    G:=SR2_Group(B,Q);
    if not SR2_IsSRData(B,Q) or not SR2_CharacterCheck(G) then
        Error("streamed candidate fails definition-level SR test");
    fi;
    SR2_refcount:=SR2_refcount+1;
end);;
if SR2_stream<>4 or SR2_refcount<>4 then Error("refinement streaming failure"); fi;
Print("PASS: ",SR2_tested," exact-character comparisons; ",SR2_roundtrips,
      " reconstruction checks; 4 streamed refinements.\n");
QUIT;
