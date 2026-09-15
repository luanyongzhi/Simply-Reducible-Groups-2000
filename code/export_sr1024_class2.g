# export_sr1024_class2.g
# Reconstruct and independently validate all 221 class-two SR representatives
# of order 1024. This exports supplied exact orbit representatives; it does
# not repeat the expensive orbit classification or deduplicate by invariants.
# Run from the directory containing this script; see README_sr1024_class2.md.

SR1024Class2_FindFile := function(name)
    local path;
    for path in [name, Concatenation("../data/",name),
                 Concatenation("sr_review/",name),
                 Concatenation("../upload/",name), Concatenation("upload/",name)] do
        if IsExistingFile(path) then return path; fi;
    od;
    Error("Missing input file: ",name);
end;

Read(SR1024Class2_FindFile("sr_class2.g"));;
Read(SR1024Class2_FindFile("sr1024_special_representatives.g"));;
if LoadPackage("smallgrp")=fail then Error("SmallGrp is required"); fi;
if not IsBound(SR1024Class2_Input512) then
    SR1024Class2_Input512:=SR1024Class2_FindFile("SR_groups_results512.csv");
fi;
if not IsBound(SR1024Class2_Output) then
    SR1024Class2_Output:="sr1024_class2.csv";
fi;

# Only the first two CSV fields are read. Later fields may contain quoted
# commas in StructureDescription and are deliberately not parsed here.
SR1024Class2_Read512Ids := function(path)
    local text,lines,line,fields,ids,id;
    text:=StringFile(path);
    if text=fail then Error("Cannot read order-512 census"); fi;
    lines:=SplitString(text,"\n",""); ids:=[];
    for line in lines do
        fields:=SplitString(line,",","");
        if Length(fields)>=2 and Int(NormalizedWhitespace(fields[1]))=512 then
            id:=Int(NormalizedWhitespace(fields[2]));
            if id=fail then Error("Invalid order-512 SmallGrp identifier"); fi;
            Add(ids,id);
        fi;
    od;
    if Length(ids)<>317 or Length(Set(ids))<>317 then
        Error("Expected 317 distinct records in the order-512 census");
    fi;
    return Set(ids);
end;

SR1024Class2_Export := function()
    local ids,parents,ea512,id,K,R,G,H,D,Z,X,Y,m,r,hasc2,code,oldcode,
          data,lines,rows,num,nonspecial,special,classes,k,squares,s2,s3,
          inv,porder,pid,source,sourceid,family,counts,item,expected;
    ids:=SR1024Class2_Read512Ids(SR1024Class2_Input512);
    parents:=[]; ea512:=fail;
    for id in ids do
        K:=SmallGroup(512,id);
        if IsElementaryAbelian(K) then ea512:=id; fi;
        if NilpotencyClassOfGroup(K)=2 then Add(parents,id); fi;
    od;
    if Length(parents)<>95 or ea512=fail then
        Error("Expected 95 class-two parents and one elementary abelian group");
    fi;
    if Length(SR2_special_representatives)<>126 then
        Error("Expected 126 classified special representatives");
    fi;
    # Keep the existing 803-row catalogue's class-two order: padding, special.
    data:=[];
    for id in parents do
        Add(data,rec(source:="class2_padding",sourceid:=id,
            family:="C2_direct_product",special:=false));
    od;
    for R in SR2_special_representatives do
        Add(data,rec(source:="class2_special",sourceid:=R.number,
            family:=R.family,special:=true,representative:=R));
    od;
    lines:=[Concatenation(
        "No,ParentOrder,ParentId,Step,Class,Exponent,CenterSize,",
        "DerivedLength,NrConjClasses,HasC2DirectFactor,Code,",
        "Order,DerivedSize,Involutions,m,r,IsSpecial,Source,SourceID,",
        "Family,WignerS2,WignerS3\n")];
    rows:=[]; num:=0; nonspecial:=0; special:=0; counts:=[];
    for item in data do
        source:=item.source; sourceid:=item.sourceid; family:=item.family;
        if item.special then
            R:=item.representative;
            G:=SR2_Group(R.B,R.Q);
            if not SR2_IsSRData(R.B,R.Q) then Error("Input form is not SR"); fi;
        else
            G:=DirectProduct(SmallGroup(512,sourceid),CyclicGroup(2));
        fi;
        if Size(G)<>1024 or NilpotencyClassOfGroup(G)<>2 then
            Error("Unexpected order or nilpotency class");
        fi;
        D:=DerivedSubgroup(G); Z:=Centre(G); X:=SR2_Extract(G);
        if X=fail or not SR2_IsSRData(X.B,X.Q) then
            Error("Extracted class-two SR criterion failed");
        fi;
        m:=X.m; r:=X.r; hasc2:=Size(Z)>Size(D);
        if item.special<>(not hasc2) or m+r<>10 then
            Error("Special/non-special split or dimensions failed");
        fi;
        if item.special and (m<>R.m or r<>R.r) then
            Error("Input representative dimension mismatch");
        fi;
        # An independent SR check using exact square-root moments.
        classes:=ConjugacyClasses(G); k:=Length(classes);
        squares:=Collected(List(Elements(G),x->x^2));
        s2:=Sum(squares,x->x[2]^2); s3:=Sum(squares,x->x[2]^3);
        if s2<>1024*k or s3<>Sum(classes,C->1024^2/Size(C)) then
            Error("Independent Wigner SR identities failed");
        fi;
        inv:=Number(Elements(G),x->x^2=One(G))-1;
        code:=CodePcGroup(G);
        if item.special and code<>R.code then
            Error("Special pc presentation differs from saved representative");
        fi;
        # PcGroupCode decodes this particular presentation. CodePcGroup is
        # not an isomorphism-canonical group identifier.
        H:=PcGroupCode(code,1024); Y:=SR2_Extract(H);
        if CodePcGroup(H)<>code or Size(H)<>1024 or
           NilpotencyClassOfGroup(H)<>2 or Y=fail or
           not SR2_IsSRData(Y.B,Y.Q) or
           Size(Centre(H))<>Size(Z) or Size(DerivedSubgroup(H))<>Size(D) then
            Error("Saved presentation reconstruction failed");
        fi;
        # Canonical lower exponent-2 parent is G/G' = C2^m; this is not
        # the order-512 direct-product source recorded in SourceID.
        porder:=2^m;
        if porder=512 then pid:=ea512;
        else pid:=IdGroup(ElementaryAbelianGroup(porder))[2]; fi;
        num:=num+1;
        if item.special then special:=special+1; Add(counts,[m,r]);
        else nonspecial:=nonspecial+1; fi;
        Add(rows,[num,porder,pid,r,2,Exponent(G),Size(Z),DerivedLength(G),
            k,hasc2,code,1024,Size(D),inv,m,r,item.special,source,sourceid,
            family,s2,s3]);
        Print("PASS ",num," ",source," ",sourceid," m=",m," r=",r,
              " involutions=",inv,"\n");
    od;
    expected:=[[[6,4],28],[[7,3],80],[[8,2],18]];
    if num<>221 or special<>126 or nonspecial<>95 or
       Collected(counts)<>expected then Error("Final catalogue tally failed"); fi;
    if Length(Set(List(rows,row->row[11])))<>221 then
        Error("Duplicate presentation code in supplied representatives");
    fi;
    for R in rows do
        Add(lines,Concatenation(JoinStringsWithSeparator(List(R,String),","),"\n"));
    od;
    # Publish the CSV only after every record has passed all checks.
    SR2_PrintTo(SR1024Class2_Output,Concatenation(lines));
    Print("COMPLETE 221 class-two representatives; 126 special (28,80,18), ",
        "95 non-special; all form tests, Wigner identities, and pc-code ",
        "reconstructions passed.\n");
    return rows;
end;

SR1024Class2_Rows:=SR1024Class2_Export();;
