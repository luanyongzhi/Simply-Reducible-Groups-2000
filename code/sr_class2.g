# sr_class2.g -- exact tests for class-at-most-two exponent-four 2-groups.
# Target: GAP 4.15.1; core functions require no optional package.
#
# B = [B_1,...,B_r], independent alternating m by m matrices over GF(2).
# Q is m by r and its i-th row is q(e_i).
# q(v)_a = sum_i v_i Q[i][a] + sum_(i<j) v_i v_j B[a][i][j].
# The group has V=G/G' of dimension m and W=G' of dimension r.
# Common radical zero means the group is stem; central C2 factors are allowed.
#
# This file DOES NOT enumerate GL(m,2)-orbits of alternating-form spaces.
# It solves all reality constraints for each supplied B and tests SR exactly.
# Multiple refinements can define isomorphic groups. Candidate counts are NOT
# isomorphism-class counts. See the accompanying research note for completeness.

# Avoid GAP's pretty-printer inserting physical newlines into CSV rows or
# wrapping a long comment into executable text in generated .g files.
SR2_PrintTo := function(arg)
    local stream;
    stream:=OutputTextFile(arg[1],false);
    SetPrintFormattingStatus(stream,false);
    CallFuncList(PrintTo,Concatenation([stream],arg{[2..Length(arg)]}));
    CloseStream(stream);
end;
SR2_AppendTo := function(arg)
    local stream;
    stream:=OutputTextFile(arg[1],true);
    SetPrintFormattingStatus(stream,false);
    CallFuncList(PrintTo,Concatenation([stream],arg{[2..Length(arg)]}));
    CloseStream(stream);
end;

SR2_NormalizeB := function(B)
    local F, C, m, M, i;
    F := GF(2);
    if Length(B)=0 then Error("B must contain at least one alternating form"); fi;
    C := List(B,M -> List(M,row -> List(row,x -> x*One(F))));
    m := Length(C[1]);
    for M in C do
        if Length(M)<>m or not ForAll(M,row -> Length(row)=m) then
            Error("all forms must have the same square dimension");
        fi;
        if M<>TransposedMat(M) or not ForAll([1..m],i -> IsZero(M[i][i])) then
            Error("forms must be alternating over GF(2)");
        fi;
    od;
    if RankMat(List(C,Flat))<>Length(C) then
        Error("forms must be linearly independent (W is the derived subgroup)");
    fi;
    return C;
end;

SR2_Bits := function(n,d)
    return List([0..d-1],i -> (QuoInt(n,2^i) mod 2)*One(GF(2)));
end;

SR2_PolarData := function(input)
    local B,m,r,lambdas,forms,ranks,radicals,k,a,M,F;
    B:=SR2_NormalizeB(input); m:=Length(B[1]); r:=Length(B); F:=GF(2);
    lambdas:=List([0..2^r-1],k -> SR2_Bits(k,r));
    forms:=[]; ranks:=[]; radicals:=[];
    for k in [1..2^r] do
        M:=NullMat(m,m,F);
        for a in [1..r] do M:=M+lambdas[k][a]*B[a]; od;
        Add(forms,M); Add(ranks,RankMat(M)); Add(radicals,NullspaceMat(M));
    od;
    return rec(B:=B,m:=m,r:=r,lambdas:=lambdas,forms:=forms,
        ranks:=ranks,radicals:=radicals,
        stem:=RankMat(Concatenation(B))=m);
end;

SR2_MFPolar := function(P)
    local i,j,k,codim,e;
    for i in [2..Length(P.lambdas)] do
        for j in [i+1..Length(P.lambdas)] do
            k:=Position(P.lambdas,P.lambdas[i]+P.lambdas[j]);
            codim:=RankMat(Concatenation(P.forms[i],P.forms[j]));
            e:=(P.ranks[i]+P.ranks[j]+P.ranks[k])/2-codim;
            if e<>0 then
                return rec(ok:=false,lambda:=P.lambdas[i],mu:=P.lambdas[j],
                    multiplicity:=2^e,exponent:=e);
            fi;
        od;
    od;
    return rec(ok:=true);
end;

SR2_Q0Scalar := function(M,v)
    local s,i,j;
    s:=Zero(GF(2));
    for i in [1..Length(v)] do
        for j in [i+1..Length(v)] do s:=s+v[i]*v[j]*M[i][j]; od;
    od;
    return s;
end;

# All ambivalent Q form an affine solution space: particular + span(kernel).
# Returns fail if there is no solution. Flat coordinates are row-major in Q.
SR2_RealityAffine := function(P)
    local A,b,k,v,row,i,a,solution,kernel,n;
    A:=[]; b:=[]; n:=P.m*P.r;
    for k in [2..Length(P.lambdas)] do
        for v in P.radicals[k] do
            row:=[];
            for i in [1..P.m] do
                for a in [1..P.r] do Add(row,v[i]*P.lambdas[k][a]); od;
            od;
            Add(A,row); Add(b,SR2_Q0Scalar(P.forms[k],v));
        od;
    od;
    if Length(A)=0 then
        solution:=List([1..n],i -> Zero(GF(2)));
        kernel:=IdentityMat(n,GF(2));
    else
        solution:=SolutionMat(TransposedMat(A),b);
        if solution=fail then return fail; fi;
        kernel:=NullspaceMat(TransposedMat(A));
    fi;
    return rec(particular:=solution,kernel:=kernel,dimension:=Length(kernel),
        number:=2^Length(kernel),equations:=A,rhs:=b,m:=P.m,r:=P.r);
end;

SR2_IsAmbivalentData := function(P,Q)
    local k,v;
    if Length(Q)<>P.m or not ForAll(Q,row -> Length(row)=P.r) then
        Error("Q must be m by r");
    fi;
    Q:=List(Q,row -> List(row,x -> x*One(GF(2))));
    for k in [2..Length(P.lambdas)] do
        for v in P.radicals[k] do
            if not IsZero((v*Q)*P.lambdas[k]+SR2_Q0Scalar(P.forms[k],v)) then
                return false;
            fi;
        od;
    od;
    return true;
end;

SR2_IsSRData := function(B,Q)
    local P;
    P:=SR2_PolarData(B);
    return SR2_MFPolar(P).ok and SR2_IsAmbivalentData(P,Q);
end;

SR2_Unflatten := function(x,m,r)
    return List([1..m],i -> x{[(i-1)*r+1..i*r]});
end;

SR2_QuadraticValue := function(B,Q,v)
    local value,a;
    value:=v*Q;
    for a in [1..Length(B)] do
        value[a]:=value[a]+SR2_Q0Scalar(B[a],v);
    od;
    return value;
end;

# callback(B,Q) is called once for every SR refinement of this fixed B.
# This is a streaming enumeration, not a list of GL-orbit representatives.
SR2_ForEachRefinement := function(B,callback)
    local P,A,n,c,x,i;
    P:=SR2_PolarData(B);
    if not SR2_MFPolar(P).ok then return 0; fi;
    A:=SR2_RealityAffine(P);
    if A=fail then return 0; fi;
    for n in [0..A.number-1] do
        c:=SR2_Bits(n,A.dimension); x:=ShallowCopy(A.particular);
        for i in [1..A.dimension] do x:=x+c[i]*A.kernel[i]; od;
        callback(P.B,SR2_Unflatten(x,P.m,P.r));
    od;
    return A.number;
end;

# Construct a consistent power-conjugate presentation on x_1,...,x_m,z_1,...,z_r.
SR2_Group := function(input,Q)
    local B,m,r,F,g,c,i,j,a,w;
    B:=SR2_NormalizeB(input); m:=Length(B[1]); r:=Length(B);
    Q:=List(Q,row -> List(row,x -> x*One(GF(2))));
    if Length(Q)<>m or not ForAll(Q,row -> Length(row)=r) then Error("invalid Q"); fi;
    F:=FreeGroup(m+r); g:=GeneratorsOfGroup(F);
    c:=SingleCollector(F,List([1..m+r],i -> 2));
    for i in [1..m] do
        w:=One(F);
        for a in [1..r] do if not IsZero(Q[i][a]) then w:=w*g[m+a]; fi; od;
        SetPower(c,i,w);
        for j in [i+1..m] do
            w:=g[j];
            for a in [1..r] do if not IsZero(B[a][i][j]) then w:=w*g[m+a]; fi; od;
            SetConjugate(c,j,i,w);
        od;
    od;
    UpdatePolycyclicCollector(c);
    return GroupByRws(c);
end;

# Extract B,Q using W=G'. Returns fail when the required exponent/class
# hypotheses do not hold. Abelian groups are deliberately handled separately.
SR2_Extract := function(G)
    local D,H,pcD,pcH,epi,x,m,r,B,Q,i,j,a,v,F;
    D:=DerivedSubgroup(G);
    if Size(D)=1 or not IsElementaryAbelian(D) or not IsSubgroup(Centre(G),D) then
        return fail;
    fi;
    epi:=NaturalHomomorphismByNormalSubgroup(G,D); H:=Image(epi);
    if not IsElementaryAbelian(H) then return fail; fi;
    pcD:=Pcgs(D); pcH:=Pcgs(H); m:=Length(pcH); r:=Length(pcD); F:=GF(2);
    x:=List(pcH,h -> PreImagesRepresentative(epi,h));
    B:=List([1..r],a -> NullMat(m,m,F));
    Q:=List(x,g -> List(ExponentsOfPcElement(pcD,g^2),a -> a*One(F)));
    for i in [1..m] do
        for j in [i+1..m] do
            v:=ExponentsOfPcElement(pcD,Comm(x[i],x[j]));
            for a in [1..r] do B[a][i][j]:=v[a]*One(F); B[a][j][i]:=v[a]*One(F); od;
        od;
    od;
    return rec(B:=B,Q:=Q,m:=m,r:=r,stem:=Size(Centre(G))=Size(D));
end;

# Independent definition-level validator using exact ordinary characters.
SR2_CharacterCheck := function(G)
    local T,irr,i,j,chi,product,mults;
    T:=CharacterTable(G); irr:=Irr(T);
    if not ForAll(irr,chi -> ForAll(ValuesOfClassFunction(chi),
        x -> x=ComplexConjugate(x))) then return false; fi;
    for i in [1..Length(irr)] do
        for j in [i..Length(irr)] do
            # Tensoring by a degree-one character permutes Irr(G).
            if irr[i][1]=1 or irr[j][1]=1 then continue; fi;
            product:=irr[i]*irr[j];
            mults:=List(irr,chi -> ScalarProduct(T,product,chi));
            if not ForAll(mults,n -> n=0 or n=1) then return false; fi;
        od;
    od;
    return true;
end;
