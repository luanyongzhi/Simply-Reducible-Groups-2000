## makeCSV.g -- rebuild all 581 SR groups of order 1024 from pcgs codes,
##              re-certify each, export sr1024_class3plus.csv

IsSimplyReducibleFast := function(G)
  local n, els, cnt, x, p, s2, s3, cc, k;
  n := Size(G);
  if not ForAll(AbelianInvariants(G), a -> a = 2) then return false; fi;
  els := AsSSortedList(G);
  cnt := ListWithIdenticalEntries(n, 0);
  for x in els do
    p := PositionSorted(els, x^2);
    cnt[p] := cnt[p] + 1;
  od;
  cc := ConjugacyClasses(G);
  k  := Length(cc);
  s2 := Sum(cnt, v -> v^2);
  if s2 <> n * k then return false; fi;
  s3 := Sum(cnt, v -> v^3);
  return s3 = Sum(cc, c -> Size(c) * (n / Size(c))^2);
end;;

HasC2DirectFactor := function(G)
  return not IsSubset(FrattiniSubgroup(G), Omega(Centre(G), 2));
end;;

s   := StringFile("scriptG_results.log");;
out := OutputTextFile("sr1024_class3plus.csv", false);;
SetPrintFormattingStatus(out, false);
AppendTo(out, "No,ParentOrder,ParentId,Step,Class,Exponent,CenterSize,",
              "DerivedLength,NrConjClasses,HasC2DirectFactor,Code\n");

n := 0;;  nC2 := 0;;  perstep := [0,0,0,0,0,0,0];;
for l in SplitString(s, "\n") do
  f := Filtered(SplitString(l, " "), x -> x <> "");
  if Length(f) = 5 and f[1] = "SRHIT" then
    k    := Int(f[2]);;  id := Int(f[3]);;
    st   := Int(f[4]);;  code := Int(f[5]);;
    G := PcGroupCode(code, 1024);;
    if Size(G) <> 1024 then Error("wrong order: ", l); fi;
    if NilpotencyClassOfGroup(G) < 3 then Error("class < 3: ", l); fi;
    if not IsSimplyReducibleFast(G) then Error("SR re-test FAILED: ", l); fi;
    n := n + 1;;  perstep[st] := perstep[st] + 1;;
    c2 := HasC2DirectFactor(G);;  if c2 then nC2 := nC2 + 1; fi;
    AppendTo(out, n, ",", 2^k, ",", id, ",", st, ",",
      NilpotencyClassOfGroup(G), ",", Exponent(G), ",",
      Size(Centre(G)), ",", DerivedLength(G), ",",
      NrConjugacyClasses(G), ",", c2, ",", code, "\n");
  fi;
od;
CloseStream(out);;

Print("rows written   : ", n, "  (expect 581)\n");
Print("per-step tally : ", perstep, "  (expect [ 474, 105, 2, 0, 0, 0, 0 ])\n");
Print("with C2 factor : ", nC2, "  (expect 221)\n");
QuitGap(0);
