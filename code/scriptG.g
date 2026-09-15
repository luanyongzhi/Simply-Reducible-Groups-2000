#############################################################################
##
##  scriptG.g  --  self-contained Stage 3
##
##  Enumerates all simply reducible (SR) groups of order 1024 with
##  2-class >= 3, as immediate ANUPQ descendants (step s) of SR parents
##  of order 2^(10-s).  Completeness by heredity: SR is quotient-closed,
##  and parent(G) = G/(last lower exponent-2 central term) is a quotient.
##
##  EXCLUDED by design: class-2 groups of order 1024 (descendants of
##  elementary abelian parents) -- Stage 4 handles those.
##
##  Usage:   gap -o 2300g -q scriptG.g        (adjust -o to your RAM)
##  Input:   SR_groups_results512.csv  (CSV; first two fields are Order,SmallGroup_ID)
##  Output:  scriptG_results.log  -- checkpoint + results; rerun-safe.
##
#############################################################################

if LoadPackage("anupq") <> true then
  Error("ANUPQ (with compiled pq binary) is required");
fi;
SetInfoLevel(InfoANUPQ, 0);
## ---------------------------------------------------------------------
## Chunked replacement for ANUPQ's PqList: statement-by-statement Read()
## instead of whole-file-as-one-function parse (dies at ~4.7GB files).
## Format verified against job [128,2163]: entries do "Add( L, F )".
## ---------------------------------------------------------------------
ANUPQmagic := fail;;  ANUPQautos := [];;  ANUPQgroups := [];;
MakeReadWriteGlobal( "PqList" );
PqList := function( file )
  local L, out, i, res;
  ANUPQmagic  := fail;
  ANUPQautos  := [];
  ANUPQgroups := [];
  if not IsReadableFile( file ) then return fail; fi;
  Read( file );                    # per-statement parse: no size limit
  if ANUPQmagic = fail then return fail; fi;
  ANUPQautos := [];                # autos never used by Script G; free now
  L := [];  out := [];
  for i in [ 1 .. Length( ANUPQgroups ) ] do
    if IsBound( ANUPQgroups[i] ) then
      res := CALL_WITH_CATCH( ANUPQgroups[i], [ L ] );
      if res[1] <> true or Length( L ) = 0 then
        Error( "chunked PqList: bad entry ", i, " in ", file );
      fi;
      Add( out, L[ Length( L ) ] );
      L := [];                     # reset appender (entries use Add(L, F))
      Unbind( ANUPQgroups[i] );    # free the parsed function as we go
      if Length( out ) mod 50000 = 0 then
        Print( "    [chunked read] ", Length( out ), " descendants built\n" );
      fi;
    fi;
  od;
  ANUPQgroups := [];
  return out;
end;;

SR512_CSV    := "SR_groups_results512.csv";;
if not IsExistingFile(SR512_CSV) then
  SR512_CSV := "../data/SR_groups_results512.csv";
fi;
RESULTFILE   := "scriptG_results.log";;
DO_STEP1     := true;;    # include s=1 jobs (order-512 parents; Stage-2 rerun)
VERIFY_TESTS := true;;    # self-audit the SR test before the real runs
PQWORKSPACE  := 10^8;;    # pq workspace (words); raise and rerun if jobs FAIL

## ------------------------------------------------- crash-safe line logger
LogG := function(line)
  local o;
  o := OutputTextFile(RESULTFILE, true);   # append; open/close per line
  SetPrintFormattingStatus(o, false);      # NO line-wrapping of big ints
  WriteLine(o, line);
  CloseStream(o);
end;;

## ------------------------------------------------- the SR test (character-free)
IsSimplyReducibleFast := function(G)
  local n, els, cnt, x, p, s2, s3, cc, k;
  n := Size(G);
  # necessary: G^ab elementary abelian (ambivalence forces exponent 2)
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
  if s2 <> n * k then return false; fi;                    # not ambivalent
  s3 := Sum(cnt, v -> v^3);
  return s3 = Sum(cc, c -> Size(c) * (n / Size(c))^2);     # mult-free test
end;;

## ------------------------------------------------- slow cross-check version
IsSimplyReducibleCT := function(G)
  local t, irr, m, i, j, chi;
  t := CharacterTable(G);  irr := Irr(t);  m := Length(irr);
  if not ForAll(irr, x -> x = ComplexConjugate(x)) then return false; fi;
  for i in [1..m] do for j in [i..m] do for chi in irr do
    if ScalarProduct(t, irr[i]*irr[j], chi) > 1 then return false; fi;
  od; od; od;
  return true;
end;;

## ------------------------------------------------- cheap direct-factor flag
HasC2DirectFactor := function(G)
  return not IsSubset(FrattiniSubgroup(G), Omega(Centre(G), 2));
end;;

## ------------------------------------------------- self-tests
if VERIFY_TESTS then
  Print("[selftest] fast vs character-table test, all orders 2..64 ...\n");
  for k in [1..6] do
    for i in [1..NumberSmallGroups(2^k)] do
      G := SmallGroup(2^k, i);
      if IsSimplyReducibleFast(G) <> IsSimplyReducibleCT(G) then
        Error("SR test mismatch at [", 2^k, ",", i, "]");
      fi;
    od;
  od;
  if not IsSimplyReducibleFast(SmallGroup(128,928)) then Error("anchor 928"); fi;
  if not IsSimplyReducibleFast(SmallGroup(128,937)) then Error("anchor 937"); fi;
  if IsSimplyReducibleFast(SmallGroup(128,931))     then Error("anchor 931"); fi;
  if IsSimplyReducibleFast(SmallGroup(16,8))        then Error("anchor SD16"); fi;
  Print("[selftest] passed\n");
fi;

## ------------------------------------------------- SR ids, orders 2..256
SR128_REFERENCE := [161,163,444,445,928,937,978,979,984,991,993,1578,1599,
  1600,1740,1743,1746,1749,1755,1876,1877,1891,2011,2018,2024,2025,2110,
  2114,2140,2142,2147,2149,2163,2172,2173,2194,2198,2209,2216,2218,2230,
  2231,2235,2306,2308,2313,2315,2320,2321,2323,2324,2326,2327,2328];;

SRids := [];;
for k in [1..8] do
  SRids[k] := [];
  for i in [1..NumberSmallGroups(2^k)] do
    if IsSimplyReducibleFast(SmallGroup(2^k, i)) then Add(SRids[k], i); fi;
    if i mod 5000 = 0 then
      Print("    order ", 2^k, ": ", i, " scanned, ",
            Length(SRids[k]), " SR so far\n");
    fi;
  od;
  Print("[bootstrap] order ", 2^k, ": ", Length(SRids[k]), " SR groups\n");
od;
if SRids[7] <> SR128_REFERENCE then
  Error("recomputed order-128 SR list disagrees with the certified CSV!");
fi;
Print("[bootstrap] order-128 list matches certified 54  -- audit passed\n");

## ------------------------------------------------- order-512 list from CSV
LoadSR512 := function(fname)
  local s, l, f, ids;
  s := StringFile(fname);
  if s = fail then
    Error("cannot read ", fname,
          " -- export your certified order-512 SR list to this file");
  fi;
  ids := [];
  for l in SplitString(s, "\n") do
    f := Filtered(SplitString(l, ",\t "), x -> x <> "");
    if Length(f) >= 2 and Int(f[1]) = 512 and Int(f[2]) <> fail then
      AddSet(ids, Int(f[2]));
    fi;
  od;
  return ids;
end;;

if DO_STEP1 then
  SRids[9] := LoadSR512(SR512_CSV);
  Print("[input] order 512: ", Length(SRids[9]),
        " SR ids loaded (expected: 317)\n");
fi;

## ------------------------------------------------- job list
JOBS := [];;
for s in [7,6,5,4,3,2] do                    # Script G proper first ...
  for id in SRids[10-s] do
    if not IsElementaryAbelian(SmallGroup(2^(10-s), id)) then
      Add(JOBS, [10-s, id, s]);
    fi;
  od;
od;
Print("[jobs] s >= 2 jobs: ", Length(JOBS), "  (earlier run: 223)\n");
if DO_STEP1 then                             # ... then the s=1 stratum
  for id in SRids[9] do
    if not IsElementaryAbelian(SmallGroup(512, id)) then
      Add(JOBS, [9, id, 1]);
    fi;
  od;
fi;
Print("[jobs] total jobs: ", Length(JOBS), "\n");

## ------------------------------------------------- checkpoint reader
LoadCheckpoint := function(fname)
  local s, l, f, done, hits;
  done := [];  hits := [];
  s := StringFile(fname);
  if s = fail then return rec(done := [], hits := []); fi;
  for l in SplitString(s, "\n") do
    f := Filtered(SplitString(l, " "), x -> x <> "");
    if Length(f) >= 4 and f[1] = "DONE" then
      AddSet(done, JoinStringsWithSeparator(f{[2..4]}, "_"));
    elif Length(f) >= 5 and f[1] = "SRHIT" then
      AddSet(hits, List(f{[2..5]}, Int));    # [k, parent-id, s, pc-code]
    fi;
  od;
  return rec(done := done, hits := hits);
end;;

## ------------------------------------------------- one job, fully isolated
RunJob := function(k, id, s)
  local P, res, desc, D, nf, code;
  P := SmallGroup(2^k, id);
  res := CALL_WITH_CATCH(function()
           return PqDescendants(P : StepSize := s, Workspace := PQWORKSPACE);
         end, []);
  if res[1] <> true then                     # fallback: standard presentation
    res := CALL_WITH_CATCH(function()
             local S;
             S := PcGroupFpGroup(StandardPresentation(P));
             return PqDescendants(S : StepSize := s, Workspace := PQWORKSPACE);
           end, []);
  fi;
  if res[1] <> true then return fail; fi;
  desc := res[2];  nf := 0;
  for D in desc do
    if Size(D) = 1024 and IsSimplyReducibleFast(D) then
      nf := nf + 1;
      code := CodePcGroup(D);
      LogG(Concatenation("SRHIT ", String(k), " ", String(id), " ",
                        String(s), " ", String(code)));
    fi;
  od;
  LogG(Concatenation("DONE ", String(k), " ", String(id), " ", String(s),
                    " ", String(nf), " ", String(Length(desc))));
  return rec(nfound := nf, ndesc := Length(desc));
end;;

## ------------------------------------------------- main loop (restart-safe)
ck := LoadCheckpoint(RESULTFILE);;
FAILED := [];;
for job in JOBS do
  key := Concatenation(String(job[1]), "_", String(job[2]), "_", String(job[3]));
  if key in ck.done then continue; fi;
  Print("JOB [", 2^job[1], ",", job[2], "] step ", job[3], " ... \c");
  r := RunJob(job[1], job[2], job[3]);
  if r = fail then
    Add(FAILED, job);
    LogG(Concatenation("FAILED ", String(job[1]), " ", String(job[2]),
                      " ", String(job[3])));
    Print("FAILED\n");
  else
    Print(r.nfound, " SR of ", r.ndesc, " descendants\n");
  fi;
  GASMAN("collect");
  Exec(Concatenation("rm -rf ", Filename(ANUPQData.tmpdir, ""),
                     "* 2>/dev/null"));
od;

## ------------------------------------------------- summary + audits
final := LoadCheckpoint(RESULTFILE);;
Print("\n================= SUMMARY =================\n");
for s in [1..7] do
  Print("s = ", s, " : ", Number(final.hits, h -> h[3] = s), " SR groups\n");
od;
n2 := Number(final.hits, h -> h[3] >= 2);;
Print("s>=2 total (Script G proper)        : ", n2, "\n");
Print("s=1  total (Stage-2 rerun; was 474) : ",
      Number(final.hits, h -> h[3] = 1), "\n");
Print("ALL class>=3 SR groups of order 1024: ", Length(final.hits), "\n");
Print("failed jobs: ", FAILED, "\n");

# Audit: D(C8^3), Q(C8^3) in the s=3 stratum over parents [128,1599/1600]
CountInv := G -> Number(AsSSortedList(G), x -> Order(x) = 2);;
hits1599 := Filtered(final.hits, h -> h[1] = 7 and h[2] = 1599 and h[3] = 3);;
Print("Audit s=3 over [128,1599]: involution counts = ",
      SortedList(List(hits1599, h ->
        Number(AsSSortedList(PcGroupCode(h[4], 1024)), x -> Order(x) = 2))),
      "  (predict [ 7, 519 ] = Q(C8^3), D(C8^3);",
      " Q descends from the dihedral seed, not from 1600)\n");
Print("Decomposability split of hits (Omega_1(Z) vs Frattini):\n");
Print("  with C2 direct factor: ", Number(final.hits,
      h -> HasC2DirectFactor(PcGroupCode(h[4], 1024))), " / ",
      Length(final.hits), "\n");

Print("=================  END  ===================\n");
