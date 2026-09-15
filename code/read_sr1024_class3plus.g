# Read exact presentation codes from the published class >= 3 CSV.
# Usage from the reproducibility root:
#   Read("code/read_sr1024_class3plus.g");;
#   rows := SR1024ReadClass3Plus("data/sr1024_class3plus.csv");;
#   G := SR1024GroupFromClass3Row(rows[1]);;
# This reader does not run ANUPQ or require a SmallGroup ID at order 1024.
# Code is a presentation integer, not a canonical isomorphism identifier.

SR1024ReadClass3Plus := function(path)
  local text, lines, header, expected, rows, line, fields, values, i;
  text := StringFile(path);
  if text = fail then Error("Cannot read CSV: ", path); fi;
  lines := Filtered(SplitString(text,"\n",""),
                   x -> Length(NormalizedWhitespace(x)) > 0);
  if Length(lines) < 1 then Error("Empty CSV: ",path); fi;
  expected := ["No","ParentOrder","ParentId","Step","Class","Exponent",
               "CenterSize","DerivedLength","NrConjClasses",
               "HasC2DirectFactor","Code"];
  header := List(SplitString(lines[1],",",""),NormalizedWhitespace);
  if header <> expected then Error("Unexpected class >= 3 CSV schema"); fi;
  rows := [];
  for i in [2..Length(lines)] do
    line := lines[i];
    fields := List(SplitString(line,",",""),NormalizedWhitespace);
    if Length(fields) <> 11 then Error("Invalid CSV row: ",i); fi;
    values := List(fields{[1..9]},Int);
    if fail in values or Int(fields[11]) = fail then
      Error("Non-integer CSV field on row: ",i);
    fi;
    if not fields[10] in ["true","false"] then
      Error("Invalid direct-factor flag on row: ",i);
    fi;
    if values[1] <> Length(rows)+1 or values[5] < 3
       or Int(fields[11]) <= 0 then Error("Invalid row metadata: ",i); fi;
    Add(rows,rec(number:=values[1],parentOrder:=values[2],
      parentId:=values[3],step:=values[4],class:=values[5],
      exponent:=values[6],centerSize:=values[7],derivedLength:=values[8],
      nrConjClasses:=values[9],hasC2DirectFactor:=fields[10]="true",
      code:=Int(fields[11])));
  od;
  return rows;
end;;

SR1024GroupFromClass3Row := function(row)
  local G;
  G := PcGroupCode(row.code,1024);
  if Size(G) <> 1024 or NilpotencyClassOfGroup(G) <> row.class then
    Error("Reconstruction disagrees with order/class metadata at No=",row.number);
  fi;
  return G;
end;;

# Optional exact SR validation for a reconstructed group: Wigner moments.
SR1024CheckSR := function(G)
  local classes, counts, n;
  n := Size(G);
  classes := ConjugacyClasses(G);
  counts := List(Collected(List(Elements(G),x -> x^2)),x -> x[2]);
  return Sum(counts,x -> x^2) = n*Length(classes)
     and Sum(counts,x -> x^3) = Sum(classes,C -> n^2/Size(C));
end;;
