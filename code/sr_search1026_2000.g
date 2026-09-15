IsSimplyReducibleFast := function(G)
    local classes, c, rep, sum_nu2, els, sqs, sum_xi3, len, i, count;

    classes := ConjugacyClasses(G);


    for c in classes do
        rep := Representative(c);
        if not rep^-1 in c then
            return false;
        fi;
    od;


    sum_nu2 := 0;
    for c in classes do
        rep := Representative(c);
        sum_nu2 := sum_nu2 + Size(c) * (Size(Centralizer(G, rep))^2);
    od;


    els := AsList(G);
    sqs := List(els, x -> x^2);


    Sort(sqs);

    sum_xi3 := 0;
    len := Length(sqs);
    i := 1;
    while i <= len do
        count := 1;
        while i + count <= len and sqs[i+count] = sqs[i] do
            count := count + 1;
        od;
        sum_xi3 := sum_xi3 + count^3;
        i := i + count;
    od;


    return sum_xi3 = sum_nu2;
end;



SearchAllSRGroups := function(start_order, end_order, filename)
    local n, num_groups, i, G, out_str;


    AppendTo(filename, "Order, SmallGroup_ID, StructureDescription, CenterSize\n");


    for n in [start_order .. end_order] do
        if n mod 2 <> 0 then
            continue;
        fi;


        if n = 512 or n = 1024 then
            Print("WARNING: Skipping order ", n, " due to extreme group counts. Please run these via ParGAP.\n");
            continue;
        fi;

        num_groups := NrSmallGroups(n);
        Print("Processing Order ", n, " (Total Groups: ", num_groups, ")...\n");

        for i in [1 .. num_groups] do
            G := SmallGroup(n, i);


            if IsSimplyReducibleFast(G) then

                out_str := Concatenation(
                    String(n), ", ",
                    String(i), ", ",
                    "\"", StructureDescription(G), "\", ",
                    String(Size(Center(G))), "\n"
                );
                AppendTo(filename, out_str);
            fi;
        od;
    od;

    Print("\nSearch Completed! Results saved to ", filename, "\n");
end;

SearchAllSRGroups(1026,2000, "SR_groups_results1026_2000.csv");
