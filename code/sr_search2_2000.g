# 1. 极其高效的纯元素层面单约化判定算法
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


# 2. 带有消除自动换行、清洗反斜杠功能的写出算法
SearchAllSRGroups := function(start_order, end_order, filename)
    local n, num_groups, i, G, out_str, desc;

    # 强制将 GAP 虚拟屏幕调到宽屏模式，防止长结构描述被截断
    SizeScreen([4096, 4096]);

    AppendTo(filename, "Order, SmallGroup_ID, StructureDescription, CenterSize\n");

    for n in [start_order .. end_order] do
        if n mod 2 <> 0 then
            continue;
        fi;

        if n = 512 or n = 1024 or n = 1536 then
            Print("WARNING: Skipping order ", n, " due to extreme group counts. \n");
            continue;
        fi;

        num_groups := NrSmallGroups(n);
        Print("Processing Order ", n, " (Total Groups: ", num_groups, ")...\n");

        for i in [1 .. num_groups] do
            G := SmallGroup(n, i);

            if IsSimplyReducibleFast(G) then

                # 提取结构名字并彻底清洗残留的格式符
                desc := StructureDescription(G);
                desc := ReplacedString(desc, "\n", "");
                desc := ReplacedString(desc, "\\", "");

                # 【新增防切割绝杀】：把群名字里的逗号替换成分号，彻底断绝 Excel 误判！
                desc := ReplacedString(desc, ",", ";");

                # 【严格 CSV 格式】：去掉逗号周围的所有空格
                out_str := Concatenation(
                    String(n), ",",
                    String(i), ",\"",
                    desc, "\",",
                    String(Size(Center(G))), "\n"
                );
                AppendTo(filename, out_str);
            fi;
        od; # 结束 i 循环
    od; # <--- 就是这里！您刚才漏掉了结束 n 循环的 od;

    Print("\nSearch Completed! Results saved to ", filename, "\n");
end;

# 3. 触发运行
SearchAllSRGroups(2, 2000, "SR_groups_results2_2000.csv");
