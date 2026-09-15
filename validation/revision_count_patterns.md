# 第一篇论文的计数函数：审查、证明与实算结果

已完成的英文正文片段为 `worknotes/revision_count_patterns.tex`。没有改动主稿或 BibTeX；没有假称读过不存在的 `SRgroup.tex`。实际读取了 `upload/OrderSRgp.tex`、`upload/SRgroup2000.tex` 的原始 patterns 段、当前主稿的结构结果与 `output/f_of_n.csv`。

## 建议纳入第一篇的精确结果

令 a(m) 为奇数 m 阶 abelian groups 的同构类型数，p(e) 为整数分拆数。

1. **全范围显式公式：** 对每个奇数 m≥1，
   f(2m)=a(m)=∏_{p^e∥m}p(e)。
   而且完整同构分类为 Dih(A)，A 跑遍 m 阶 abelian groups。
2. **全范围的 4m 阶结构分类：** 任意 4m 阶 SR 群都有唯一 normal odd Hall subgroup M，且 G/M≅V4。M abelian 当且仅当 G≅Dih(A)×Dih(B)，其 unordered pair {A,B} 唯一。
3. 因此 abelian Hall 部分的精确公式是
   h(m)=½[∑_{d|m}a(d)a(m/d)+1_{m square}a(√m)]。
   全范围恒等式 f(4m)=h(m)+e(m)，其中 e(m) 是 nonabelian odd Hall 部分的真实分类计数。**这不是将未知量 e(m) 消掉的全范围显式解。**
4. **本次完整有限范围公式：** 奇数 m、4m≤2000 时，e(m)=1 当且仅当 m∈{75,375,405}，其余为零。三个例外分别是 [300,25]、[1500,37]、[1620,422]。
5. **全阶精确递推：** c(n) 为无 C2 direct factor 的 SR 群数，则 c(n)=f(n)−f(n/2)，并且 f(2^k m)=∑_{j=0}^k c(2^j m)。这里 f(noninteger)=0，f(1)=1。由此 f(2n)≥f(n)。
6. **class-two 精确增量：** N2(2^k)−N2(2^(k−1))=s_k。依据中心分裂定理，无 C2 direct factor 正好是 special。8,…,1024 阶的 N2 为 2,2,7,10,20,42,95,221，增量 s_k 为 2,0,5,3,10,22,53,126。

第 1、2、3、5、6 项均有完整手证。第 4 项是以主稿的 finite census completeness 为输入的 computer-assisted corollary，有额外可重复记录审计，不被表格拟合替代。可以把前两条描述为论文证明的结果；未经更全面文献查新，不宜自行加上“首次”或“前人未知”。

## 证明中必须保留的细节

### 对原 2m 阶证明缺口的修补

OrderSRgp 从“逆元共轭须使用 M 外的元素”直接跳到“同一个 involution t 逐个取逆所有元素”，这是不成立的推理步骤：共轭者会随元素变化。

正文新证法如下。

- 任意 2m 阶群的 regular permutation sign homomorphism 非平凡：involution 左乘排列由 m 个 transpositions 组成。因此 kernel M 有奇数阶 m，且 M 恰是所有奇数阶元素的集合。这一步甚至不需要 Burnside transfer 或 SR。
- 若非平凡 x∈C_M(t)，SR 要求将 x 送到 x⁻¹ 的共轭者必须在 M 外，写成 yt。但 t centralizes x，这就会使 y∈M 本身将 x 送到 x⁻¹，与奇数阶群不含非平凡 real elements 矛盾。
- 因此 t 在 M 上 fixed-point-free。设 α=conjugation by t，则 x↦x⁻¹α(x) 为 injection 从而 surjection；α 把所有像取逆，所以在所有 M 上就是 inversion，因此 M abelian。
- 反向用 χλχμ=χλμ+χλμ⁻¹（χ1=1+ε）逐一排除碰撞，而不是仅说所有非线性次数为 2 就自动 multiplicity-free。

### 4m 阶 normal Hall 与两个 sign spaces

- Kazarin–Chankov solvability 定理保证非平凡 G/G′；它又是 elementary abelian 2，所以有正常 index-two H。
- |H|=2m，前述 sign-kernel 论证给 characteristic odd Hall M⊂H；于是 M normal G。
- G/M 是 4 阶 SR 群，故 V4；Sylow 2-subgroup 就是一个 complement。
- M abelian 时，V4 在奇数阶 M 上可按四个 ±1 characters 作同时分解。trivial component central，故为零。
- 三个 nontrivial sign components 若同时非零，选择每个中非零元素之和，要取逆就需要三个 signs 同时 −1；而三个 nontrivial characters 之积恒为 +1。这与 ambivalence 矛盾。
- 至多两个 sign components，选 dual basis 即得到两个 generalized dihedral 因子。
- 每个 Dih(A) 的阶只被 2 整除一次，不能分解为两个非平凡 SR 直因子（两者都须偶数阶）。Krull–Remak–Schmidt 给 unordered pair 的唯一性，Burnside 的二元交换 orbit formula 得 h(m)。

## 实际执行的验证

`worknotes/verify_count_formulas.py output/f_of_n.csv --output worknotes/count_formulas_certificate.json`：

- 检查全部 500 个 n≡2 (mod 4) 条目，全部等于 a(n/2)。
- 检查全部 250 个 n≡4 (mod 8) 条目，仅三条比 h(n/4) 大，且各大 1。
- 此 progression 总共 654 个记录，h(m) 加总为 651。
- 检查当前表中全部可比较的 f(2n)≥f(n)，全部通过。

`worknotes/audit_count_patterns.g` 实际读取上述 progression 对应的 654 个 SmallGroups identifier，逐个构造 GAP 群、计算 Hall subgroup、核验 normality 和 V4 quotient。651 个 Hall subgroup abelian，恰有下列三例 nonabelian：

| Identifier | Hall order | GAP 结构显示 | 独立 ordinary-character SR test |
|---|---:|---|---|
| [300,25] | 75 | (C5×C5)⋊C3 | true |
| [1500,37] | 375 | ((C5×C5)⋊C5)⋊C3 | true |
| [1620,422] | 405 | C3^4⋊C5 | true |

最后还实际检查 A5 的 SR test=false，以及 D16 的 exponent=8、SR test=true。

干净最终运行日志为 `worknotes/audit_count_patterns.log`，末尾 `COMPLETE 654 records`，`PROCESS_EXIT=0`，实际约 1.62 秒。最初运行发现本地 GAP 基础文件 clasmax.grp 截断导致启动语法错误；负责 runtime 的同事修复后已重跑，最终日志干净。初始错误日志只作为内部诊断保留，不是最终计算证据。

上述程序审计的是 supplied complete census 的记录性质；它本身不重新扫描此 progression 中所有非 SR 的 SmallGroups。正文明确将其作为主普查上的结构审计，不能说这是新的独立 exhaustive enumeration。

程序输入 `worknotes/count_v2_2_ids.g` 从 `upload/SR_groups_results2_2000.csv` 精确解析取得。打包时应保留三个文件或将该输入内嵌进 GAP 文件；当前 GAP audit 从项目根目录运行，读取 `worknotes/count_v2_2_ids.g` 与 `output/sr_review/sr_class2.g`。

## OrderSRgp.tex 中不采用的错误论断

- “f(n) 严格由 2-adic valuation 决定”错误：例如 f(6)=1、f(18)=2，但 v2 都是 1。正确关系同时依赖奇部的 prime-exponent partitions。
- “A5 是 60 阶第二个 SR 群”错误。A5 的次数 5 的不可约特征标满足 χ5²=1+χ3+χ3′+2χ4+2χ5。实际两类为 C2×Dih(C15) 与 Dih(C3)×Dih(C5)=S3×D10。
- “4m 阶 SR 群会因为 simple sections 而无正常 2-complement”错误。SR 全部 solvable；本节上述定理更证明这些群全部都有正常 odd Hall complement。真正的例外是 odd Hall **非交换**，不是没有 normal complement。
- “SR 2-groups exponent≤4”错误：D16 已有 exponent 8 且 SR，一般 dihedral 2-groups 给任意大 2-power exponent。quadratic 模型的 exponent≤4 限制来自本文专门的 **class-two** 情形，不得套到整个 2-group 类。
- “nilpotency 自动支持 Wigner 条件”“SR 必须由 involutions 生成”等说法没有有效论证；Q8 是 SR 却不由 involutions 生成。
- 2-power counts 的末项 483 已过时；当前 census 总数为 803，并应保留正文已经说明的计算完整性依赖。
- “绝大多数 SR groups 位于 2-power orders”在当前 ≤2000 的按群计数中也不对：2-power 部分为 1351/7889，约 17.1%。若只说图中某些 2-power orders 有显著峰值则可以，但它不是分类定理。
- 原文中的 PORC、submanifold、bounded derived length 等语句没有建立明确计数结论，不宜放在这篇 classification 论文末节。

## 文献核验与后续界限

本文新增证明仅用现成引用 KC10 的 solvability；原始来源页面已核实：
https://www.mathnet.ru/eng/sm7540
其英文题名为 *Finite simply reducible groups are soluble*，Sb. Math. 201:5 (2010), 655–668，DOI 10.1070/SM2010v201n05ABEH004087。当前主稿 BibTeX 使用俄文刊名页码配英文 DOI，主编可酌情统一英文版本信息。

Burnside normal complement 定理也在 David Craven 的原作者讲义核对过，但正文现已通过 regular sign 给出自含证明，不需要为这一步增加引用：
https://web.mat.bham.ac.uk/D.A.Craven/docs/lectures/finitegroups2010.pdf

现有 norm-circle family 的全 q 手证能够解释 [300,25] 与 [1620,422] 的来源，并表明超出 2000 后仍存在额外 e(m)>0 的例子，但把整段长 family proof 重新加入最后一节会分散 classification 主线。建议只保留主稿所需结构审计或简短联系，不以 three exceptions 推断无限范围例外列表。完整全偶阶 f(n) 显式公式没有在本次工作中得到，不应虚构一个。
