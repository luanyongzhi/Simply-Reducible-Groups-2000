# 新增普查源程序与最终 checkpoint 审计（2026-09-08）

本次仅审计原始程序、最终日志和已有CSV，不重跑539个ANUPQ任务，也不重跑512或1536的完整SmallGroups遍历。可复核脚本为 `worknotes/audit_census_checkpoint.py`，运行通过；机器可读结果为 `worknotes/revision_census_audit.json`。

## 结论及稿件措辞

新上传的原始checkpoint弥补了旧稿唯一明确的1024类≥3证据缺口：539个预期任务全部DONE，0 FAILED，581个正记录与现有CSV完全一致。可将正文的581和由此产生的803结论写为 **computer-assisted classification result**，不必继续保留“没有完整日志，故只能conditional”的旧限定。

数学上的完整性来自三环：完整的小阶SR输入；SR商封闭保证canonical-parent覆盖；ANUPQ的allowable-subgroup全自同构轨道代表恰对应一次每个immediate descendant同构类型。实际执行证据为完整且逐任务闭合的checkpoint。不同canonical parent的后代不能同构，因此无需额外581两两同构检验来证明ANUPQ输出的非重复性。不要将“581个不同PC编码”当成非同构证明。

这属于通常意义的计算机辅助分类，仍使用原始SmallGroups普查及ANUPQ的正确实现。不能改写成“本次独立重跑了全部普查/539任务”。512和1536新上传的是遍历源程序及既有正结果，不是每个被拒绝ID的独立复算日志；用户明确报告其原始大规模计算在七台机器上执行三周，正文可作为作者的计算记录陈述。

## 实际检查结果

- `sr512.csv` 与 `SR_groups_results512.csv` 不仅317行逐条一致，文件SHA256也完全相同。
- 从已交付8、16、32、64、128、256、512阶SR名单删除唯一初等交换群，构造预期 `(k,parentId,10-k)` 集合。539个键与checkpoint的DONE键完全相同，无遗漏、无额外键。
- `scriptG_results.FINAL.bak` 严格由539条DONE和581条SRHIT组成，没有FAILED，没有重复DONE，没有重复SRHIT。
- 每条DONE记录的接受数与对应任务的SRHIT个数一致；每个SRHIT均属于一个DONE任务。
- 581个 `(2^k,parentId,step,Code)` 与 `sr1024_class3plus.csv` 按顺序逐条一致。
- 127任务有SR输出，412任务输出零；两者均有完整DONE证据。
- 所有任务累计返回3,167,745个后代。按step分布如下。

| Parent order | Step | DONE jobs | Descendants | SR hits |
|---:|---:|---:|---:|---:|
|512|1|316|51,886|474|
|256|2|129|798,367|105|
|128|3|53|2,317,492|2|
|64|4|24|0|0|
|32|5|11|0|0|
|16|6|4|0|0|
|8|7|2|0|0|
|Total||539|3,167,745|581|

- 控制台显示selftest通过：所有2、4、8、16、32、64阶群，快速SR判定与ordinary character-table判定一致。注意这是原始日志的实测，不是本次重新运行。
- bootstrap计数依次为1、1、3、5、12、25、54、130；order128的54个ID另有固定名单比对。
- console保存的是最后恢复的一次运行，其正文主要列出316个step1任务；早期223个任务保存在完整checkpoint中。因此不能仅按console中JOB行数判定缺失。
- mainCSV6210、512CSV317、1536CSV559，总计7086个不同(order,ID)键。
- 本次GAP只读调用确认库中order512共有10,494,213群，1536共有408,641,062群。

## 原始源程序的筛选与分区

`sr_search2_2000.g` 跳过奇数及512、1024、1536；对其余每个even order的所有SmallGroup(n,i)遍历，先按每个共轭类检查inverse，再以square-root fibres计算三次矩并与中央化子表达式相等比较。此算法精确等价于文中SR判定。

`sr_search512.g` 和 `sr_search1536.g` 使用完全相同的筛选，分别遍历对应单个阶。输出有(order,ID,StructureDescription,CenterSize)；结构描述只帮助阅读，不参与同构判定。文中使用GAP不带版本号，软件版本保留在参考文献/代码档案。

`sr_search1026_2000.g` 是upper-range替代遍历，**并未排除1536**。它与mainCSV对应的上半区及独立1536任务均可能重叠。正文不可称四个程序天然形成不相交的partition；合并时必须以(order,ID)键去重。最终统计的三份CSV已经无重复，不存在本次发现的计数错误。

各搜索程序用AppendTo累积输出，重跑同一输出文件可重复附加header和正记录。原始数据目前无重复；投稿复现时从空文件开始或按精确(order,ID)键合并。

## 1024程序的数学范围

`scriptG.g` 调用PqDescendants(P : StepSize:=s)。官方默认ClassBound为PClass(P)+1，所以是immediate descendants；默认同时构造capable和terminal descendants，没有遗漏终端群。

这不是单纯central-C2 extension sweep。最后一个非平凡lower exponent-2 central term可有2^s阶，因此需要s=1至7。父群最小非交换阶为8，故k=3至9。

SR的交换化和中心均初等交换。类2 SR群满足P1=G'、P2=1，canonical parent初等交换。反之初等交换父群的immediate descendants有nilpotency class≤2。因此只跳过初等交换父群恰好排除待由quadratic-map方法处理的类2层，不损失class≥3 SR群。

ANUPQ在一个父群内通过自同构轨道给出每个同构型一次；canonical parent是同构不变量，跨父群不重叠。pc-code仅用于保存/重建，不负责去重。

`makeCSV.g` 逐个重建PC群、断言阶1024和class≥3、重新跑SR判定，并写出分类元数据。这解释了现有581CSV的产生方式。

## 可复现性上的最小修正

原始 `scriptG.g` 两处写 `Workspace := PQWORKSPACE`，但ANUPQ非交互式选项的文档名称是 `PqWorkspace`。这是资源参数拼写问题，不改变StepSize等数学筛选条件。原始文件保留不动；已写 `output/scriptG_reproduce.g`，仅更正两处拼写，并在头部说明用途及从空checkpoint运行。

原始checkpoint策略一般不保证任意中断都完全restart-safe：某任务已经写入若干SRHIT，但尚未写DONE便中断，重新运行该任务会再写其正记录。最终交付checkpoint没有任何重复，故不影响本次581结论。新复现若中断，应将未DONE任务原有临时SRHIT清除后重跑，或实现每任务事务输出；不能靠不相等PC编码做同构去重。

原始程序的chunked PqList只改变大文件读取/释放策略。正文无需把该工程措施提升为新的分类定理；档案保存原始代码即可。

## 已核对的官方primary sources

1. [SmallGrp manual, Chapter 1](https://gap-packages.github.io/smallgrp/doc/chap1.html)：≤2000除1024均有完整、无重复同构型代表列表；NumberSmallGroups计数可用不意味着SmallGroup对象存在。
2. [ANUPQ manual, Chapter 2](https://gap-packages.github.io/anupq/doc/chap2.html)：§2.1-4说明allowable-subgroup自同构轨道代表的quotients给出每个immediate-descendant同构型恰一次；§2.3说明p-group generation算法。
3. [ANUPQ manual, Chapter 4](https://gap-packages.github.io/anupq/doc/chap4.html)：§4.4-1为PqDescendants函数。
4. [ANUPQ manual, Chapter 6](https://gap-packages.github.io/anupq/doc/chap6.html)：默认ClassBound=PClass(parent)+1；StepSize是父子阶数之比的p对数；默认包含capable与terminal；PqWorkspace是正确资源参数名。

正文草稿使用已有bib键ANUPQ、SmallGrp、GAP即可，不需要为每页另建参考文献。`worknotes/revision_census.tex` 具有两section，分别对应普通阶普查和1024canonical-parent层；`worknotes/revision_census_algorithms.tex` 给出四个附录伪代码及相互引用。
