# ksw2

记录了对双序列比对的学习，主要是参照原作者的代码风格，实现对从 NW 到 SW 的改写。 
学习对象：https://github.com/lh3/ksw2 

## 引用文件介绍

* `ksw2.h`：头文件，包含主要数据结构与接口声明。
* `ksw2_gg.c`：全局比对；Green 的标准 DP 形式。
* `ksw2_gg2.c`：全局比对；Suzuki 的对角线形式。
* `ksw2_gg2_sse.c`：全局比对；基于 SSE intrinsics 的 Suzuki 对角线向量化实现。

## 改写文件介绍

* `ksw2_sw.c`：基于 `ksw2_gg.c` 改写得到的标准 DP 版 SW（**！！！**已完成，但仍有以下问题：由于直接调用了全局比对的回溯函数，所以返回的 CIGAR 会呈现出从最优分数位置到（0，0）起点，目前尚未解决**！！！**）
* `ksw2_sw2.c`：基于 `ksw2_gg2.c` 改写得到的对角线版 SW（尚未完成）
* `ksw2_sw2_sse.c`：基于 `ksw2_gg2_sse.c` 改写得到的 SSE 向量化 SW（尚未完成）

（以下内容是自己的一些理解）

## 名词解释

* **匹配（match）**：两条序列在某一位匹配成功（字符相同），会得到加分奖励
* **错配（mismatch）**：两条序列在某一位匹配失败（字符不同），会得到扣分惩罚
* **gap open**：新添加一个空位的惩罚 q（类似起步价，一般比较大）
* **gap extension**：每再添加一个空位的惩罚 r（一般比较小）
* **全局比对（global alignment）**：比较两条完整的序列
* **局部比对（local alignment）**：比较两条序列的某一片段
* **traceback**：回溯，在完成 DP 后，通过终点进行回溯到原点，可以得到一个最优添加空位的序列方案
* **CIGAR**：通过字符串记录回溯得到的序列（如 6M2I2M），M（match/mismatch）比对的字符、I（insertion）target 插入的空格、D（deletion）query插入的空格，矩阵中纵轴 i 表示 target，横轴 j 表示 Query，回溯的时候往上走就多个 D，往左走就多个 I。用这个 CIGAR 字符串来看参考序列，就可以知道在查询序列相对于参考序列发生了哪些变化。如果是以 Target 为标准，Query 在这里多长了字符，CIGAR 就记为 I；Query 在这里多加了空位，CIGAR 就记为 D。

------

## 在 `ksw2_gg.c` 中

### 主要的函数：

```c
int ksw_gg(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
```

#### 函数中传入的参数

| 参数名                           | 角色                                              | 作用                                       |
| -------------------------------- | ------------------------------------------------- | :----------------------------------------- |
| `km`                             | 内存池地址                                        | 负责算法内部的高速内存分配                 |
| `qlen`, `query`                  | `Query` 串的长度及其本身地址                      | 纵轴 $i$，待检测的高通量 `Reads`           |
| `tlen`, `target`                 | `Target` 串的长度及其本身地址                     | 横轴 $j$，权威的基因组参考标准             |
| `m`, `mat`                       | 打分矩阵（维度与指针）                            | 规定 $5 \times 5$ 的 `Match/Mismatch` 规则 |
| `gapo`, `gape`                   | 空位罚分（开启惩罚和延伸惩罚）                    | 开启和延续 `Gap` 的扣分底牌                |
| `w`                              | 带宽                                              | 框定 DP 的带状剪枝边界 `[beg, end)`        |
| `m_cigar_`, `n_cigar_`, `cigar_` | 输出的结果`CIGAR`字符串（大小、片段数、二维指针） | 回溯路径（M/I/D）的信息                    |

### 状态定义

| 状态      | 物理意义                                                     | 视觉表现            | 状态转移（外层 $i$ 代表 Target 字符，内层 $j$ 代表 Query 字符） |
| --------- | ------------------------------------------------------------ | ------------------- | ------------------------------------------------------------ |
| $H(i, j)$ | **主状态 (Match / Mismatch) ；**当前 Target 字符与 Query 字符强行对齐的最优评分。 | Target: A Query: A  | 必须从**左上角对角线** $(i-1, j-1)$ 转移过来，加上当前碱基对齐的打分。 |
| $E(i, j)$ | **水平空位状态 (Insertion)；**相对于 Target 来说，Query 插入了字符，Target 插入空位。 | Target: - Query: A  | 必须从**左边格子** $(i, j-1)$ 转移过来。意味着在横向移动，持续往 Target 里灌空位。 |
| $F(i, j)$ | **垂直空位状态 (Deletion)；**相对于 Target 来说，Query 缺失了字符，Query 插入空位。 | Target: A Query:  - | 必须从**上边格子** $(i-1, j)$ 转移过来。意味着在纵向移动，持续往 Query 里灌空位。 |

### 结构体定义：

```c
typedef struct { int32_t h, e; } eh_t; 
```

* 用 $eh\_t.h$ 存储 $H(i)(j-1)$, $eh\_t.e$ 存储 $E(i+1)(j)$
  由于为了节省空间，用了滚动优化，所以当遍历到下一层i时，可以从 $eh$ 中直接读取到 $H(i-1)(j-1)$ 和 $E(i)(j)$ 以便后续状态转移

### 状态转移方程（NW）

$$
H(i,j) = \max \left\{ H(i-1,j-1) + S(i,j),\; E(i,j),\; F(i,j) \right\} \\
E(i+1,j) = \max \left\{ H(i,j) - gapo,\; E(i,j) \right\} - gape = \max \left\{ H(i,j) - gapoe,\; E(i,j) - gape \right\} \\
F(i,j+1) = \max \left\{ H(i,j) - gapo,\; F(i,j) \right\} - gape = \max \left\{ H(i,j) - gapoe,\; F(i,j) - gape \right\}
$$

## 在 `ksw2_sw.c` 中

### 状态转移方程（SW）

$$
H(i,j) = \max \left\{ 0,\; H(i-1,j-1) + S(i,j),\; E(i,j),\; F(i,j) \right\} \\
E(i+1,j) = \max \left\{ 0,\; H(i,j) - gapo,\; E(i,j) - gape \right\} \\
F(i,j+1) = \max \left\{ 0,\; H(i,j) - gapo,\; F(i,j) - gape \right\}
$$
