/*
本代码由ai完成，由于还是维护了绝对E、F、H，感觉没有体现做差分的意义，目前没有想到很好的手法在差分写法下维护 SW 的全局最优值，因此这份代码先暂时放这
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ksw2.h"
#include "ksw_sw_backtrack.h"

int ksw_sw2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, 
            int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, 
            int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int r, t, n_col, *off = 0;
	int8_t *s;
	uint8_t *p = 0;

	int max_score = 0;
	int end_r = -1, end_t = -1;
	int qe = q + e;
	
	int alloc_len = tlen + 2;

	// 绝对打分滚动矩阵
	int32_t *H0, *H1, *H2;
	int32_t *E1, *E0;
	int32_t *F1, *F0;

	// 专门用于维护回溯状态位 d 的差分数组（与绝对值计算同步更新）
	int8_t *x, *y;

	s = (int8_t*)kmalloc(km, tlen);
	if (w < 0) w = tlen > qlen ? tlen : qlen;
	n_col = w + 1 < tlen ? w + 1 : tlen;

	H0 = (int32_t*)kmalloc(km, alloc_len * sizeof(int32_t));
	H1 = (int32_t*)kmalloc(km, alloc_len * sizeof(int32_t));
	H2 = (int32_t*)kmalloc(km, alloc_len * sizeof(int32_t));
	E1 = (int32_t*)kmalloc(km, alloc_len * sizeof(int32_t));
	E0 = (int32_t*)kmalloc(km, alloc_len * sizeof(int32_t));
	F1 = (int32_t*)kmalloc(km, alloc_len * sizeof(int32_t));
	F0 = (int32_t*)kmalloc(km, alloc_len * sizeof(int32_t));

	x  = (int8_t*)kcalloc(km, alloc_len, 1);
	y  = (int8_t*)kcalloc(km, alloc_len, 1);

	// 初始化状态
	for (t = 0; t < alloc_len; ++t) {
		H0[t] = H1[t] = H2[t] = KSW_NEG_INF;
		E1[t] = E0[t] = KSW_NEG_INF;
		F1[t] = F0[t] = KSW_NEG_INF;
		x[t] = y[t] = 0;
	}

	if (m_cigar_ && n_cigar_ && cigar_) {
		p = (uint8_t*)kcalloc(km, (size_t)(qlen + tlen) * n_col, 1);
		off = (int*)kmalloc(km, (qlen + tlen) * sizeof(int));
	}

	// 主轴波面对角线扫描 (r = i + j)
	for (r = 0; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1;
		
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r - w + 1) >> 1) st = (r - w + 1) >> 1;
		if (en > (r + w) >> 1) en = (r + w) >> 1;

		if (p) off[r] = st;

		for (t = st; t <= en; ++t) {
			s[t] = mat[target[t] * m + query[r - t]];
		}

		uint8_t *pr = p ? p + (size_t)r * n_col : 0;

		// 边界安全清洗
		if (st > 0) {
			H0[st - 1] = E0[st - 1] = F0[st - 1] = KSW_NEG_INF;
		}
		H0[en + 1] = E0[en + 1] = F0[en + 1] = KSW_NEG_INF;

		for (t = st; t <= en; ++t) {
			int32_t score_match = s[t]; 
			int32_t score_del = KSW_NEG_INF;
			int32_t score_ins = KSW_NEG_INF;

			// ===== 1. 匹配分支 (Match) =====
			if (r >= 2 && t > 0) {
				if (H2[t - 1] > KSW_NEG_INF) {
					int32_t from_match = H2[t - 1] + s[t];
					if (from_match > score_match) score_match = from_match;
				}
			}

			// ===== 2. 删除分支 (E) =====
			if (r >= 1 && t > 0) {
				int32_t from_H = (H1[t - 1] > 0) ? H1[t - 1] - qe : KSW_NEG_INF;
				int32_t from_E = (E1[t - 1] > KSW_NEG_INF) ? E1[t - 1] - e : KSW_NEG_INF;
				
				E0[t] = (from_E > from_H) ? from_E : from_H;
				score_del = E0[t];
			} else {
				E0[t] = KSW_NEG_INF;
			}

			// ===== 3. 插入分支 (F) =====
			if (r >= 1) {
				int32_t from_H = (H1[t] > 0) ? H1[t] - qe : KSW_NEG_INF;
				int32_t from_F = (F1[t] > KSW_NEG_INF) ? F1[t] - e : KSW_NEG_INF;

				F0[t] = (from_F > from_H) ? from_F : from_H;
				score_ins = F0[t];
			} else {
				F0[t] = KSW_NEG_INF;
			}

			// ===== 4. SW 局部比对决策 =====
			int32_t H_val = 0; 
			uint8_t d = 4; 

			if (score_match > H_val) { H_val = score_match; d = 0; }
			if (score_del > H_val)   { H_val = score_del;   d = 1; }
			if (score_ins > H_val)   { H_val = score_ins;   d = 2; }

			H0[t] = H_val;

			if (H_val > max_score) {
				max_score = H_val;
				end_r = r;
				end_t = t;
			}

			// ===== 5. 差分生成与回溯标志维护 =====
			// 即使 H_val 截断为 0，我们也需要就地维护 x[t] 和 y[t] 的差分值，供下一轮边界判定
			if (H1[t - 1] > KSW_NEG_INF && E0[t] > KSW_NEG_INF) {
				int32_t diff_x = H1[t - 1] - E0[t];
				x[t] = diff_x > 127 ? 127 : (diff_x < -128 ? -128 : diff_x);
			} else {
				x[t] = 127; // 趋于无穷远，代表不从此处延续
			}

			if (H1[t] > KSW_NEG_INF && F0[t] > KSW_NEG_INF) {
				int32_t diff_y = H1[t] - F0[t];
				y[t] = diff_y > 127 ? 127 : (diff_y < -128 ? -128 : diff_y);
			} else {
				y[t] = 127;
			}

			if (p) {
				if (H_val == 0) {
					d = 0xff; // SW Restart 状态斩断
				} else {
					// 通过上一轮遗留下来的差分数组直接判定当前状态是否由 Gap 延续而来
					if (r >= 1 && t > 0 && x[t - 1] <= q) {
						d |= 0x08; // 差分判定：上一个格子传过来的 E 属于延续状态
					}
					if (r >= 1 && y[t] <= q) {
						d |= 0x10; // 差分判定：上一个格子传过来的 F 属于延续状态
					}
				}
				pr[t - st] = d;
			}
		}

		if (st > 0) H1[st - 1] = E1[st - 1] = F1[st - 1] = KSW_NEG_INF;
		H1[en + 1] = E1[en + 1] = F1[en + 1] = KSW_NEG_INF;

		// 滚动轮替
		int32_t *tmp_H = H2; H2 = H1; H1 = H0; H0 = tmp_H;
		int32_t *tmp_E = E1; E1 = E0; E0 = tmp_E;
		int32_t *tmp_F = F1; F1 = F0; F0 = tmp_F;
	}

	kfree(km, s); kfree(km, x); kfree(km, y);
	kfree(km, H0); kfree(km, H1); kfree(km, H2);
	kfree(km, E1); kfree(km, E0);
	kfree(km, F1); kfree(km, F0);

	if (m_cigar_ && n_cigar_ && cigar_) {
		if (max_score > 0 && end_r >= 0 && end_t >= 0) {
			int end_q = end_r - end_t;
			if (end_t >= tlen) end_t = tlen - 1;
			if (end_q < 0) end_q = 0;
			if (end_q >= qlen) end_q = qlen - 1;
			ksw_sw_backtrack(km, 1, 0, 0, p, off, 0, n_col, end_t, end_q, m_cigar_, n_cigar_, cigar_);
		} else {
			*n_cigar_ = 0;
			*cigar_ = NULL;
		}
		kfree(km, p); kfree(km, off);
	}

	return max_score;
}
