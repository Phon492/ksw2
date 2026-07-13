#include <stdio.h> 
#include "ksw2.h"
#include "ksw_sw_backtrack.h"

#ifdef __SSE2__
#include <emmintrin.h>

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

int ksw_sw2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int r, t, n_col, n_col_, *off, tlen_, last_st, last_en;
	uint8_t *qr, *mem, *mem2;
	__m128i *u, *v, *x, *y, *s, *p;
	__m128i q_, qe2_, zero_, flag1_, flag2_, flag8_, flag16_;

    int32_t *H0 = 0, *H1 = 0, *H2 = 0;
    int max_score = 0, max_r = 0, max_t = 0;

	zero_   = _mm_set1_epi8(0);
	q_      = _mm_set1_epi8(q);
	qe2_    = _mm_set1_epi8((q + e) * 2);
	flag1_  = _mm_set1_epi8(1);
	flag2_  = _mm_set1_epi8(2);
	flag8_  = _mm_set1_epi8(0x08);
	flag16_ = _mm_set1_epi8(0x10);

	if (w < 0) w = tlen > qlen? tlen : qlen;
	n_col = w + 1 < tlen? w + 1 : tlen;
	tlen_ = (tlen + 15) / 16;
	n_col_ = (n_col + 15) / 16 + 1;
	n_col = n_col_ * 16;

	mem = (uint8_t*)kcalloc(km, tlen_ * 5 + 1, 16);
	u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16 字节对齐的位运算技巧
	v = u + tlen_, x = v + tlen_, y = x + tlen_, s = y + tlen_;
	qr = (uint8_t*)kcalloc(km, qlen, 1);
	mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
	p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
	off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int));

    H0 = (int32_t*)kcalloc(km, tlen + 1, sizeof(int32_t));
    H1 = (int32_t*)kcalloc(km, tlen + 1, sizeof(int32_t));
    H2 = (int32_t*)kcalloc(km, tlen + 1, sizeof(int32_t));

	for (t = 0; t < qlen; ++t)
		qr[t] = query[qlen - 1 - t];

	for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1, st0, en0, st_, en_;
		int8_t x1, v1;
		__m128i x1_, v1_, *pr;
		
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-w+1)>>1) st = (r-w+1)>>1;
		if (en > (r+w)>>1) en = (r+w)>>1;
		st0 = st, en0 = en;
		st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
		off[r] = st;
		
		if (st > 0) {
			if (st - 1 >= last_st && st - 1 <= last_en)
				x1 = ((uint8_t*)x)[st - 1], v1 = ((uint8_t*)v)[st - 1];
			else x1 = v1 = 0;
		} else x1 = 0, v1 = 0;
		if (en >= r) ((uint8_t*)y)[r] = 0, ((uint8_t*)u)[r] = 0;
		
		for (t = st0; t <= en0; ++t)
			((uint8_t*)s)[t] = mat[target[t] * m + qr[t + qlen - 1 - r]];
		
		x1_ = _mm_cvtsi32_si128(x1);
		v1_ = _mm_cvtsi32_si128(v1);
		st_ = st>>4, en_ = en>>4;
		pr = p + (size_t)r * n_col_ - st_;
		for (t = st_; t <= en_; ++t) {
            __m128i d, z, a, b, xt1, vt1, ut, tmp;
            int base = t << 4;   // 当前块对应的起始 t 索引

            // SW：z 从 0 开始（允许重新开始）
            z = zero_;

            // s_val = S + 2q + 2e
            __m128i s_val = _mm_add_epi8(_mm_load_si128(&s[t]), qe2_);

            // 移位 x：模拟 x(r-1, t-1)
            xt1 = _mm_load_si128(&x[t]);
            tmp = _mm_srli_si128(xt1, 15);      // 保存最后一个元素
            xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); // 左移 + 拼接
            x1_ = tmp;

            // 移位 v：模拟 v(r-1, t-1)
            vt1 = _mm_load_si128(&v[t]);
            tmp = _mm_srli_si128(vt1, 15);
            vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_);
            v1_ = tmp;

            // 计算 a = x+v, b = y+u
            a = _mm_add_epi8(xt1, vt1);
            ut = _mm_load_si128(&u[t]);
            b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);

            // SW：z = max(0, s_val, a, b)
            d = _mm_set1_epi8(0xff);  
            // 比较 s_val 与 z（z 初始为 0）
            __m128i mask_s = _mm_cmpgt_epi8(s_val, z);
#ifdef __SSE4_1__
            z = _mm_max_epi8(z, s_val);
            d = _mm_blendv_epi8(d, _mm_set1_epi8(0), mask_s);
#else
            z = _mm_max_epu8(z, s_val); 
            d = _mm_or_si128(_mm_andnot_si128(mask_s, d),
                             _mm_and_si128(mask_s, _mm_set1_epi8(0)));
#endif
            // a 与 z
            __m128i mask_a = _mm_cmpgt_epi8(a, z);
#ifdef __SSE4_1__
            z = _mm_max_epi8(z, a);
            d = _mm_blendv_epi8(d, _mm_set1_epi8(1), mask_a);
#else
            z = _mm_max_epu8(z, a);
            d = _mm_or_si128(_mm_andnot_si128(mask_a, d),
                             _mm_and_si128(mask_a, _mm_set1_epi8(1)));
#endif
            // 比较 b 与 z
            __m128i mask_b = _mm_cmpgt_epi8(b, z);
#ifdef __SSE4_1__
            z = _mm_max_epi8(z, b);
            d = _mm_blendv_epi8(d, _mm_set1_epi8(2), mask_b);
#else
            z = _mm_max_epu8(z, b);
            d = _mm_or_si128(_mm_andnot_si128(mask_b, d),
                             _mm_and_si128(mask_b, _mm_set1_epi8(2)));
#endif

            // 更新 u 和 v
            _mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));
            _mm_store_si128(&v[t], _mm_sub_epi8(z, ut));

            // 更新 x 和 y，并设置延续位
            // x = max(0, a - z + q)，y = max(0, b - z + q)
            z = _mm_sub_epi8(z, q_);
            a = _mm_sub_epi8(a, z);
            b = _mm_sub_epi8(b, z);
            __m128i mask_a_pos = _mm_cmpgt_epi8(a, zero_);
            d = _mm_or_si128(d, _mm_and_si128(flag8_, mask_a_pos));  
            _mm_store_si128(&x[t], _mm_and_si128(a, mask_a_pos));
            __m128i mask_b_pos = _mm_cmpgt_epi8(b, zero_);
            d = _mm_or_si128(d, _mm_and_si128(flag16_, mask_b_pos)); 
            _mm_store_si128(&y[t], _mm_and_si128(b, mask_b_pos));

            _mm_store_si128(&pr[t], d);   

			for (int tt = base; tt < base + 16 && tt < tlen; ++tt) {
                int32_t h = 0;
                
                if (tt < st0 || tt > en0) {
                    H0[tt] = 0;
                    continue;
                }

                // 匹配: H(r,t) = H(r-2, t-1) + S(i,j)
                if (tt > 0 && (tt - 1) >= (r - 2 - qlen + 1)) {
                    int32_t match_val = H2[tt - 1] + ((int8_t*)s)[tt];
                    h = match_val > h? match_val : h;
                } else {
                    int32_t match_val = ((int8_t*)s)[tt];
                    h = match_val > h? match_val : h;
                }

                // 删除: H(r,t) = H(r-1, t) - q - e
                if (H1[tt] >= 0) {
                    int32_t del_val = H1[tt] - (q + e);
                    h = del_val > h? del_val : h;
                }

                // 插入: H(r,t) = H(r-1, t-1) - q - e
                if (tt > 0) {
                    int32_t ins_val = H1[tt - 1] - (q + e);
                    h = ins_val > h? ins_val : h;
                }

                H0[tt] = h;

                if (h > max_score) {
                    max_score = h;
                    max_r = r;
                    max_t = tt;
                }
            }
	        int32_t *H_ = H2;
	        H2 = H1;
	        H1 = H0;
	        H0 = H_;
	        last_st = st;
	        last_en = en;
		}
		
		//for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%x\n", r, t, ((uint8_t*)u)[t], ((uint8_t*)v)[t], ((uint8_t*)x)[t], ((uint8_t*)y)[t], ((uint8_t*)(p + r * n_col_))[t-st]); // for debugging
	}
	kfree(km, mem); kfree(km, qr);
    if (m_cigar_ && n_cigar_ && cigar_) {
        int end_i = max_t;
        int end_j = max_r - max_t;
        if (end_i < 0) end_i = 0;
        if (end_i >= tlen) end_i = tlen - 1;
        if (end_j < 0) end_j = 0;
        if (end_j >= qlen) end_j = qlen - 1;

        ksw_sw_backtrack(km, 1, 0, 0, (uint8_t*)p, off, 0, n_col,
                         end_i, end_j, m_cigar_, n_cigar_, cigar_);
        kfree(km, mem2);
        kfree(km, off);
        if (H0) kfree(km, H0);
        if (H1) kfree(km, H1);
        if (H2) kfree(km, H2);
    } else {
        kfree(km, mem2);
        kfree(km, off);
        if (H0) kfree(km, H0);
        if (H1) kfree(km, H1);
        if (H2) kfree(km, H2);
    }

    return max_score;
}
#endif // __SSE2__
