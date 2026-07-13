#include <stdio.h>
#include "ksw2.h"
#include "ksw_sw_backtrack.h"

int ksw_sw2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int qe = q + e, qe2 = qe + qe, r, t, n_col, *off = 0;
	int8_t *u, *v, *x, *y, *s;
	uint8_t *p = 0, *qr;
	int32_t *H0 = 0, *H1 = 0, *H2 = 0; // r, r-1, r-2
	int max_score = 0;
	int end_r = -1, end_t = -1;

	u = (int8_t*)kcalloc(km, tlen + 1, 1);
	v = (int8_t*)kcalloc(km, tlen + 1, 1);
	x = (int8_t*)kcalloc(km, tlen + 1, 1);
	y = (int8_t*)kcalloc(km, tlen + 1, 1);
	s = (int8_t*)kmalloc(km, tlen);
	qr = (uint8_t*)kmalloc(km, qlen);
	if (w < 0) w = tlen > qlen? tlen : qlen;
	n_col = w + 1 < tlen? w + 1 : tlen;
	H0 = (int32_t*)kcalloc(km, tlen + 1, sizeof(int32_t));
	H1 = (int32_t*)kcalloc(km, tlen + 1, sizeof(int32_t));
	H2 = (int32_t*)kcalloc(km, tlen + 1, sizeof(int32_t));	

	if (m_cigar_ && n_cigar_ && cigar_) {
		p = (uint8_t*)kcalloc(km, (size_t)(qlen + tlen) * n_col, 1);
		off = (int*)kmalloc(km, (qlen + tlen) * sizeof(int));
	}

	for (t = 0; t < qlen; ++t)
		qr[t] = query[qlen - 1 - t];

	for (r = 0; r < qlen + tlen - 1; ++r) {
		int st = 0, en = tlen - 1;
		int8_t x1, v1;
		
		if (st < r - qlen + 1) st = r - qlen + 1;
		if (en > r) en = r;
		if (st < (r-w+1)>>1) st = (r-w+1)>>1; 
		if (en > (r+w)>>1) en = (r+w)>>1; 
		
		if (st != 0) {
			if (r > st + st + w - 1) x1 = v1 = 0;
			else x1 = x[st-1], v1 = v[st-1]; 
		} else x1 = 0, v1 = 0;
		if (en != r) {
			if (r < en + en - w - 1) y[en] = u[en] = 0; 
		} else y[r] = 0, u[r] = 0;
		
		for (t = st; t <= en; ++t)
			s[t] = mat[target[t] * m + qr[t + qlen - 1 - r]];
		
		if (m_cigar_ && n_cigar_ && cigar_) {
			uint8_t *pr = p + (size_t)r * n_col;
			off[r] = st;
			for (t = st; t <= en; ++t) {
				int32_t h = 0;
				// 匹配：H(r,t) = H(r-2, t-1) + S(i,j)
				if (t - 1 >= 0 && (t - 1) >= ((r - 2) - qlen + 1)) {  // r - qlen + 1 <= t <= r; t <- t-1, r <- r-2
					h = H2[t - 1] + s[t] > h ? H2[t - 1] + s[t] : h;
				}
                // 删除：H(r,t) = H(r-1, t) - q - e
                if (H1[t] >= 0) {
                    h = H1[t] - q - e > h? H1[t] - qe : h;
                }
                // 插入：H(r,t) = H(r, t-1) - q - e
                if (t - 1 >= 0) {
                    h = H0[t - 1] - q - e > h? H0[t - 1] - qe : h;
                }

				if (h > max_score) {
					max_score = h;
					end_t = t;
					end_r = r;
				}

				H0[t] = h;

				uint8_t d;
				if(h == 0) d = 0xff;
				
				int8_t u1;
				int8_t z = s[t] + qe2;
				int8_t a = x1   + v1;
				int8_t b = y[t] + u[t];
				d = a > z? 1 : 0; 
				z = a > z? a : z;
				d = b > z? 2 : d;
				z = b > z? b : z;

				d = z > 0? d : 0xff;
				
				u1 = u[t];            
				u[t] = z - v1;     
				v1 = v[t];            
				v[t] = z - u1;       
				z -= q;
				a -= z;
				b -= z;
				x1 = x[t];           
				d   |= a > 0? 0x08 : 0;
				x[t] = a > 0? a    : 0;
				d   |= b > 0? 0x10 : 0;
				y[t] = b > 0? b    : 0; 
				pr[t - st] = d;
			}
			int32_t *tmp = H2;
			H2 = H1;
			H1 = H0;
			H0 = tmp;
		} else {
			for (t = st; t <= en; ++t) {
				int32_t h = 0;
				// 匹配：H(r,t) = H(r-2, t-1) + S(i,j)
				if (t - 1 >= 0 && (t - 1) >= ((r - 2) - qlen + 1)) {  
					h = H2[t - 1] + s[t] > h ? H2[t - 1] + s[t] : h;
				}
                // 删除：H(r,t) = H(r-1, t) - q - e
                if (H1[t] >= 0) {
                    h = H1[t] - q - e > h? H1[t] - qe : h;
                }
                // 插入：H(r,t) = H(r, t-1) - q - e
                if (t - 1 >= 0) {
                    h = H0[t - 1] - q - e > h? H0[t - 1] - qe : h;
                }

				if (h > max_score) {
					max_score = h;
					end_t = t;
					end_r = r;
				}

				H0[t] = h;

				int8_t u1;
				int8_t z = s[t] + qe2;
				int8_t a = x1   + v1;
				int8_t b = y[t] + u[t];
				z = a > z? a : z;
				z = b > z? b : z;
				
				u1 = u[t];            
				u[t] = z - v1;     
				v1 = v[t];            
				v[t] = z - u1;       
				z -= q;
				a -= z;
				b -= z;
				x1 = x[t];           
				x[t] = a > 0? a    : 0;
				y[t] = b > 0? b    : 0; 
			}
			int32_t *tmp = H2;
			H2 = H1;
			H1 = H0;
			H0 = tmp;
		}
	}
	kfree(km, u); kfree(km, v); kfree(km, x); kfree(km, y); kfree(km, s); kfree(km, qr);
	if (m_cigar_ && n_cigar_ && cigar_) {
		int end_q = end_r - end_t;
		if (end_t < 0) end_t = 0;
		if (end_t >= tlen) end_t = tlen - 1;
		if (end_q < 0) end_q = 0;
		if (end_q >= qlen) end_q = qlen - 1;
		ksw_sw_backtrack(km, 1, 0, 0, p, off, 0, n_col, end_t, end_q, m_cigar_, n_cigar_, cigar_);
		kfree(km, p); kfree(km, off);
	}
	return max_score;
}
