#include <stdio.h>
#include <ksw2.h>

typedef struct { int32_t h, e; } eh_t; 

int ksw_sw(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	eh_t *eh;
	int8_t *qp;
	int32_t i, j, k, gapoe = gapo + gape, n_col, *off = 0;
	uint8_t *z = 0;
	int max_score = 0;
	int end_t = -1, end_q = -1;

	if (w < 0) w = tlen > qlen? tlen : qlen;
	n_col = qlen < 2*w+1? qlen : 2*w+1;
	qp = (int8_t*)kmalloc(km, qlen * m);
	eh = (eh_t*)kcalloc(km, qlen + 1, 8);
	if (m_cigar_ && n_cigar_ && cigar_) {
		*n_cigar_ = 0;
		z = (uint8_t*)kmalloc(km, (size_t)n_col * tlen);
		off = (int32_t*)kcalloc(km, tlen, 4);
	}

	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	
	for (j = 0; j <= qlen; ++j) 
		eh[j].h = 0, eh[j].e = 0;

	for (i = 0; i < tlen; ++i) {
		int32_t f = 0, h1 = 0, st, en;
		int8_t *q = &qp[target[i] * qlen];

		st = i > w? i - w : 0;
		en = i + w + 1 < qlen? i + w + 1 : qlen;

		if (m_cigar_ && n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			off[i] = st;
			for (j = st; j < en; ++j){
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				uint8_t d; 
				p->h = h1;
				h += q[j];

				d = h >= e? 0 : 1;
				h = h >= e? h : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;

				if (h > max_score) {
					max_score = h;
					end_t = i;
					end_q = j;
				}

				h1 = h > 0? h : 0;

				h -= gapoe;
				e -= gape;
				d |= e > h? 0x08 : 0;
				e  = e > h? e    : h; 
				e  = e > 0? e    : 0;
				p->e = e;

				f -= gape;
				f  = f > 0? f : 0;
				f  = f > h? f    : h;
				d |= f > h? 0x10 : 0;

				zi[j - st] = d; 
			} 
		} else {
			for (j = st; j < en; ++j) {
				eh_t *p = &eh[j];
				int32_t h = p->h, e = p->e;
				p->h = h1;
				h += q[j];
				
				h = h >= e? h : e;
				h = h >= f? h : f;


				if (h > max_score) {
					max_score = h;
					end_t = i;
					end_q = j;
				}

				h1 = h > 0? h : 0;

				h -= gapoe;
				e -= gape;
				e  = e > 0? e : 0;
				e  = e > h? e : h;
				p->e = e;
				
				f -= gape;
				f  = f > 0? f : 0;
				f  = f > h? f : h;
			}			
		}
		eh[en].h = h1, eh[en].e = KSW_NEG_INF;
	}
	kfree(km, qp); kfree(km, eh);
	if (m_cigar_ && n_cigar_ && cigar_) {
		ksw_backtrack(km, 0, 0, 0, z, off, 0, n_col, end_t, end_q, m_cigar_, n_cigar_, cigar_);
		kfree(km, z);
		kfree(km, off);
	}
	return max_score;
}
