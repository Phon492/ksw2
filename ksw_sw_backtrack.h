#ifndef KSW_SW_BACKTRACK_H_
#define KSW_SW_BACKTRACK_H_

#include <stdint.h>
#include "ksw2.h"

static inline void ksw_sw_backtrack(void *km, int is_rot, int is_rev, int min_intron_len, const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0,
									int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
	uint32_t *cigar = *cigar_, tmp;
	while (i >= 0 && j >= 0) { 
		int force_state = -1;
		if (is_rot) {
			r = i + j;
			if (i < off[r]) force_state = 2;
			if (off_end && i > off_end[r]) force_state = 1;
			tmp = force_state < 0? p[(size_t)r * n_col + i - off[r]] : 0;
		} else {
			if (j < off[i]) force_state = 2;
			if (off_end && j > off_end[i]) force_state = 1;
			tmp = force_state < 0? p[(size_t)i * n_col + j - off[i]] : 0;
		}
		if(tmp == 0xff) break;
		if (state == 0) state = tmp & 7; 
		else if (!(tmp >> (state + 2) & 1)) state = 0; 
		if (state == 0) state = tmp & 7;
		if (force_state >= 0) state = force_state;
		if (state == 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_MATCH, 1), --i, --j;
		else if (state == 1 || (state == 3 && min_intron_len <= 0)) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_DEL, 1), --i;
		else if (state == 3 && min_intron_len > 0) cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_N_SKIP, 1), --i;
		else cigar = ksw_push_cigar(km, &n_cigar, &m_cigar, cigar, KSW_CIGAR_INS, 1), --j;
	}
	// 删去了跑到边界的步骤，从全局比对变为局部比对
	if (!is_rev)
		for (i = 0; i < n_cigar>>1; ++i) 
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
	*m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

#endif
