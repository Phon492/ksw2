#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// 声明 ksw_sw2，测试其他 SW 下代码时，修改此条声明即可
int ksw_sw2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target,
            int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w,
            int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

void print_cigar(uint32_t *cigar, int n_cigar) {
    if (n_cigar <= 0 || !cigar) {
        printf("CIGAR: *\n");
        return;
    }
    printf("CIGAR: ");
    for (int i = 0; i < n_cigar; ++i) {
        int len = cigar[i] >> 4;
        int op = cigar[i] & 0xf;
        char op_char = "MIDNSH"[op];
        printf("%d%c", len, op_char);
    }
    printf("\n");
}

// 辅助函数：运行单个测试
void run_test(const char *desc,
              const uint8_t *query, int qlen,
              const uint8_t *target, int tlen,
              int8_t q, int8_t e,
              int expected_score, const char *expected_cigar) {
    int8_t m = 5;
    int8_t match = 2, mismatch = -4;
    int w = 10;

    // 构建打分矩阵
    int8_t mat[25];
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (i == 4 || j == 4) mat[i*5 + j] = 0;
            else mat[i*5 + j] = (i == j) ? match : mismatch;
        }
    }

    printf("===== %s =====\n", desc);
    printf("query:  ");
    for (int i = 0; i < qlen; ++i) printf("%d ", query[i]);
    printf("\n");
    printf("target: ");
    for (int i = 0; i < tlen; ++i) printf("%d ", target[i]);
    printf("\n");
    printf("罚分: open=%d, ext=%d\n", q, e);

    // 纯算分
    int score = ksw_sw2(NULL, qlen, query, tlen, target, m, mat, q, e, w, NULL, NULL, NULL);
    printf("纯算分得分: %d  (期望 %d)\n", score, expected_score);

    // 带回溯
    int m_cigar = 0, n_cigar = 0;
    uint32_t *cigar = NULL;
    int score_b = ksw_sw2(NULL, qlen, query, tlen, target, m, mat, q, e, w,
                          &m_cigar, &n_cigar, &cigar);
    printf("带回溯得分: %d  (期望 %d)\n", score_b, expected_score);
    print_cigar(cigar, n_cigar);
    printf("期望 CIGAR: %s\n\n", expected_cigar);
    if (cigar) free(cigar);
}

int main() {
    // ---------- 测试 1：局部性（两端截断） ----------
    {
        // query: AAAAA, target: TTAAAAATT
        uint8_t query1[]  = {0, 0, 0, 0, 0};
        uint8_t target1[] = {3, 3, 0, 0, 0, 0, 0, 3, 3};
        run_test("局部性（两端截断）",
                 query1, sizeof(query1),
                 target1, sizeof(target1),
                 3, 1,   // 罚分任意，因为不涉及 gap
                 10, "5M");
    }

    // ---------- 测试 2：仿射 gap 罚分（连续缺失） ----------
    {
        // query: AACCGT, target: AAGT
        // 使用低罚分使 gap 比对更优：q=1, e=1
        uint8_t query2[]  = {0, 0, 1, 1, 2, 3};
        uint8_t target2[] = {0, 0, 2, 3};
        run_test("仿射 gap 罚分（连续缺失）",
                 query2, sizeof(query2),
                 target2, sizeof(target2),
                 1, 1,    // q=1, e=1
                 5, "2M2I2M");
    }

    // 可选：添加一个独立 gap 对比测试，展示仿射优势
    // 这里不赘述。

    return 0;
}
