/*
 * Copyright 2020 Amazon.com, Inc. or its affiliates. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License").
 * You may not use this file except in compliance with the License.
 * A copy of the License is located at
 *
 * http://aws.amazon.com/apache2.0
 *
 * or in the "license" file accompanying this file. This file is distributed
 * on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing
 * permissions and limitations under the License.
 * The license is detailed in the file LICENSE.md, and applies to this file.
 *
 * Written by Nir Drucker, Shay Gueron, and Dusan Kostic,
 * AWS Cryptographic Algorithms Group.
 * (ndrucker@amazon.com, gueron@amazon.com, dkostic@amazon.com)
 *
 * [1] The optimizations are based on the description developed in the paper:
 *     Drucker, Nir, and Shay Gueron. 2019. “A Toolbox for Software Optimization
 *     of QC-MDPC Code-Based Cryptosystems.” Journal of Cryptographic Engineering,
 *     January, 1–17. https://doi.org/10.1007/s13389-018-00200-4.
 *
 * [2] The decoder algorithm is the Black-Gray decoder in
 *     the early submission of CAKE (due to N. Sandrier and R Misoczki).
 *
 * [3] The analysis for the constant time implementation is given in
 *     Drucker, Nir, Shay Gueron, and Dusan Kostic. 2019.
 *     “On Constant-Time QC-MDPC Decoding with Negligible Failure Rate.”
 *     Cryptology EPrint Archive, 2019. https://eprint.iacr.org/2019/1289.
 *
 * [4] it was adapted to BGF in:
 *     Drucker, Nir, Shay Gueron, and Dusan Kostic. 2019.
 *     “QC-MDPC decoders with several shades of gray.”
 *     Cryptology EPrint Archive, 2019. To be published.
 *
 * [5] Chou, T.: QcBits: Constant-Time Small-Key Code-Based Cryptography.
 *     In: Gier-lichs, B., Poschmann, A.Y. (eds.) Cryptographic Hardware
 *     and Embedded Systems– CHES 2016. pp. 280–300. Springer Berlin Heidelberg,
 *     Berlin, Heidelberg (2016)
 *
 * [6] The rotate512_small funciton is a derivative of the code described in:
 *     Guimarães, Antonio, Diego F Aranha, and Edson Borin. 2019.
 *     “Optimized Implementation of QC-MDPC Code-Based Cryptography.”
 *     Concurrency and Computation: Practice and Experience 31 (18):
 *     e5089. https://doi.org/10.1002/cpe.5089.
 */

#include "decode.h"
#include "gf2x.h"
#include "sampling.h"
#include "utilities.h"
#include <math.h>
#include <string.h>
#include <time.h>

// Decoding (bit-flipping) parameter
#ifdef BG_DECODER
#  if(LEVEL == 1)
#    define MAX_IT 3
#  elif(LEVEL == 3)
#    define MAX_IT 4
#  elif(LEVEL == 5)
#    define MAX_IT 7
#  else
#    error "Level can only be 1/3/5"
#  endif
#elif defined(BGF_DECODER)
#  if(LEVEL == 1)
#    define MAX_IT 5
#  elif(LEVEL == 3)
#    define MAX_IT 6
#  elif(LEVEL == 5)
#    define MAX_IT 7
#  else
#    error "Level can only be 1/3/5"
#  endif
#endif

#define EQ_COLUMN         201 // 索引矩阵列数
#define ROW               R_BITS
#define X                 EQ_COLUMN - 1
#define N                 2 * R_BITS
#define GUSS_INDEX_COLUMN 100
#define GUSS_BLOCK        8

// 定义是否包括 stepII 和 stepIII 的集合, 0 代表不包括, 1 代表包括
#define STEP23_ON 0
// 定义解方程的 delta 大小
#define DELTA_EQ 6
// 定义 stepII 和 stepIII 中取的集合 delta 大小
#define DELTA_STEP23 0

// 0 使用拟合方法，1 使用论文方法
#define TH_SELECT 0

// 利用论文中的方法计算 th
_INLINE_ uint8_t
compute_th_R(IN uint16_t sk_wlist_all_0[][DV],
             IN uint16_t sk_wlist_all_1[][DV],
             // IN const split_e_t  *e,
             IN const uint16_t    T,
             IN const split_e_t  *R_e,
             IN const syndrome_t *s)
{
  // ---- 1.计算 s 的重量 ----
  uint16_t s_weight = r_bits_vector_weight((const r_t *)s->qw);

  // ---- 2. 计算 X ----
  // 根据 H 每一行的索引去找 R_e 中是 1 的位置，计算重量
  double_t x          = 0;
  uint16_t tmp_weight = 0;
  uint16_t L[R_BITS]  = {0};
  uint8_t  mask_1     = 1;

  for(uint16_t r_R_BIT = 0; r_R_BIT < R_BITS; r_R_BIT++)
  {
    for(uint8_t i_DV = 0; i_DV < DV; i_DV++)
    {
      // 将索引位置除 8 找 e 对应字节，用 mod 8 找 对应位置
      if((R_e->val[0].raw[(sk_wlist_all_0[r_R_BIT][i_DV] / 8)] &
          (mask_1 << (sk_wlist_all_0[r_R_BIT][i_DV] % 8))) != 0)
      {
        // 如果与结果不为 0 则重量加 1
        tmp_weight += 1;
      }
      if((R_e->val[1].raw[(sk_wlist_all_1[r_R_BIT][i_DV] / 8)] &
          (mask_1 << (sk_wlist_all_1[r_R_BIT][i_DV] % 8))) != 0)
      {
        // 如果与结果不为 0 则重量加 1
        tmp_weight += 1;
      }
    }
    // 将当前行的重量作为下标保存到 L 中
    L[tmp_weight] += 1;
    // tmp_weight 清 0
    tmp_weight = 0;
  }

  // 将 L 中奇数索引进行运算
  for(uint16_t i_l = 1; i_l < R_BITS; i_l = i_l + 2)
  {
    if(L[i_l] == 0)
    {
      break;
    }
    double_t A = 0;
    double_t B = 0;
    double_t C = 0;

    // 奇数运算
    if(i_l % 2 == 1)
    {
      // C(w,l)
      for(int i_c = 0; i_c < i_l; i_c++)
      {
        A += log10((double_t)(2 * DV - i_c) / (i_l - i_c));
      }
      // C(n-w,t-l)
      for(int i_c = 0; i_c < T - i_l; i_c++)
      {
        B += log10((double_t)(2 * R_BITS - 2 * DV - i_c) / (T - i_l - i_c));
      }
      // C(n,t)
      for(int i_c = 0; i_c < T; i_c++)
      {
        C += log10((double_t)(2 * R_BITS - i_c) / (T - i_c));
      }
      // 求当前 x
      x += (i_l - 1) * R_BITS * pow(10, (A + B - C));
    }
  }

  // ---- 3. 计算 T ----
  double_t pai_0 = (double_t)(2 * DV * s_weight - x) / ((2 * R_BITS - T) * DV);
  double_t pai_1 = (double_t)(s_weight + x) / (T * DV);

  double_t th = (log10((double_t)(2 * R_BITS - T) / T) +
                 DV * log10((1 - pai_0) / (1 - pai_1))) /
                (log10(pai_1 / pai_0) + log10((1 - pai_0) / (1 - pai_1)));

  // ---- test ---- 打印参数
  // printf("x = %f\n", x);
  // printf("s_weight = %u\n", s_weight);
  // printf("pai_0 = %f\n", pai_0);
  // printf("pai_1 = %f\n", pai_1);
  // printf("th = %f\n", th);
  // printf("th_up = %u\n",(uint8_t)th);
  return (uint8_t)th + 1;
}

// 用于交换两个数组
_INLINE_ void
swap(OUT uint8_t *a, OUT uint8_t *b, uint32_t eq_j, uint32_t guss_j_num)
{
  uint8_t tmp_guss[guss_j_num];
  for(uint16_t change_i = eq_j; change_i < guss_j_num; change_i++)
  {
    tmp_guss[change_i] = a[change_i];
    a[change_i]        = b[change_i];
    b[change_i]        = tmp_guss[change_i];
  }
}

// 64 位异或
_INLINE_ ret_t
xor_8(OUT uint8_t      *res,
      IN const uint8_t *a,
      IN const uint8_t *b,
      IN const uint64_t bytelen,
      IN const uint64_t r_bytelen)
{
  for(uint64_t i = r_bytelen; i < bytelen; i++)
  {
    res[i] = a[i] ^ b[i];
  }
  return SUCCESS;
}

// 利用解出来的 b 和 ct 还原 fm(ct_verify)
_INLINE_ void
solving_equations_mf(IN OUT ct_t *ct_verify, IN uint16_t b[])
{
  // 放 0 用 '与', 放 1 用 '或'
  // 定义 11111111 和 00000001 用于计算
  uint8_t mask_255 = 255;
  uint8_t mask_1   = 1;
  int     bit_u    = 8;
  // 对第一组操作
  for(int i_v = 0; i_v < R_BITS; i_v++)
  {
    if(b[i_v] != 0)
    {
      b[i_v] = b[i_v] % 2;
      if(b[i_v] == 0)
      {
        // 用与操作
        ct_verify->val[0].raw[i_v / bit_u] =
            (mask_255 ^ (mask_1 << (i_v % bit_u))) &
            ct_verify->val[0].raw[i_v / bit_u];
      }
      else
      {
        // 用或操作
        ct_verify->val[0].raw[i_v / bit_u] =
            (mask_1 << (i_v % bit_u)) | ct_verify->val[0].raw[i_v / bit_u];
      }
    }
  }
  // 对第二组操作
  for(int i_v = R_BITS; i_v < 2 * R_BITS; i_v++)
  {
    if(b[i_v] != 0)
    {
      b[i_v] = b[i_v] % 2;
      if(b[i_v] == 0)
      {
        // 用与操作
        ct_verify->val[1].raw[(i_v - R_BITS) / bit_u] =
            (mask_255 ^ (mask_1 << ((i_v - R_BITS) % bit_u))) &
            ct_verify->val[1].raw[(i_v - R_BITS) / bit_u];
      }
      else
      {
        // 用或操作
        ct_verify->val[1].raw[(i_v - R_BITS) / bit_u] =
            (mask_1 << ((i_v - R_BITS) % bit_u)) |
            ct_verify->val[1].raw[(i_v - R_BITS) / bit_u];
      }
    }
  }
}

// 将 增广常数数组 传递给 equations 的最后一列
_INLINE_ void
term_to_equations(OUT uint16_t         equations[][EQ_COLUMN],
                  IN const syndrome_t *pad_constant_term)
{
  // 处理前 11776 位
  for(uint8_t i = 0; i < R_QW - 1; i++)
  {
    for(uint64_t index = 0, location = 1; location != 0; location <<= 1)
    {
      if((pad_constant_term->qw[i] & location) != 0)
      {
        equations[64 * i + index][EQ_COLUMN - 1] = 1;
      }
      index++;
    }
  }
  // 处理最后三位
  for(uint64_t index = 0, location = 1; location <= MASK(LAST_R_QW_LEAD);
      location <<= 1)
  {
    if((pad_constant_term->qw[R_QW - 1] & location) != 0)
    {
      equations[64 * (R_QW - 1) + index][EQ_COLUMN - 1] = 1;
    }
    index++;
  }
}

// 方程组求解算法
_INLINE_ void
solving_equations(OUT uint16_t     *b,
                  IN uint16_t       equations[][EQ_COLUMN],
                  IN const uint16_t e_num)
{
  // 结果被保存在 b[23558] 中, 0 被保存为 2, 1 被保存为 1
  uint16_t M = e_num;
  uint16_t i, j;
  uint16_t y = 0;
  uint16_t t = 0;
  uint16_t c = 0;
  while(t < X)
  {
    for(i = 0; i < ROW; i++)
    {
      for(j = 0; j < X; j++)
      {
        if(equations[i][j] != 0)
          y++;
        else
          break;
      }
      if(y == 0)
      {
        equations[i][X] = 0;
        continue;
      } // 过滤全0
      for(j = 0; j < X; j++)
      {
        if(equations[i][j] != 0)
        {
          if(b[equations[i][j] - 1] != 0)
          {
            equations[i][X] = (equations[i][X] + b[equations[i][j] - 1]) %
                              2; // 直接减去索引的值
            equations[i][j] = 0;
          }
        }
      }
      if(y == 1)
      {
        for(j = 0; j < X; j++)
        {
          if(equations[i][j] != 0)
          {
            if(equations[i][X] == 0)
              b[equations[i][j] - 1] = 2;
            else
              b[equations[i][j] - 1] = 1;
            equations[i][j] = 0;
            equations[i][X] = 0;
          }
        }
      }
      y = 0;
    }
    for(j = 0; j < N; j++)
      if(b[j] != 0)
        c++;
    if(c == M)
      t += EQ_COLUMN;
    else
      t++;
  }
  // for(j = 0; j < N; j++)
  // {
  //   if(b[j] != 0)
  //   {
  //     printf("%d %d\n", j + 1, (b[j] % 2));
  //   }
  // }
}

// 对 qw[2 * R_QW] 循环右移一位, 仅 [185]-[369] 有效
_INLINE_ void
rotate_right_one(OUT single_h_t *out, IN const single_h_t *in)
{
  for(size_t i = 2 * R_QW - 1; i > 0; i--)
  {
    out->qw[i] = (in->qw[i] << 1) | (in->qw[i - 1] >> 63);
  }
  out->qw[0] = in->qw[0] << 1;
}

// 对 bytelen 长字节流, a 取反并和 b 与 (res = ~a & b)
_INLINE_ ret_t
negate_and(OUT uint8_t      *res,
           IN const uint8_t *a,
           IN const uint8_t *b,
           IN const uint64_t bytelen)
{
  for(uint64_t i = 0; i < bytelen; i++)
  {
    res[i] = (~a[i]) & b[i];
  }
  return SUCCESS;
}

// 对长字节流, a 和 b 或的值再和 res 或, 保存在 res 中 (res = (a | b) | res)
_INLINE_ ret_t
gf2x_or(OUT uint8_t      *res,
        IN const uint8_t *a,
        IN const uint8_t *b,
        IN const uint64_t bytelen)
{
  for(uint64_t i = 0; i < bytelen; i++)
  {
    res[i] = a[i] | b[i] | res[i];
  }
  return SUCCESS;
}

// 对 bytelen 长字节流, a 和 b '与', 索引存储在 res 行
// 注意：此处的索引和 wlist 不同的是从 1 开始
_INLINE_ ret_t
and_index(OUT uint16_t     *res,
          IN OUT uint8_t   *eq_index,
          IN const uint8_t *a,
          IN const uint8_t *b,
          IN const uint64_t bytelen,
          IN uint32_t       i_N0,
          IN uint16_t       i_eq)
{
  // tmp 用于暂存'与'信息
  uint8_t tmp[R_SIZE] = {0};

  // // count 用于记录列位置 (modify -> eq_index[i_eq])
  // uint8_t count = 0;

  // 定义低三位 mask_3 = 00000111 用于取值
  uint8_t mask_3 = MASK(R_BITS & MASK(3));
  // 对前 1472 个字节依次比对
  for(uint64_t i = 0; i < bytelen - 1; i++)
  {
    tmp[i] = a[i] & b[i];
    // 若存在重合则将重合索引添加到 res 中
    if(tmp[i] != 0)
    {
      // 将 location = 00000001 依次左移 和 tmp[i] 与运算找到重合位置
      for(uint8_t index = 1, location = 1; location != 0; location <<= 1)
      {
        if((location & tmp[i]) != 0)
        {
          res[eq_index[i_eq]] = i * 8 + index + i_N0 * R_BITS;
          (eq_index[i_eq])++;
        }
        index++;
      }
    }
  }
  // 对最后 3 位单独比对
  tmp[bytelen - 1] = a[bytelen - 1] & b[bytelen - 1] & mask_3;
  for(uint8_t index_2 = 1, location_2 = 1; location_2 <= mask_3; location_2 <<= 1)
  {
    if((location_2 & tmp[bytelen - 1]) != 0)
    {
      res[eq_index[i_eq]] = (bytelen - 1) * 8 + index_2 + i_N0 * R_BITS;
      (eq_index[i_eq])++;
    }
    index_2++;
  }
  return SUCCESS;
}

// Duplicates the first R_BITS of the syndrome three times
// |------------------------------------------|
// |  Third copy | Second copy | first R_BITS |
// |------------------------------------------|
// This is required by the rotate functions.
_INLINE_ void
dup(IN OUT syndrome_t *s)
{
  // R_QW = 185
  // qw[0] = (s63,s62,...,s1,s0)
  // qw[183] = (s11775,s11774,...,s11713,s11712)
  // 将最后一个 64 位 的空余位添加 s0 -- s60
  // qw[R_QW - 1] = (s60,s59,...,s0,s11778,s11777,s11776) = qw[184]
  s->qw[R_QW - 1] =
      (s->qw[0] << LAST_R_QW_LEAD) | (s->qw[R_QW - 1] & LAST_R_QW_MASK);

  // qw[R_QW] = (s124,s123,...,s63,s62,s61) = qw[185]
  // qw[R_QW + 1] = (s188,s187,...,s127,s126,s125) = qw[186]
  // ...
  // qw[R_QW + 182] = (s11772,s11771,...,s11711,s11710,s11709) = qw[367]
  // qw[R_QW + 183] = (s57,s56,...,s0,s11775,s11774,s11773) = qw[368]
  // qw[R_QW + 184] = (s121,s120,...,s60,s59,s58) = qw[369]
  // ...
  // qw[R_QW + 367] = (s54,s53,...,s0,s11775,s11774,s11773,s11772,s11771,s11770) =
  // qw[552] qw[R_QW + 368] = (s118,s117,...,s57,s56,s55) = qw[553]
  for(size_t i = 0; i < (2 * R_QW) - 1; i++)
  {
    s->qw[R_QW + i] =
        (s->qw[i] >> LAST_R_QW_TRAIL) | (s->qw[i + 1] << LAST_R_QW_LEAD);
  }
}

// Duplicates the first R_BITS of the syndrome two times
// |-----------------------------|
// |  first R_BITS | Second copy |
// |-----------------------------|
// 进入的 h 当从 185 位开始存放 (h0--h11779)，前 0-184 应为 0
// 此处减少一次复制，并为循环右移做准备
_INLINE_ void
dup_two(IN OUT single_h_t *h)
{
  // R_QW = 185
  // qw[369] = (0,0,...,h11778,h11777,h11776)
  // qw[185] = (h63,h62,...,h1,h0)
  // qw[184] = (h11778,h11777,...,h11716,h11715)
  if(LAST_R_QW_LEAD <= 32)
  {
    h->qw[0] = (h->qw[R_QW] << LAST_R_QW_TRAIL) |
               (h->qw[2 * R_QW - 1] << (LAST_R_QW_TRAIL_2 % 64)) |
               (h->qw[2 * R_QW - 2] >> (LAST_R_QW_LEAD_2 % 64));
  }
  else
  {
    h->qw[0] = (h->qw[R_QW] << LAST_R_QW_TRAIL) |
               (h->qw[2 * R_QW - 1] >> (LAST_R_QW_TRAIL_3 % 64));
  }

  // qw[0] = (h2,h1,h0,h11778,...,h11718)

  for(size_t i = R_QW - 1; i > 0; i--)
  {
    h->qw[i] = (h->qw[R_QW + i] << LAST_R_QW_TRAIL) |
               (h->qw[R_QW + i - 1] >> LAST_R_QW_LEAD);
  }
}

ret_t
compute_syndrome(OUT syndrome_t *syndrome, IN const ct_t *ct, IN const sk_t *sk)
{
  // gf2x_mod_mul requires the values to be 64bit padded and extra (dbl) space
  // for the results
  // gf2x_mod_mul 要求值是 64 位填充和额外 (dbl) 空间用于结果
  DEFER_CLEANUP(dbl_pad_syndrome_t pad_s, dbl_pad_syndrome_cleanup);
  DEFER_CLEANUP(pad_sk_t pad_sk = {0}, pad_sk_cleanup);
  pad_sk[0].val = sk->bin[0];
  pad_sk[1].val = sk->bin[1];

  DEFER_CLEANUP(pad_ct_t pad_ct = {0}, pad_ct_cleanup);
  pad_ct[0].val = ct->val[0];
  pad_ct[1].val = ct->val[1];

  // Compute s = c0*h0 + c1*h1:
  GUARD(gf2x_mod_mul((uint64_t *)&pad_s[0], (uint64_t *)&pad_ct[0],
                     (uint64_t *)&pad_sk[0]));
  GUARD(gf2x_mod_mul((uint64_t *)&pad_s[1], (uint64_t *)&pad_ct[1],
                     (uint64_t *)&pad_sk[1]));

  GUARD(gf2x_add(pad_s[0].val.raw, pad_s[0].val.raw, pad_s[1].val.raw, R_SIZE));

  // ---- test ---- 打印 pad_s 的值
  // for(uint16_t i_pad_s = 0; i_pad_s < 1473; i_pad_s++)
  // {
  //   printf("第 %u 个 s 的值为: %u\n", i_pad_s, pad_s[0].val.raw[i_pad_s]);
  // }

  // 将 s 按照每 8 位存储在 1473 个字节中
  // 复制 1473 个字节到 qw 的前 185 个 64 位整型中
  memcpy((uint8_t *)syndrome->qw, pad_s[0].val.raw, R_SIZE);

  // ---- test ---- 打印 syndrome->qw 中的值
  // for(uint16_t i_qw = 0; i_qw < 555; i_qw++)
  // {
  //   printf("第 %u 个 syndrome->qw 的值为: %lu\n", i_qw, syndrome->qw[i_qw]);
  // }

  dup(syndrome);
  // ---- test ---- 打印复制后 syndrome->qw 中的值
  // for(uint16_t i_qw = 0; i_qw < 555; i_qw++)
  // {
  //   printf("第 %u 个 syndrome->qw 的值为: %lu\n", i_qw, syndrome->qw[i_qw]);
  // }

  return SUCCESS;
}

_INLINE_ ret_t
recompute_syndrome(OUT syndrome_t     *syndrome,
                   IN const ct_t      *ct,
                   IN const sk_t      *sk,
                   IN const split_e_t *splitted_e)
{
  ct_t tmp_ct = *ct;

  // Adapt the ciphertext
  // 更新 c+e
  GUARD(gf2x_add(tmp_ct.val[0].raw, tmp_ct.val[0].raw, splitted_e->val[0].raw,
                 R_SIZE));
  GUARD(gf2x_add(tmp_ct.val[1].raw, tmp_ct.val[1].raw, splitted_e->val[1].raw,
                 R_SIZE));

  // Recompute the syndrome
  // 计算更新后的 c 与 H 计算 s
  GUARD(compute_syndrome(syndrome, &tmp_ct, sk));

  return SUCCESS;
}

_INLINE_ uint8_t
get_threshold(IN const syndrome_t *s)
{
  bike_static_assert(sizeof(*s) >= sizeof(r_t), syndrome_is_large_enough);

  const uint32_t syndrome_weight = r_bits_vector_weight((const r_t *)s->qw);

  // The equations below are defined in BIKE's specification:
  // https://bikesuite.org/files/round2/spec/BIKE-Spec-Round2.2019.03.30.pdf
  // Page 20 Section 2.4.2
  const uint8_t threshold =
      THRESHOLD_COEFF0 + (THRESHOLD_COEFF1 * syndrome_weight);

  DMSG("    Thresold: %d\n", threshold);

  if(threshold > THRESHOLD_MAX)
  {
    return threshold;
  }
  else
  {
    return THRESHOLD_MAX;
  }
}

// Use half-adder as described in [5].
_INLINE_ void
bit_sliced_adder(OUT upc_t         *upc,
                 IN OUT syndrome_t *rotated_syndrome,
                 IN const size_t    num_of_slices)
{
  // From cache-memory perspective this loop should be the outside loop
  for(size_t j = 0; j < num_of_slices; j++)
  {
    for(size_t i = 0; i < R_QW; i++)
    {
      const uint64_t carry = (upc->slice[j].u.qw[i] & rotated_syndrome->qw[i]);
      upc->slice[j].u.qw[i] ^= rotated_syndrome->qw[i];
      rotated_syndrome->qw[i] = carry;
    }
  }
}

_INLINE_ void
bit_slice_full_subtract(OUT upc_t *upc, IN uint8_t val)
{
  // Borrow
  uint64_t br[R_QW] = {0};

  for(size_t j = 0; j < SLICES; j++)
  {

    const uint64_t lsb_mask = 0 - (val & 0x1);
    val >>= 1;

    // Perform a - b with c as the input/output carry
    // br = 0 0 0 0 1 1 1 1
    // a  = 0 0 1 1 0 0 1 1
    // b  = 0 1 0 1 0 1 0 1
    // -------------------
    // o  = 0 1 1 0 0 1 1 1
    // c  = 0 1 0 0 1 1 0 1
    //
    // o  = a^b^c
    //            _     __    _ _   _ _     _
    // br = abc + abc + abc + abc = abc + ((a+b))c

    for(size_t i = 0; i < R_QW; i++)
    {
      const uint64_t a      = upc->slice[j].u.qw[i];
      const uint64_t b      = lsb_mask;
      const uint64_t tmp    = ((~a) & b & (~br[i])) | ((((~a) | b) & br[i]));
      upc->slice[j].u.qw[i] = a ^ b ^ br[i];
      br[i]                 = tmp;
    }
  }
}

// Calculate the Unsatisfied Parity Checks (UPCs) and update the errors
// vector (e) accordingy. In addition, update the black and gray errors vector
// with the relevant values.
// 计算未满足的奇偶校验 (UPC) 并相应地更新错误向量 (e)。
// 此外，用相关值更新黑色和灰色误差向量。 此部分对应 procedure BitFlipIter(s, e,
// th, H)
_INLINE_ void
find_err1(OUT split_e_t                  *e,
          OUT split_e_t                  *black_e,
          OUT split_e_t                  *gray_e,
          IN const syndrome_t            *syndrome,
          IN const compressed_idx_dv_ar_t wlist,
          IN const uint8_t                threshold,
          IN const uint8_t                delat)
{
  // This function uses the bit-slice-adder methodology of [5]:
  // QcBits: Constant-Time Small-Key Code-Based Cryptography
  // 此函数使用 [5] 中的 bit-slice-adder 方法：
  DEFER_CLEANUP(syndrome_t rotated_syndrome = {0}, syndrome_cleanup);
  DEFER_CLEANUP(upc_t upc, upc_cleanup);

  for(uint32_t i = 0; i < N0; i++)
  {
    // UPC must start from zero at every iteration
    memset(&upc, 0, sizeof(upc));

    // 1) Right-rotate the syndrome for every secret key set bit index
    //    Then slice-add it to the UPC array.
    // 对每个密钥集位索引的校正子进行右循环
    // 然后将其切片添加到 UPC 数组中
    for(size_t j = 0; j < DV; j++)
    {
      // 向右旋转 syndrome 的第一个 R_BITS
      // 假设：syndrome 包含三个 R_BITS 重复
      // 输出校验子仅包含一个 R_BITS 旋转，其他 (2 * R_BITS) 位未定义
      rotate_right(&rotated_syndrome, syndrome, wlist[i].val[j]);
      bit_sliced_adder(&upc, &rotated_syndrome, LOG2_MSB(j + 1));
    }

    // 2) Subtract the threshold from the UPC counters
    // 从 UPC 计数器中减去阈值
    bit_slice_full_subtract(&upc, threshold);

    // 3) Update the errors and the black errors vectors.
    //    The last slice of the UPC array holds the MSB of the accumulated values
    //    minus the threshold. Every zero bit indicates a potential error bit.
    //    The errors values are stored in the black array and xored with the
    //    errors Of the previous iteration.
    // 3) 更新错误和黑色错误向量。
    // UPC 数组的最后一个切片保存累积值的 MSB 减去阈值。
    // 每个零位表示一个潜在的错误位。
    // 错误值存储在黑色数组中，并与上一次迭代的错误进行异或
    const r_t *last_slice = &(upc.slice[SLICES - 1].u.r.val);
    for(size_t j = 0; j < R_SIZE; j++)
    {
      const uint8_t sum_msb  = (~last_slice->raw[j]);
      black_e->val[i].raw[j] = sum_msb;
      e->val[i].raw[j] ^= sum_msb;
    }

    // Ensure that the padding bits (upper bits of the last byte) are zero so
    // they will not be included in the multiplication and in the hash function.
    // 确保填充位（最后一个字节的高位）为零，因此它们不会包含在乘法和散列函数中。
    e->val[i].raw[R_SIZE - 1] &= LAST_R_BYTE_MASK;

    // 4) Calculate the gray error array by adding "DELTA" to the UPC array.
    //    For that we reuse the rotated_syndrome variable setting it to all "1".
    // 通过将 “DELTA” δ 添加到 UPC 数组来计算灰度误差数组。
    // 为此，我们重用 rotate_syndrome 变量，将其设置为全“1”。
    for(size_t l = 0; l < delat; l++)
    {
      memset((uint8_t *)rotated_syndrome.qw, 0xff, R_SIZE);
      bit_sliced_adder(&upc, &rotated_syndrome, SLICES);
    }

    // 5) Update the gray list with the relevant bits that are not
    //    set in the black list.
    // 用黑-名单中没有设置的相关位更新灰-名单。
    for(size_t j = 0; j < R_SIZE; j++)
    {
      const uint8_t sum_msb = (~last_slice->raw[j]);
      gray_e->val[i].raw[j] = (~(black_e->val[i].raw[j])) & sum_msb;
    }
  }
}

// Recalculate the UPCs and update the errors vector (e) according to it
// and to the black/gray vectors.
// 重新计算 UPC 并根据它和黑色/灰色向量组更新错误向量 (e)。
_INLINE_ void
find_err2(OUT split_e_t                  *e,
          IN split_e_t                   *pos_e,
          IN const syndrome_t            *syndrome,
          IN const compressed_idx_dv_ar_t wlist,
          IN const uint8_t                threshold)
{
  DEFER_CLEANUP(syndrome_t rotated_syndrome = {0}, syndrome_cleanup);
  DEFER_CLEANUP(upc_t upc, upc_cleanup);

  for(uint32_t i = 0; i < N0; i++)
  {
    // UPC must start from zero at every iteration
    memset(&upc, 0, sizeof(upc));

    // 1) Right-rotate the syndrome for every secret key set bit index
    //    Then slice-add it to the UPC array.
    // 对每个密钥集位索引的校正子进行右循环
    // 然后将其切片添加到 UPC 数组中。
    for(size_t j = 0; j < DV; j++)
    {
      rotate_right(&rotated_syndrome, syndrome, wlist[i].val[j]);
      bit_sliced_adder(&upc, &rotated_syndrome, LOG2_MSB(j + 1));
    }

    // 2) Subtract the threshold from the UPC counters
    // 从 UPC 计数器中减去阈值
    bit_slice_full_subtract(&upc, threshold);

    // 3) Update the errors vector.
    //    The last slice of the UPC array holds the MSB of the accumulated values
    //    minus the threshold. Every zero bit indicates a potential error bit.
    // 更新错误向量。
    // UPC 数组的最后一个切片保存累积值减去阈值的 MSB。
    // 每个零位表示一个潜在的错误位。
    const r_t *last_slice = &(upc.slice[SLICES - 1].u.r.val);
    for(size_t j = 0; j < R_SIZE; j++)
    {
      const uint8_t sum_msb = (~last_slice->raw[j]);
      e->val[i].raw[j] ^= (pos_e->val[i].raw[j] & sum_msb);
    }

    // Ensure that the padding bits (upper bits of the last byte) are zero so
    // they will not be included in the multiplication and in the hash function.
    // 确保填充位（最后一个字节的高位）为零，因此它们不会包含在乘法和散列函数中。
    e->val[i].raw[R_SIZE - 1] &= LAST_R_BYTE_MASK;
  }
}

// 此译码算法依据 QC-MDPC decoders with several shades of gray 中第 4 页
ret_t
decode(OUT split_e_t       *black_or_gray_e_out,
       OUT split_e_t       *e,
       IN OUT uint8_t      *flag,
       IN const split_e_t  *R_e,
       IN const syndrome_t *original_s,
       IN const ct_t       *ct,
       IN const sk_t       *sk,
       IN const uint8_t     delat,
       IN OUT uint32_t     *decoder_error_count,
       IN OUT uint32_t     *equations_error_count)
{
  // 初始化黑灰数组
  split_e_t  black_e           = {0};
  split_e_t  gray_e            = {0};
  split_e_t  black_e_eq        = {0};
  split_e_t  gray_e_eq         = {0};
  split_e_t  e_eq              = {0};
  split_e_t  black_or_gray_e   = {0};
  split_e_t  fixed_e           = {0};
  ct_t       ct_remove_BG      = {0};
  ct_t       ct_pad            = {0};
  ct_t       ct_verify         = {0};
  sk_t       sk_transpose      = {0};
  syndrome_t pad_constant_term = {0};
  syndrome_t s;
  uint8_t    delat_eq        = DELTA_EQ;
  uint8_t    delat_eq_step23 = DELTA_STEP23;

  // 定义 11779 行方程组, 前 EQ_COLUMN-1 个元素用于保存索引, 第 EQ_COLUMN
  // 个用于存放增广常数
  uint16_t equations[R_BITS][EQ_COLUMN] = {0};

  // 构建出循环矩阵的索引 h_matrix 方便后面使用
  uint16_t sk_wlist_all_0[R_BITS][DV] = {0};
  uint16_t sk_wlist_all_1[R_BITS][DV] = {0};

  // 填充对应的索引值
  // 填充第一行
  for(uint16_t i_DV = 0; i_DV < DV; i_DV++)
  {
    sk_wlist_all_0[0][i_DV] = sk->wlist[0].val[i_DV];
    sk_wlist_all_1[0][i_DV] = sk->wlist[1].val[i_DV];
  }

  // 填充后 2 - 11779 行
  for(uint16_t i_r = 1; i_r < R_BITS; i_r++)
  {
    for(uint16_t i_DV = 0; i_DV < DV; i_DV++)
    {
      sk_wlist_all_0[i_r][i_DV] = (sk_wlist_all_0[i_r - 1][i_DV] + 1) % R_BITS;
      sk_wlist_all_1[i_r][i_DV] = (sk_wlist_all_1[i_r - 1][i_DV] + 1) % R_BITS;
    }
  }

  // 初始化 fixed_e 为 R_e
  fixed_e.val[0] = R_e->val[0];
  fixed_e.val[1] = R_e->val[1];

  // Reset (init) the error because it is xored in the find_err funcitons.
  // 初始化 e
  memset(e, 0, sizeof(*e));
  s = *original_s;
  dup(&s);

  // 进入大迭代过程(for itr in 1...XBG do:)
  for(uint32_t iter = 0; iter < MAX_IT; iter++)
  {
    // 将 fixed_e 和 求出来的 e 异或
    GUARD(gf2x_add((uint8_t *)&fixed_e.val[0].raw, R_e->val[0].raw, e->val[0].raw,
                   R_SIZE));
    GUARD(gf2x_add((uint8_t *)&fixed_e.val[1].raw, R_e->val[1].raw, e->val[1].raw,
                   R_SIZE));

    // 解码器使用阈值(th)来决定某个位是否为错误位
    // 该位确是错误位的概率随着间隙(upc[i] - th)的增加而增加
    // 该算法记录黑/灰掩码中有小间隙的位，以便后续步骤II和步骤III可以使用掩码，以获得翻转位的更多信息
    // 22: th = computeThreshold(s)
    // 参: Bit Flipping Key Encapsulation(v2.1) 17页，Threshold Selection Rule
    // printf("\n---->当前迭代阶段: %d<----\n", iter);

    // 获取当前 fixed_e 的重量
    uint16_t fixed_e_weight = r_bits_vector_weight(&fixed_e.val[0]) +
                              r_bits_vector_weight(&fixed_e.val[1]);

    uint8_t threshold = 0;
    // 选择使用何种方法计算 th
    if(TH_SELECT == 0)
    {
      threshold = get_threshold(&s);
    }
    else
    {
      threshold = compute_th_R(sk_wlist_all_0, sk_wlist_all_1, fixed_e_weight,
                               &fixed_e, &s);
    }

    // ---- test ---- 查看 th
    // printf("threshold = %u\n", threshold);
    // printf("MY_threshold = %u\n\n", threshold_2);

    DMSG("    Iteration: %d\n", iter);
    DMSG("    Weight of e: %lu\n",
         r_bits_vector_weight(&e->val[0]) + r_bits_vector_weight(&e->val[1]));
    DMSG("    Weight of syndrome: %lu\n", r_bits_vector_weight((r_t *)s.qw));

    // 进入选取解方程黑灰集合
    find_err1(&e_eq, &black_e_eq, &gray_e_eq, &s, sk->wlist, threshold, delat_eq);

    // 23:  (s, e, black, gray) = BitFlipIter(s, e, th, H) . Step I
    // H -- sk->wlist
    // 进入 procedure BitFlipIter(s, e, th, H)
    find_err1(e, &black_e, &gray_e, &s, sk->wlist, threshold, delat);

    for(uint8_t i = 0; i < N0; i++)
    {
      if(iter <= 0)
      {
        // 将黑灰集合'或'运算(black_e | gray_e) 存放于
        // black_or_gray_e，即所有未知数位
        GUARD(gf2x_or((uint8_t *)&black_or_gray_e.val[i].raw,
                      black_e_eq.val[i].raw, gray_e_eq.val[i].raw, R_SIZE));
      }
      else
      {
        GUARD(gf2x_or((uint8_t *)&black_or_gray_e.val[i].raw,
                      black_e_eq.val[i].raw, gray_e_eq.val[i].raw, R_SIZE));
      }
    }

    // 10:  s = H(cT + eT ) . 更新校验子 syndrome
    GUARD(recompute_syndrome(&s, ct, sk, e));

// 此处代码中在 iter >= 1 时候去除了 Step II 和 Step III (BGF)
// 相当于只进行了一轮黑灰迭代后进行了多轮比特位反转(step I)
#ifdef BGF_DECODER
    if(iter >= 1)
    {
      continue;
    }
#endif
    DMSG("    Weight of e: %lu\n",
         r_bits_vector_weight(&e->val[0]) + r_bits_vector_weight(&e->val[1]));
    DMSG("    Weight of syndrome: %lu\n", r_bits_vector_weight((r_t *)s.qw));

    if(STEP23_ON == 1)
    {
      // 进入选取解方程黑灰集合
      find_err1(&e_eq, &black_e_eq, &gray_e_eq, &s, sk->wlist, ((DV + 1) / 2) + 1,
                delat_eq_step23);

      for(uint8_t i = 0; i < N0; i++)
      {
        // 将黑灰集合'或'运算(black_e | gray_e) 存放于
        // black_or_gray_e，即所有未知数位
        GUARD(gf2x_or((uint8_t *)&black_or_gray_e.val[i].raw,
                      black_e_eq.val[i].raw, gray_e_eq.val[i].raw, R_SIZE));
      }
    }

    // 24:  (s, e) = BitFlipMaskedIter(s, e, black, ((d + 1)/2), H) . Step II
    // procedure BitFlipMaskedIter(s, e, mask, th, H)
    find_err2(e, &black_e, &s, sk->wlist, ((DV + 1) / 2) + 1);
    GUARD(recompute_syndrome(&s, ct, sk, e));

    DMSG("    Weight of e: %lu\n",
         r_bits_vector_weight(&e->val[0]) + r_bits_vector_weight(&e->val[1]));
    DMSG("    Weight of syndrome: %lu\n", r_bits_vector_weight((r_t *)s.qw));

    // 25:  (s, e) = BitFlipMaskedIter(s, e, gray, ((d + 1)/2), H) . Step III
    // procedure BitFlipMaskedIter(s, e, mask, th, H)
    find_err2(e, &gray_e, &s, sk->wlist, ((DV + 1) / 2) + 1);

    GUARD(recompute_syndrome(&s, ct, sk, e));
  }

  // ================================================================

  // // 计算求解的 未知数 总个数(black_or_gray_e 的重量)
  // uint16_t x_weight_1 = r_bits_vector_weight((r_t *)black_or_gray_e.val[0].raw)
  // +
  //                     r_bits_vector_weight((r_t *)black_or_gray_e.val[1].raw);
  // printf("x_weight_1: %u\n",x_weight_1);

  // // 将未知数个数填充满
  // if(x_weight_1 < R_BITS)
  // {
  //   uint16_t pad_num_1 = (R_BITS - x_weight_1) / 2;
  //   uint16_t pad_num_2 = (R_BITS - x_weight_1 - pad_num_1);

  //   for(uint16_t test_i = 0; test_i < R_SIZE; test_i++)
  //   {
  //     for(uint8_t mask_1 = 1; mask_1 != 0; mask_1 <<= 1)
  //     {
  //       if((mask_1 & black_or_gray_e.val[0].raw[test_i]) == 0)
  //       {
  //         black_or_gray_e.val[0].raw[test_i] += mask_1;
  //         pad_num_1 -= 1;
  //       }
  //       if(pad_num_1 == 0)
  //       {
  //         break;
  //       }
  //     }
  //     if(pad_num_1 == 0)
  //     {
  //       break;
  //     }
  //   }
  //   for(uint16_t test_i = 0; test_i < R_SIZE; test_i++)
  //   {
  //     for(uint8_t mask_1 = 1; mask_1 != 0; mask_1 <<= 1)
  //     {
  //       if((mask_1 & black_or_gray_e.val[1].raw[test_i]) == 0)
  //       {
  //         black_or_gray_e.val[1].raw[test_i] += mask_1;
  //         pad_num_2 -= 1;
  //       }
  //       if(pad_num_2 == 0)
  //       {
  //         break;
  //       }
  //     }
  //     if(pad_num_2 == 0)
  //     {
  //       break;
  //     }
  //   }
  // }

  // 计算求解的 未知数 总个数(black_or_gray_e 的重量)
  uint16_t x_weight = r_bits_vector_weight((r_t *)black_or_gray_e.val[0].raw) +
                      r_bits_vector_weight((r_t *)black_or_gray_e.val[1].raw);

  // 将 black_or_gray_e 传递出去比较是否包含所有错误向量
  for(uint16_t i = 0; i < N0; i++)
  {
    black_or_gray_e_out->val[i] = black_or_gray_e.val[i];
  }

  // 创建失败 flag 1 成功 0 失败
  uint16_t flag_BG = 1;

  //  26: if (wt(s) != 0) then
  //  27:     return ⊥(ERROR)
  if(r_bits_vector_weight((r_t *)s.qw) > 0)
  {
    FILE *fp_3;
    fp_3 = fopen("weight_bad.txt", "a");
    fprintf(fp_3, "DELAT: %d 黑灰译码失败\n", delat);
    fclose(fp_3);
    *flag = 1;
    DMSG("s 重量不为 0...");
    flag_BG = 0;
    *decoder_error_count += 1;
  }

  if(flag_BG == 0)
  {

    // ================> 增加方程组求解算法(当 s 不为 0) <================
    // =================================================================
    // --------------------- 1.构建方程组 ---------------------

    // 构建方程组计时
    double start_construct = clock();

    for(uint32_t i = 0; i < N0; i++)
    {
      // 获取 ct 的值
      ct_pad.val[i] = ct->val[i];

      // 构造 sk 转置 sk_transpose
      // 获取 sk 转置的首行索引
      // 𝜑(A)' = a0 + ar-1X + ar-2X^2 ...
      for(uint8_t i_DV = 0; i_DV < DV; i_DV++)
      {
        if(sk->wlist[i].val[i_DV] != 0)
        {
          sk_transpose.wlist[i].val[i_DV] = R_BITS - sk->wlist[i].val[i_DV];
        }
        else
        {
          sk_transpose.wlist[i].val[i_DV] = sk->wlist[i].val[i_DV];
        }
      }

      // 去除 c 中的未知数位，将 black_or_gray_e 取反后与 c 做与操作
      GUARD(negate_and(ct_remove_BG.val[i].raw, black_or_gray_e.val[i].raw,
                       ct_pad.val[i].raw, R_SIZE));
    }

    // 对方程组未知数进行构建，将 x0-xall 的对应关系列出来
    // black_or_gray_e 的每个位置对应 旋转 h 的位置满足 (e+r-h) % r
    // 对每个 black_or_gray_e 进行 and 寻找是否存在未知数
    // guss_j_num 最后一个字用来存储 b
    uint16_t guss_j_num = 0;
    if(x_weight % GUSS_BLOCK == 0)
    {
      guss_j_num = x_weight / GUSS_BLOCK + 1;
    }
    else
    {
      guss_j_num = x_weight / GUSS_BLOCK + 2;
    }
    uint8_t equations_guss_byte[R_BITS][guss_j_num];
    memset(equations_guss_byte, 0, sizeof(equations_guss_byte));

    uint8_t  mask_e       = 1;
    uint8_t  mask_e_byte  = 1;
    uint16_t e_count      = 0;
    uint16_t e_index      = 0;
    uint16_t e_index_byte = 0;
    uint16_t x_arr[x_weight];
    memset(x_arr, 0, sizeof(x_arr));
    for(uint8_t i = 0; i < N0; i++)
    {
      for(uint16_t i_e_x = 0; i_e_x < R_BITS; i_e_x++)
      {
        if(i_e_x % GUSS_BLOCK == 0)
        {
          mask_e  = 1;
          e_index = i_e_x / GUSS_BLOCK;
        }
        if((mask_e & black_or_gray_e.val[i].raw[e_index]) != 0)
        {
          if(e_count % GUSS_BLOCK == 0)
          {
            mask_e_byte  = 1;
            e_index_byte = e_count / GUSS_BLOCK;
          }
          uint16_t e_add_R = i_e_x + R_BITS;
          x_arr[e_count]   = i_e_x + i * R_BITS;
          e_count += 1;
          // 根据 e 的和 h 的位置来确定 equations_guss_byte 的构建 (e+r-h) % r
          for(uint16_t wlist_i = 0; wlist_i < DV; wlist_i++)
          {
            equations_guss_byte[(e_add_R - sk_transpose.wlist[i].val[wlist_i]) %
                                R_BITS][e_index_byte] += mask_e_byte;
          }
          mask_e_byte <<= 1;
        }
        mask_e <<= 1;
      }
    }

    // 将 ct_remove_BG 和 H 相乘, 使用 gf2x_mod_mul(), 得到结果 constant_term
    // 这里计算方式与 compute_syndrome() 计算方式一致, 可调用此函数构建
    GUARD(compute_syndrome(&pad_constant_term, &ct_remove_BG, sk));

    // ---- test ---- 打印 pad_constant_term 的值
    print("\npad_constant_term: \n", (uint64_t *)pad_constant_term.qw, R_BITS);

    // 将增广常数 pad_constant_term 赋值给 equations[i][EQ_COLUMN]
    term_to_equations(equations, (syndrome_t *)&pad_constant_term);

    // equations_guss_byte 加入常数列
    for(uint16_t i_b = 0; i_b < R_BITS; i_b++)
    {
      if(equations[i_b][EQ_COLUMN - 1] == 1)
      {
        equations_guss_byte[i_b][guss_j_num - 1] = 1;
      }
    }

    double end_construct = clock();
    printf("建立方程 took %lfs\n",
           ((double)(end_construct - start_construct) / CLOCKS_PER_SEC));

    // ==========================================================

    // 计时
    double start = clock();
    // 设置 x 主元表
    uint8_t guss_x_main[R_BITS] = {0};
    // 开始消元
    for(uint16_t guss_j = 0; guss_j < x_weight; guss_j++)
    {
      uint8_t  mask_1    = 1;
      uint8_t  mask_guss = (mask_1 << (guss_j % GUSS_BLOCK));
      uint16_t eq_j      = guss_j / GUSS_BLOCK;
      for(uint16_t guss_i = guss_j; guss_i < R_BITS; guss_i++)
      {
        if((mask_guss & equations_guss_byte[guss_i][eq_j]) != 0)
        {
          if(guss_x_main[guss_j] == 0)
          {
            // 如果此列没有主元优先挑选主元
            // 将此行作为当前列主元，交换第一行并继续向后消元
            guss_x_main[guss_j] = 1;
            swap(equations_guss_byte[guss_j], equations_guss_byte[guss_i], eq_j,
                 guss_j_num);
          }
          else
          {
            // 使用第 guss_j 行消此行
            GUARD(xor_8(equations_guss_byte[guss_i], equations_guss_byte[guss_i],
                        equations_guss_byte[guss_j], guss_j_num, eq_j));
          }
        }
      }
    }
    // 倒着求解
    for(int guss_j = x_weight - 1; guss_j >= 0; guss_j--)
    {
      uint16_t eq_j = guss_j / GUSS_BLOCK;
      for(uint16_t guss_i = guss_j; guss_i > 0; guss_i--)
      {
        if((equations_guss_byte[guss_j][eq_j] &
            equations_guss_byte[guss_i - 1][eq_j]) != 0)
        {
          equations_guss_byte[guss_i - 1][eq_j] ^=
              equations_guss_byte[guss_j][eq_j];
          equations_guss_byte[guss_i - 1][guss_j_num - 1] ^=
              equations_guss_byte[guss_j][guss_j_num - 1];
        }
      }
    }
    double end = clock();
    printf("guss took %lfs\n", ((double)(end - start) / CLOCKS_PER_SEC));

    // ========================================================================

    // 构造高斯消元解数组
    uint16_t b[N] = {0};
    for(uint16_t i = 0; i < x_weight; i++)
    {
      if(equations_guss_byte[i][guss_j_num - 1] == 0)
      {
        b[x_arr[i]] = 2;
      }
      else
      {
        b[x_arr[i]] = 1;
      }
    }

    // // 解方程计时
    // double start_solve = clock();
    // // --------------------- 2.解方程函数 ---------------------
    // 结果被保存在 b[23558] 中, 0 被保存为 2, 1 被保存为 1
    // uint16_t b[N] = {0};
    // solving_equations((uint16_t *)&b, equations, x_weight);

    // 检验解方程的正确性, 将 ct 对应位置放上解方程结果 b, 还原 fm 加真实 e 和 ct
    // 比较
    ct_verify.val[0] = ct->val[0];
    ct_verify.val[1] = ct->val[1];
    solving_equations_mf((ct_t *)&ct_verify, b);

    // double end_solve = clock();
    // printf("-- 解方程时间 --: %lfs\n",
    //        ((double)(end_solve - start_solve) / CLOCKS_PER_SEC));

    // 将 ct_verify = mf 和真实 e 异或后再异或 ct 检查重量
    uint8_t verify_weight[2] = {0};
    for(uint8_t i = 0; i < N0; i++)
    {
      GUARD(gf2x_add((uint8_t *)&ct_verify.val[i].raw, ct_verify.val[i].raw,
                     R_e->val[i].raw, R_SIZE));
      GUARD(gf2x_add((uint8_t *)&ct_verify.val[i].raw, ct_verify.val[i].raw,
                     ct->val[i].raw, R_SIZE));
      verify_weight[i] = r_bits_vector_weight((r_t *)ct_verify.val[i].raw);
    }
    print("ct_verify.val[0]: ", (uint64_t *)ct_verify.val[0].raw, R_BITS);
    print("ct_verify.val[1]: ", (uint64_t *)ct_verify.val[1].raw, R_BITS);

    if(verify_weight[0] || verify_weight[1] != 0)
    {
      printf("DELAT: %d 需求解未知数: %u 解方程失败\n", delat, x_weight);
      FILE *fp_2;
      fp_2 = fopen("weight_bad.txt", "a");
      fprintf(fp_2, "DELAT: %d 当前未知数: %u 解方程失败\n", delat, x_weight);
      fclose(fp_2);
      *flag = 1;
      *equations_error_count += 1;
    }
    else
    {
      printf("DELAT: %d 需求解未知数: %u\n", delat, x_weight);
    }

    if(flag_BG == 0)
    {
      BIKE_ERROR(E_DECODING_FAILURE);
    }
  }

  return SUCCESS;
}
