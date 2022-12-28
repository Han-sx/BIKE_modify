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
#    define MAX_IT 5
#  elif(LEVEL == 5)
#    define MAX_IT 5
#  else
#    error "Level can only be 1/3/5"
#  endif
#endif

// 0 使用代码原始线性拟合计算, 1 使用论文方法计算
#define TH_SELECT 0

#define EQ_COLUMN 101 // 索引矩阵列数
#define ROW       R_BITS
#define X         EQ_COLUMN - 1
#define N         2 * R_BITS

// 用于计算出upc切片的值并保存在文件中
_INLINE_ void
compute_upc_and_save_test(IN upc_t upc)
{
  // ---- test ---- 将 upc 切片的值计算出来并保存
  uint64_t mask_1 = 1;
  // 处理前 R_QW - 1 位
  // 将每层累计得到的 upc_i 写入文件
  FILE *fp_2;
  fp_2 = fopen("iter_data_all.txt", "a");
  for(uint16_t i_upc = 0; i_upc < R_QW - 1; i_upc++)
  {
    for(uint64_t location = 1; location != 0; location <<= 1)
    {
      // 用于保存每个upc[i]的值
      uint16_t upc_i = 0;
      for(uint8_t location_s = 1, i_upc_s = 0; i_upc_s < SLICES - 1;
          i_upc_s++, location_s <<= 1)
      {
        if((upc.slice[i_upc_s].u.qw[i_upc] & location) != 0)
        {
          upc_i += location_s;
        }
      }
      fprintf(fp_2, "%u ", upc_i);
    }
  }
  // 处理最后 R_BITS - (R_QW - 1) * 64 位
  for(uint64_t location = 1; location < (mask_1 << (R_BITS - (R_QW - 1) * 64));
      location <<= 1)
  {
    // 用于保存每个upc[i]的值
    uint16_t upc_i = 0;
    for(uint8_t location_s = 1, i_upc_s = 0; i_upc_s < SLICES - 1;
        i_upc_s++, location_s <<= 1)
    {
      if((upc.slice[i_upc_s].u.qw[R_QW - 1] & location) != 0)
      {
        upc_i += location_s;
      }
    }
    fprintf(fp_2, "%u ", upc_i);
  }
  fclose(fp_2);
}

// 用于计算出upc切片的值并保存在文件中
_INLINE_ void
compute_upc_and_save(IN upc_t upc)
{
  // ---- test ---- 将 upc 切片的值计算出来并保存
  uint64_t mask_1 = 1;
  // 处理前 R_QW - 1 位
  // 将每层累计得到的 upc_i 写入文件
  FILE *fp_2;
  fp_2 = fopen("iter_data.txt", "a");
  for(uint16_t i_upc = 0; i_upc < R_QW - 1; i_upc++)
  {
    for(uint64_t location = 1; location != 0; location <<= 1)
    {
      // 用于保存每个upc[i]的值
      uint16_t upc_i = 0;
      for(uint8_t location_s = 1, i_upc_s = 0; i_upc_s < SLICES - 1;
          i_upc_s++, location_s <<= 1)
      {
        if((upc.slice[i_upc_s].u.qw[i_upc] & location) != 0)
        {
          upc_i += location_s;
        }
      }
      fprintf(fp_2, "%u ", upc_i);
    }
  }
  // 处理最后 R_BITS - (R_QW - 1) * 64 位
  for(uint64_t location = 1; location < (mask_1 << (R_BITS - (R_QW - 1) * 64));
      location <<= 1)
  {
    // 用于保存每个upc[i]的值
    uint16_t upc_i = 0;
    for(uint8_t location_s = 1, i_upc_s = 0; i_upc_s < SLICES - 1;
        i_upc_s++, location_s <<= 1)
    {
      if((upc.slice[i_upc_s].u.qw[R_QW - 1] & location) != 0)
      {
        upc_i += location_s;
      }
    }
    fprintf(fp_2, "%u ", upc_i);
  }
  fclose(fp_2);
}

// 利用论文中的方法计算 th
_INLINE_ double_t
compute_th_R(IN OUT double_t *x_out,
             IN OUT double_t *p_p,
             IN uint16_t      sk_wlist_all_0[][DV],
             IN uint16_t      sk_wlist_all_1[][DV],
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
      // 求一个 概率验证数 p_p
      *p_p += (double_t)pow(10, (A + B - C));
    }
  }

  // ---- 3. 计算 T ----
  double_t pai_0 = (double_t)(2 * DV * s_weight - x) / ((2 * R_BITS - T) * DV);
  double_t pai_1 = (double_t)(s_weight + x) / (T * DV);
  // 取出 pai_1
  *x_out = x;

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
  return th;
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
  for(uint32_t i = 0; i < R_QW - 1; i++)
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
  for(uint64_t index = 0, location = 1; location < 8; location <<= 1)
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
  uint32_t M = e_num;
  uint32_t i, j;
  uint32_t y = 0;
  uint32_t t = 0;
  uint32_t c = 0;
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
      t += 100;
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
  for(size_t i = 369; i > 0; i--)
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
  uint8_t mask_3 = 7;
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
  for(uint8_t index_2 = 1, location_2 = 1; location_2 < 8; location_2 <<= 1)
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
// _INLINE_ void
// dup_two(IN OUT single_h_t *h)
// {
//   // R_QW = 185
//   // qw[369] = (0,0,...,h11778,h11777,h11776)
//   // qw[185] = (h63,h62,...,h1,h0)
//   // qw[184] = (h11778,h11777,...,h11716,h11715)
//   h->qw[0] = (h->qw[185] << LAST_R_QW_TRAIL) | (h->qw[369] <<
//   LAST_R_QW_TRAIL_2) |
//              (h->qw[368] >> LAST_R_QW_LEAD_2);
//   // qw[0] = (h2,h1,h0,h11778,...,h11718)

//   for(size_t i = 184; i > 0; i--)
//   {
//     h->qw[i] = (h->qw[R_QW + i] << LAST_R_QW_TRAIL) |
//                (h->qw[R_QW + i - 1] >> LAST_R_QW_LEAD);
//   }
// }

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
find_err1(OUT upc_all_t                  *upc_out,
          OUT split_e_t                  *e,
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

    // // ---- test ---- 将 upc 切片的值计算出来并保存
    // compute_upc_and_save(upc);

    // 保存 upc 到 upc_out
    for(uint8_t slice_i = 0; slice_i < SLICES; slice_i++)
    {
      upc_out->val[i].slice[slice_i] = upc.slice[slice_i];
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
find_err2(OUT upc_all_t                  *upc_out,
          OUT split_e_t                  *e,
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

    // // 保存 upc
    // compute_upc_and_save(upc);

    // 保存 upc 到 upc_out
    for(uint8_t slice_i = 0; slice_i < SLICES; slice_i++)
    {
      upc_out->val[i].slice[slice_i] = upc.slice[slice_i];
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
decode(OUT split_e_t       *e,
       IN const split_e_t  *R_e,
       IN const syndrome_t *original_s,
       IN const ct_t       *ct,
       IN const sk_t       *sk,
       IN const uint8_t     delat)
{
  // 初始化黑灰数组
  split_e_t  black_e = {0};
  split_e_t  gray_e  = {0};
  // split_e_t  fixed_e = {0};
  syndrome_t s;

  // // 构建出循环矩阵的索引 h_matrix 方便后面使用
  // uint16_t sk_wlist_all_0[R_BITS][DV] = {0};
  // uint16_t sk_wlist_all_1[R_BITS][DV] = {0};

  // // 填充对应的索引值
  // // 填充第一行
  // for(uint16_t i_DV = 0; i_DV < DV; i_DV++)
  // {
  //   sk_wlist_all_0[0][i_DV] = sk->wlist[0].val[i_DV];
  //   sk_wlist_all_1[0][i_DV] = sk->wlist[1].val[i_DV];
  // }

  // // 填充后 2 - 11779 行
  // for(uint16_t i_r = 1; i_r < R_BITS; i_r++)
  // {
  //   for(uint16_t i_DV = 0; i_DV < DV; i_DV++)
  //   {
  //     sk_wlist_all_0[i_r][i_DV] = (sk_wlist_all_0[i_r - 1][i_DV] + 1) % R_BITS;
  //     sk_wlist_all_1[i_r][i_DV] = (sk_wlist_all_1[i_r - 1][i_DV] + 1) % R_BITS;
  //   }
  // }

  // // 初始化 fixed_e 为 R_e
  // fixed_e.val[0] = R_e->val[0];
  // fixed_e.val[1] = R_e->val[1];

  // Reset (init) the error because it is xored in the find_err funcitons.
  // 初始化 e
  memset(e, 0, sizeof(*e));
  s = *original_s;
  dup(&s);

  // // ---- test ----
  // uint8_t test_res = compute_th_R(sk_wlist_all_0, sk_wlist_all_1, R_e, &s);
  // printf("test_res = %u\n", test_res);

  // -- test --
  // for(uint16_t i_test_s = 0; i_test_s < 555; i_test_s++)
  // {
  //   printf("第 %u 个 syndrome->qw 的值为: %lu\n", i_test_s, s.qw[i_test_s]);
  // }

  // 记录顺序 R_e0 R_e1 断行

  // 用于保存 s 的重量
  uint32_t s_array[MAX_IT + 2] = {0};
  // 用于保存 th
  uint32_t th_array[MAX_IT + 2] = {0};
  // 用于保存 upc 的值
  upc_all_t upc_array[MAX_IT + 2] = {0};
  // 记录 e
  split_e_t e_array[MAX_IT + 2] = {0};
  // 记录 black_e
  split_e_t black_e_array[MAX_IT + 2] = {0};
  // 记录 gray_e
  split_e_t gray_e_array[MAX_IT + 2] = {0};

  // 将 p_p x 加入文件尾部, 并断行
  double_t p_p_array[MAX_IT + 2] = {0};
  double_t x_array[MAX_IT + 2]   = {0};

  // // 记录R_e0
  // fprintf_LE((const uint64_t *)R_e->val[0].raw, R_BITS);
  // // 记录R_e1
  // fprintf_LE((const uint64_t *)R_e->val[1].raw, R_BITS);

  // // 断行
  // FILE *fp_2;
  // fp_2 = fopen("iter_data.txt", "a");
  // fprintf(fp_2, "\n");
  // fclose(fp_2);

  // 进入大迭代过程(for itr in 1...XBG do:)
  for(uint32_t iter = 0; iter < MAX_IT; iter++)
  {
    // // 将 fixed_e 和 求出来的 e 异或
    // GUARD(gf2x_add((uint8_t *)&fixed_e.val[0].raw, R_e->val[0].raw, e->val[0].raw,
    //                R_SIZE));
    // GUARD(gf2x_add((uint8_t *)&fixed_e.val[1].raw, R_e->val[1].raw, e->val[1].raw,
    //                R_SIZE));

    // 解码器使用阈值(th)来决定某个位是否为错误位
    // 该位确是错误位的概率随着间隙(upc[i] - th)的增加而增加
    // 该算法记录黑/灰掩码中有小间隙的位，以便后续步骤II和步骤III可以使用掩码，以获得翻转位的更多信息
    // 22: th = computeThreshold(s)
    // 参: Bit Flipping Key Encapsulation(v2.1) 17页，Threshold Selection Rule
    // printf("\n---->当前迭代阶段: %d<----\n", iter);

    // // 获取当前 fixed_e 的重量
    // uint16_t fixed_e_weight = r_bits_vector_weight(&fixed_e.val[0]) +
    //                           r_bits_vector_weight(&fixed_e.val[1]);

    // // ---- test ---- 使用论文方法计算的 th, 同时获取一个概率验证 p_p 和 x
    double_t       p_p         = 0;
    double_t       x           = 0;
    // const double_t threshold_2 = compute_th_R(
    //     &x, &p_p, sk_wlist_all_0, sk_wlist_all_1, fixed_e_weight, &fixed_e, &s);

    // 获取当前 s 的重量
    uint16_t s_weight = r_bits_vector_weight((const r_t *)s.qw);
    if(iter == 0)
    {
      s_array[iter] = s_weight;
    }
    else
    {
      s_array[iter + 2] = s_weight;
    }

    // 判断使用哪种方式计算 th
    uint8_t threshold = 0;
    // FILE   *fp_1;
    // fp_1 = fopen("iter_data.txt", "a");
    // fprintf(fp_1, "%u ", s_weight);
    if(TH_SELECT == 0)
    {
      // 原线性拟合方法计算 th
      threshold = get_threshold(&s);
      // fprintf(fp_1, "%u ", threshold);

      if(iter == 0)
      {
        th_array[iter] = threshold;
      }
      else
      {
        th_array[iter + 2] = threshold;
      }
    }
    else
    {
      // 论文方法计算 th, 把 th 向上取整
      // threshold = (uint8_t)threshold_2 + 1;
      // fprintf(fp_1, "%f ", threshold_2);

      if(iter == 0)
      {
        th_array[iter] = threshold;
      }
      else
      {
        th_array[iter + 2] = threshold;
      }
    }
    // fclose(fp_1);

    DMSG("    Iteration: %d\n", iter);
    DMSG("    Weight of e: %lu\n",
         r_bits_vector_weight(&e->val[0]) + r_bits_vector_weight(&e->val[1]));
    DMSG("    Weight of syndrome: %lu\n", r_bits_vector_weight((r_t *)s.qw));

    // 23:  (s, e, black, gray) = BitFlipIter(s, e, th, H) . Step I
    // H -- sk->wlist
    // 进入 procedure BitFlipIter(s, e, th, H)
    if(iter == 0)
    {
      find_err1(&upc_array[iter], e, &black_e, &gray_e, &s, sk->wlist, threshold,
                delat);
    }
    else
    {
      find_err1(&upc_array[iter + 2], e, &black_e, &gray_e, &s, sk->wlist,
                threshold, delat);
    }

    // 记录 e black_e gray_e
    if(iter == 0)
    {
      e_array[iter].val[0]       = e->val[0];
      e_array[iter].val[1]       = e->val[1];
      black_e_array[iter].val[0] = black_e.val[0];
      black_e_array[iter].val[1] = black_e.val[1];
      gray_e_array[iter].val[0]  = gray_e.val[0];
      gray_e_array[iter].val[1]  = gray_e.val[1];
      p_p_array[iter]            = p_p;
      x_array[iter]              = x;
    }
    else
    {
      e_array[iter + 2].val[0]       = e->val[0];
      e_array[iter + 2].val[1]       = e->val[1];
      black_e_array[iter + 2].val[0] = black_e.val[0];
      black_e_array[iter + 2].val[1] = black_e.val[1];
      gray_e_array[iter + 2].val[0]  = gray_e.val[0];
      gray_e_array[iter + 2].val[1]  = gray_e.val[1];
      p_p_array[iter + 2]            = p_p;
      x_array[iter + 2]              = x;
    }

    // // 记录当前的 e0
    // fprintf_LE((uint64_t *)e->val[0].raw, R_BITS);
    // // 记录当前的 e1
    // fprintf_LE((uint64_t *)e->val[1].raw, R_BITS);
    // // 记录当前的 black_e0
    // fprintf_LE((uint64_t *)black_e.val[0].raw, R_BITS);
    // // 记录当前的 black_e1
    // fprintf_LE((uint64_t *)black_e.val[1].raw, R_BITS);
    // // 记录当前的 gray_e0
    // fprintf_LE((uint64_t *)gray_e.val[0].raw, R_BITS);
    // // 记录当前的 gray_e1
    // fprintf_LE((uint64_t *)gray_e.val[1].raw, R_BITS);

    // // 将 p_p x 加入文件尾部, 并断行
    // FILE *fp_3;
    // fp_3 = fopen("iter_data.txt", "a");
    // fprintf(fp_3, "%f %f", p_p, x);
    // fprintf(fp_3, "\n");
    // fclose(fp_3);

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

    // 记录 s 和 th
    // 获取当前 s 的重量
    uint16_t s_weight_2 = r_bits_vector_weight((const r_t *)s.qw);
    // FILE    *fp_4;
    // fp_4 = fopen("iter_data.txt", "a");
    // fprintf(fp_4, "%u ", s_weight_2);
    // fprintf(fp_4, "%u ", ((DV + 1) / 2) + 1);
    // fclose(fp_4);

    // 24:  (s, e) = BitFlipMaskedIter(s, e, black, ((d + 1)/2), H) . Step II
    // procedure BitFlipMaskedIter(s, e, mask, th, H)
    find_err2(&upc_array[1], e, &black_e, &s, sk->wlist, ((DV + 1) / 2) + 1);
    GUARD(recompute_syndrome(&s, ct, sk, e));

    DMSG("    Weight of e: %lu\n",
         r_bits_vector_weight(&e->val[0]) + r_bits_vector_weight(&e->val[1]));
    DMSG("    Weight of syndrome: %lu\n", r_bits_vector_weight((r_t *)s.qw));

    // 记录 s 和 th
    s_array[1]  = s_weight_2;
    th_array[1] = ((DV + 1) / 2) + 1;
    // 记录 e black_e gray_e
    e_array[1].val[0]       = e->val[0];
    e_array[1].val[1]       = e->val[1];
    black_e_array[1].val[0] = black_e.val[0];
    black_e_array[1].val[1] = black_e.val[1];
    gray_e_array[1].val[0]  = gray_e.val[0];
    gray_e_array[1].val[1]  = gray_e.val[1];
    p_p_array[1]            = p_p;
    x_array[1]              = x;

    // // 记录当前的 e0
    // fprintf_LE((uint64_t *)e->val[0].raw, R_BITS);
    // // 记录当前的 e1
    // fprintf_LE((uint64_t *)e->val[1].raw, R_BITS);
    // // 记录当前的 black_e0
    // fprintf_LE((uint64_t *)black_e.val[0].raw, R_BITS);
    // // 记录当前的 black_e1
    // fprintf_LE((uint64_t *)black_e.val[1].raw, R_BITS);
    // // 记录当前的 gray_e0
    // fprintf_LE((uint64_t *)gray_e.val[0].raw, R_BITS);
    // // 记录当前的 gray_e1
    // fprintf_LE((uint64_t *)gray_e.val[1].raw, R_BITS);

    // // 将 p_p 加入尾部, 并断行
    // FILE *fp_5;
    // fp_5 = fopen("iter_data.txt", "a");
    // fprintf(fp_5, "%f %f", p_p, x);
    // fprintf(fp_5, "\n");
    // fclose(fp_5);

    // 记录 s 和 th
    // 获取当前 s 的重量
    uint16_t s_weight_3 = r_bits_vector_weight((const r_t *)s.qw);
    // FILE    *fp_6;
    // fp_6 = fopen("iter_data.txt", "a");
    // fprintf(fp_6, "%u ", s_weight_3);
    // fprintf(fp_6, "%u ", ((DV + 1) / 2) + 1);
    // fclose(fp_6);

    // 25:  (s, e) = BitFlipMaskedIter(s, e, gray, ((d + 1)/2), H) . Step III
    // procedure BitFlipMaskedIter(s, e, mask, th, H)
    find_err2(&upc_array[2], e, &gray_e, &s, sk->wlist, ((DV + 1) / 2) + 1);

    // 记录 s 和 th
    s_array[2]  = s_weight_3;
    th_array[2] = ((DV + 1) / 2) + 1;
    // 记录 e black_e gray_e
    e_array[2].val[0]       = e->val[0];
    e_array[2].val[1]       = e->val[1];
    black_e_array[2].val[0] = black_e.val[0];
    black_e_array[2].val[1] = black_e.val[1];
    gray_e_array[2].val[0]  = gray_e.val[0];
    gray_e_array[2].val[1]  = gray_e.val[1];
    p_p_array[2]            = p_p;
    x_array[2]              = x;

    // // 记录当前的 e0
    // fprintf_LE((uint64_t *)e->val[0].raw, R_BITS);
    // // 记录当前的 e1
    // fprintf_LE((uint64_t *)e->val[1].raw, R_BITS);
    // // 记录当前的 black_e0
    // fprintf_LE((uint64_t *)black_e.val[0].raw, R_BITS);
    // // 记录当前的 black_e1
    // fprintf_LE((uint64_t *)black_e.val[1].raw, R_BITS);
    // // 记录当前的 gray_e0
    // fprintf_LE((uint64_t *)gray_e.val[0].raw, R_BITS);
    // // 记录当前的 gray_e1
    // fprintf_LE((uint64_t *)gray_e.val[1].raw, R_BITS);

    // // 将 p_p,x 加入尾部, 并断行
    // FILE *fp_7;
    // fp_7 = fopen("iter_data.txt", "a");
    // fprintf(fp_7, "%f %f", p_p, x);
    // fprintf(fp_7, "\n");
    // fclose(fp_7);

    GUARD(recompute_syndrome(&s, ct, sk, e));
  }

  //  26: if (wt(s) != 0) then
  //  27:     return ⊥(ERROR)
  if(r_bits_vector_weight((r_t *)s.qw) > 0)
  {
    FILE *fp_3;
    fp_3 = fopen("weight_bad.txt", "a");
    fprintf(fp_3, "DELAT: %d 黑灰译码失败\n", delat);
    fclose(fp_3);
    // *flag = 1;

    // ---- 将记录的数据写入文件 ----
    // 首先写入 R_e
    fprintf_LE_test((const uint64_t *)R_e->val[0].raw, R_BITS);
    fprintf_LE_test((const uint64_t *)R_e->val[1].raw, R_BITS);
    // 断行
    FILE *fp_LE_test_1;
    fp_LE_test_1 = fopen("iter_data_all.txt", "a");
    fprintf(fp_LE_test_1, "\n");
    fclose(fp_LE_test_1);

    // 写入迭代数据
    for(uint8_t iter_i = 0; iter_i < MAX_IT + 2; iter_i++)
    {
      // 写入 s th
      FILE *fp_iter;
      fp_iter = fopen("iter_data_all.txt", "a");
      fprintf(fp_iter, "%u %u ", s_array[iter_i], th_array[iter_i]);
      fclose(fp_iter);
      // 写入 upc
      compute_upc_and_save_test(upc_array[iter_i].val[0]);
      compute_upc_and_save_test(upc_array[iter_i].val[1]);
      // 写入 e black_e gray_e p_p x
      // 记录当前的 e0
      fprintf_LE_test((uint64_t *)e_array[iter_i].val[0].raw, R_BITS);
      // 记录当前的 e1
      fprintf_LE_test((uint64_t *)e_array[iter_i].val[1].raw, R_BITS);
      // 记录当前的 black_e0
      fprintf_LE_test((uint64_t *)black_e_array[iter_i].val[0].raw, R_BITS);
      // 记录当前的 black_e1
      fprintf_LE_test((uint64_t *)black_e_array[iter_i].val[1].raw, R_BITS);
      // 记录当前的 gray_e0
      fprintf_LE_test((uint64_t *)gray_e_array[iter_i].val[0].raw, R_BITS);
      // 记录当前的 gray_e1
      fprintf_LE_test((uint64_t *)gray_e_array[iter_i].val[1].raw, R_BITS);
      // 将 p_p,x 加入尾部, 并断行
      FILE *fp_iter_2;
      fp_iter_2 = fopen("iter_data_all.txt", "a");
      fprintf(fp_iter_2, "%f %f", p_p_array[iter_i], x_array[iter_i]);
      fprintf(fp_iter_2, "\n");
      fclose(fp_iter_2);
    }

    DMSG("s 重量不为 0...");
    BIKE_ERROR(E_DECODING_FAILURE);
  }

  return SUCCESS;
}
