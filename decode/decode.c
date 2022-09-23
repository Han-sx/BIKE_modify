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
#include "utilities.h"
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
#    define MAX_IT 6
#  elif(LEVEL == 5)
#    define MAX_IT 7
#  else
#    error "Level can only be 1/3/5"
#  endif
#endif

#define EQ_COLUMN 11 // 索引矩阵列数

// 定义一个全局变量用于记录 equations 中每一行的非零索引
uint8_t eq_index[R_BITS] = {0};

// // 用于 将一个 8 位数组转换为十进制 存放于 bin_to_uint8_tmp 中
// void binary_to_uint8_t(uint8_t a[])
// {
// 	for(uint8_t i = 0; i < 8; i++)
// 	{
// 		bin_to_uint8_tmp += a[i]*power(2, 7 - i);
// 	}
// }

// // 用于 uint8_t 转换为二进制 保存在 8 位数组中
// void
// uint8_t_toBinary(uint8_t b[], uint8_t a)
// {
//   uint8_t i;
//   uint8_t count = 0;

//   for(i = 0x80; i != 0; i >>= 1)
//   {
//     b[count] = (a & i) ? 1 : 0;
//     count++;
//   }
// }

// // 数组拼接，将 8 位长数组 b 的值连接到数组 a 后
// void
// array_concatenation(uint32_t index, uint8_t a[], uint8_t b[])
// {
//   for(uint32_t i = 0; i < 8; i++)
//   {
//     a[index] = b[i];
//     index++;
//   }
// }

// TODO 将增广常数的值赋值给索引矩阵
// _INLINE_ void
// assign_constants(OUT uint16_t *out, IN const uint8_t *in)
// {
//   for(size_t i = 0; i < R_QW - 1; i++){
//     for(uint8_t index = 0, location = 1; location != 0; location <<= 1)
//     {
//       if(location & in[i] != 0){

//       }
//       index++;
//     }
//   }
// }

// 对 qw[2 * R_QW] 循环右移一位, 仅 [185]-[369] 有效
_INLINE_ void
rotate_right_one(OUT single_h_t *out, IN const single_h_t *in)
{
  for(size_t i = 369; i > R_QW - 1; i--)
  {
    out->qw[i] = (in->qw[i] << 1) | (in->qw[i - 1] >> 63);
  }
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
    res[i] = ~a[i] & b[i];
  }
  return SUCCESS;
}

// 对 bytelen 长字节流, a 和 b '与', 索引存储在 res 行
// 注意：此处的索引和 wlist 不同的是从 1 开始
_INLINE_ ret_t
and_index(OUT uint16_t     *res,
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
          eq_index[i_eq]++;
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
      eq_index[i_eq]++;
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
  h->qw[0] = (h->qw[185] << LAST_R_QW_TRAIL) | (h->qw[369] << LAST_R_QW_TRAIL_2) |
             (h->qw[368] >> LAST_R_QW_LEAD_2);
  // qw[0] = (h2,h1,h0,h11778,...,h11718)

  for(size_t i = 184; i > 0; i--)
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

  // 打印 pad_s 的值
  // for(uint16_t i_pad_s = 0; i_pad_s < 1473; i_pad_s++)
  // {
  //   printf("第 %u 个 s 的值为: %u\n", i_pad_s, pad_s[0].val.raw[i_pad_s]);
  // }

  // 将 s 按照每 8 位存储在 1473 个字节中
  // 复制 1473 个字节到 qw 的前 185 个 64 位整型中
  memcpy((uint8_t *)syndrome->qw, pad_s[0].val.raw, R_SIZE);

  // 打印 syndrome->qw 中的值
  // for(uint16_t i_qw = 0; i_qw < 555; i_qw++)
  // {
  //   printf("第 %u 个 syndrome->qw 的值为: %lu\n", i_qw, syndrome->qw[i_qw]);
  // }

  dup(syndrome);
  // 打印复制后 syndrome->qw 中的值
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
  return threshold;
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
          IN const uint8_t                threshold)
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
    for(size_t l = 0; l < DELTA; l++)
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
decode(OUT split_e_t       *e,
       IN const syndrome_t *original_s,
       IN const ct_t       *ct,
       IN const sk_t       *sk)
{
  // 初始化黑灰数组
  split_e_t  black_e         = {0};
  split_e_t  gray_e          = {0};
  split_e_t  black_or_gray_e = {0};
  ct_t       ct_remove_BG    = {0};
  ct_t       ct_pad          = {0};
  dup_c_t    c               = {0};
  dup_c_t    rotated_c       = {0};
  dup_c_t    constant_term   = {0};
  h_t        h               = {0};
  syndrome_t s;
  // 定义 11779 行方程组, 前 10 个元素用于保存索引, 第 11 个用于存放增广常数
  uint16_t equations[R_BITS][EQ_COLUMN] = {0};
  // // 定义 求解矩阵
  // uint8_t matrix_solver[R_BITS][R_BITS + 1] = {0};
  // // 定义 结果向量
  // uint8_t x_solution[2 * R_BITS] = {0};

  // 获取 ct 的值
  ct_pad.val[0] = ct->val[0];
  ct_pad.val[1] = ct->val[1];

  // 从 sk 中获取 h 第一行的 bin
  // 复制 1473 个字节到 qw 的前 185 个 64 位整型中
  memcpy((uint8_t *)&h.val[0].qw[185], sk->bin[0].raw, R_SIZE);
  memcpy((uint8_t *)&h.val[1].qw[185], sk->bin[1].raw, R_SIZE);

  // 复制 h
  dup_two(&h.val[0]);
  dup_two(&h.val[1]);

  // Reset (init) the error because it is xored in the find_err funcitons.
  // 初始化 e
  memset(e, 0, sizeof(*e));
  s = *original_s;
  dup(&s);

  // -- test --
  // for(uint16_t i_test_s = 0; i_test_s < 555; i_test_s++)
  // {
  //   printf("第 %u 个 syndrome->qw 的值为: %lu\n", i_test_s, s.qw[i_test_s]);
  // }

  // 进入大迭代过程(for itr in 1...XBG do:)
  for(uint32_t iter = 0; iter < MAX_IT; iter++)
  {
    // 解码器使用阈值(th)来决定某个位是否为错误位
    // 该位确是错误位的概率随着间隙(upc[i] - th)的增加而增加
    // 该算法记录黑/灰掩码中有小间隙的位，以便后续步骤II和步骤III可以使用掩码，以获得翻转位的更多信息
    // 22: th = computeThreshold(s)
    // 参: Bit Flipping Key Encapsulation(v2.1) 17页，Threshold Selection Rule
    printf("\n---->当前迭代阶段: %d<----\n", iter);

    // // 优先获取此轮译码开始的初始 c 值, 获取二进制长度的 c_bin
    // // 将 1473 长度的十进制 c0 转换为 11784 长度的 c_bin
    // for(uint32_t i_c0 = 0; i_c0 < 1473; i_c0++)
    // {
    //   uint8_t_toBinary(c_pad_tmp, ct->val[0].raw[i_c0]);
    //   array_concatenation(i_c0 * 8, c_bin, c_pad_tmp);
    // }
    // // 将 1473 长度的十进制 h1 转换为 23568 长度的 c_bin
    // for(uint32_t i_c1 = 0; i_c1 < 1473; i_c1++)
    // {
    //   uint8_t_toBinary(c_pad_tmp, ct->val[1].raw[i_c1]);
    //   array_concatenation(1473 + i_c1 * 8, c_bin, c_pad_tmp);
    // }

    const uint8_t threshold = get_threshold(&s);

    DMSG("    Iteration: %d\n", iter);
    DMSG("    Weight of e: %lu\n",
         r_bits_vector_weight(&e->val[0]) + r_bits_vector_weight(&e->val[1]));
    DMSG("    Weight of syndrome: %lu\n", r_bits_vector_weight((r_t *)s.qw));

    // 23:  (s, e, black, gray) = BitFlipIter(s, e, th, H) . Step I
    // H -- sk->wlist
    // 进入 procedure BitFlipIter(s, e, th, H)
    find_err1(e, &black_e, &gray_e, &s, sk->wlist, threshold);

    // // 输出black_e
    // printf("\n第 %d 轮迭代的black_e:\n", iter);
    // print("\nblack_e0: \n", (uint64_t *)black_e.val[0].raw, R_BITS);
    // print("\nblack_e1: \n", (uint64_t *)black_e.val[1].raw, R_BITS);

    // 输出 black_e 和 gray_e 的重量
    printf("\nblack_e 的重量：%lu \n",
           (r_bits_vector_weight((r_t *)black_e.val[0].raw) +
            r_bits_vector_weight((r_t *)black_e.val[1].raw)));
    printf("\ngray_e 的重量：%lu \n",
           (r_bits_vector_weight((r_t *)gray_e.val[0].raw) +
            r_bits_vector_weight((r_t *)gray_e.val[1].raw)));

    // 输出当前迭代的第 I 步骤中的 e
    printf("\n第 %d 轮迭代的 e:\n", iter);
    print("\ntmp_find_e0: \n", (uint64_t *)e->val[0].raw, R_BITS);
    print("\ntmp_find_e1: \n", (uint64_t *)e->val[1].raw, R_BITS);

    // 10:  s = H(cT + eT ) . 更新校验子 syndrome
    GUARD(recompute_syndrome(&s, ct, sk, e));

// 此处代码中在 iter >= 1 时候去除了 Step II 和 Step III
// 相当于只进行了一轮黑灰迭代后进行了多轮比特位反转(step I)
// 在论文中没有对 iter 迭代次数的判断依据
#ifdef BGF_DECODER
    if(iter >= 1)
    {
      continue;
    }
#endif
    DMSG("    Weight of e: %lu\n",
         r_bits_vector_weight(&e->val[0]) + r_bits_vector_weight(&e->val[1]));
    DMSG("    Weight of syndrome: %lu\n", r_bits_vector_weight((r_t *)s.qw));

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

    // ---------> 增加方程组求解算法(当 s 不为 0) <---------
    // =================================================================

    for(uint32_t i = 0; i < N0; i++)
    {
      // 将黑灰集合'或'运算(black_e | gray_e) 存放于
      // black_or_gray_e，即所有未知数位
      GUARD(gf2x_add(black_or_gray_e.val[i].raw, black_e.val[i].raw,
                     gray_e.val[i].raw, R_SIZE));

      // 去除 c 中的未知数位，将 black_or_gray_e 取反后与 c 做与操作
      GUARD(negate_and(ct_remove_BG.val[i].raw, black_or_gray_e.val[i].raw,
                       ct_pad.val[i].raw, R_SIZE));

      // 将 ct_remove_BG 的 uint8 存储结构调整位 uint64
      // 调整 1473 个字节到 qw 的前 185 个 64 位整型中，并复制三份
      memcpy((uint8_t *)c.val[i].qw, ct_remove_BG.val[i].raw, R_SIZE);
      dup(&c.val[i]);

      // 对每个密钥集位索引的校正子进行右循环，这里表示 ct_remove_BG 乘 H 转置
      for(size_t j = 0; j < DV; j++)
      {
        // 输出校验子仅包含一个 R_BITS 旋转，其他 (2 * R_BITS) 位未定义
        rotate_right(&rotated_c.val[i], &c.val[i], sk->wlist[i].val[j]);

        // 将每个 rotated_c.val[i] 进行异或相加, 保存到 constant_term 中
        GUARD(gf2x_add((uint8_t *)&constant_term.val[i].qw,
                       (uint8_t *)constant_term.val[i].qw,
                       (uint8_t *)rotated_c.val[i].qw, R_SIZE));
      }

      // 对方程组未知数进行构建，两次循环的索引(从 1 开始)都存储于 equeations 中
      for(uint16_t i_eq = 0; i_eq < R_BITS; i_eq++)
      {
        // 将当前 h 与 black_or_gray_e 与运算
        // h 的有效位是 [185]-[369]
        GUARD(and_index(equations[i_eq], black_or_gray_e.val[i].raw,
                        (uint8_t *)&h.val[i].qw[R_QW], R_SIZE, i, i_eq));
        // 对 H 进行 1 bit 循环右移位
        rotate_right_one(&h.val[i], &h.val[i]);
      }
    }

    // 将 constant_term.val[0] ^ constant_term.val[1] 放入 constant_term.val[0]
    GUARD(gf2x_add((uint8_t *)&constant_term.val[0].qw,
                   (uint8_t *)constant_term.val[0].qw,
                   (uint8_t *)constant_term.val[1].qw, R_SIZE));

    // --------------------- 整合解方程函数 ---------------------
    // 将增广常数 constant_term.val[0].qw 赋值给 equations[i][EQ_COLUMN]
    // 处理前 11776 位
    for(uint8_t i = 0; i < R_QW - 1; i++)
    {
      for(uint64_t index = 0, location = 1; location != 0; location <<= 1)
      {
        if((constant_term.val[0].qw[i] & location) != 0)
        {
          equations[64 * i + index][EQ_COLUMN - 1] = 1;
        }
        index++;
      }
    }
    // 处理最后三位
    for(uint64_t index = 0, location = 1; location < 8; location <<= 1)
    {
      if((constant_term.val[0].qw[R_QW] & location) != 0)
      {
        equations[64 * (R_QW - 1) + index][EQ_COLUMN - 1] = 1;
      }
      index++;
    }

    // 给 matrix_solver 矩阵赋值
    // // TODO 第 11 个元素未赋值
    // for(uint16_t i = 0; i < R_BITS; i++)
    // {
    //   for(uint8_t j = 0; j < (EQ_COLUMN - 1); j++)
    //   {
    //     if(equations[i][j] != 0)
    //     {
    //       matrix_solver[i][equations[i][j] - 1] = 1;
    //     }
    //   }
    // }

    // // 查看需要求解的未知数个数
    // printf("\nblack_or_gray_e 的未知数个数：%lu \n",
    //        (r_bits_vector_weight((r_t *)black_or_gray_e.val[0].raw) +
    //         r_bits_vector_weight((r_t *)black_or_gray_e.val[1].raw)));

    // -- test -- 输出 equations 的值, 并保存到 data_1.txt 中
    FILE *fp;
    fp = fopen("data_1.txt", "a");
    for(uint16_t i = 0; i < 11779; i++)
    {
      // if(equations[i][0] == 0)
      // {
      //   continue;
      // }
      for(uint8_t j = 0; j < 11; j++)
      {
        if(j == 10){
          fprintf(fp, "%u\n", equations[i][j]);
          continue;
        }
        fprintf(fp, "%u,", equations[i][j]);
        // if(equations[i][j] != 0)
        // {
        //   printf("%u  ", equations[i][j]);
        // }
      }
      // printf("\n");
    }
    fclose(fp);

    // =================================================================
  }

  // 打印当前 e 查看译码结果
  DMSG("\n---->当前译码获得的错误向量如下:<----\n\n")
  print("\ndecode_e0: \n", (uint64_t *)e->val[0].raw, R_BITS);
  print("\ndecode_e1: \n", (uint64_t *)e->val[1].raw, R_BITS);

  //  26: if (wt(s) != 0) then
  //  27:     return ⊥(ERROR)
  if(r_bits_vector_weight((r_t *)s.qw) > 0)
  {
    DMSG("s 重量不为 0...");
    BIKE_ERROR(E_DECODING_FAILURE);
  }

  return SUCCESS;
}
