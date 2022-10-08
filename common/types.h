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
 * Written by Nir Drucker and Shay Gueron
 * AWS Cryptographic Algorithms Group.
 * (ndrucker@amazon.com, gueron@amazon.com)
 */

#pragma once

#include "bike_defs.h"
#include "error.h"
#include <stdint.h>

typedef struct uint128_s
{
  union
  {
    uint8_t  bytes[16];
    uint32_t dw[4];
    uint64_t qw[2];
  } u;
} uint128_t;

// Make sure no compiler optimizations.
#pragma pack(push, 1)

typedef struct seed_s
{
  uint8_t raw[32];
} seed_t;

typedef struct seeds_s
{
  seed_t seed[NUM_OF_SEEDS];
} seeds_t;

// raw[R_SIZE] = raw[1473] 1473
typedef struct r_s
{
  uint8_t raw[R_SIZE];
} r_t;

typedef struct e_s
{
  uint8_t raw[N_SIZE];
} e_t;

typedef struct generic_param_n_s
{
  r_t val[N0];
} generic_param_n_t;

typedef generic_param_n_t ct_t;
typedef generic_param_n_t pk_t;
typedef generic_param_n_t split_e_t;

typedef uint32_t idx_t;

typedef struct compressed_idx_dv_s
{
  idx_t val[DV];
} compressed_idx_dv_t;

// val[71] compressed_idx_dv_ar_t[2]
typedef compressed_idx_dv_t compressed_idx_dv_ar_t[N0];

typedef struct compressed_idx_t_t
{
  idx_t val[T1];
} compressed_idx_t_t;

// The secret key holds both representation for avoiding
// the compression in the decaps stage
//
// r_t                    bin[2];
// compressed_idx_dv_ar_t wlist;  wlist[ val[71], val[71] ]
// r_t                    sigma0;
// r_t                    sigma1;
typedef struct sk_s
{
  r_t                    bin[N0]; // bin 每个元素为 1473 int8_t 数组
  compressed_idx_dv_ar_t wlist; // 密钥汉明重位置 wlist[ val[71], val[71] ]
  r_t                    sigma0;
  r_t                    sigma1;
} sk_t;

// Pad e to the next Block
typedef ALIGN(8) struct padded_e_s
{
  e_t     val;
  uint8_t pad[N_PADDED_SIZE - N_SIZE];
} padded_e_t;

// Pad r to the next Block
typedef ALIGN(8) struct padded_r_s
{
  r_t     val;                         // raw[1473]
  uint8_t pad[R_PADDED_SIZE - R_SIZE]; // pad[575]
} padded_r_t;

typedef padded_r_t       padded_param_n_t[N0];
typedef padded_param_n_t pad_sk_t;
typedef padded_param_n_t pad_pk_t;
typedef padded_param_n_t pad_ct_t;

// Need to allocate twice the room for the results
// 需要为结果分配两倍的空间
typedef ALIGN(8) struct dbl_padded_r_s
{
  r_t     val;
  uint8_t pad[(2 * R_PADDED_SIZE) - R_SIZE]; // pad[2623]
} dbl_padded_r_t;

typedef dbl_padded_r_t       dbl_padded_param_n_t[N0];
typedef dbl_padded_param_n_t dbl_pad_pk_t;
typedef dbl_padded_param_n_t dbl_pad_ct_t;
typedef dbl_padded_param_n_t dbl_pad_syndrome_t;

typedef struct ss_s
{
  uint8_t raw[ELL_K_SIZE];
} ss_t;

// For optimization purposes
//  1- For a faster rotate we duplicate the syndrome (dup1/2)
//  2- We extend it to fit the boundary of DDQW
// 出于优化目的
// 1- 为了更快的旋转，我们复制了 s (dup1/2)
// 2- 我们扩展它以适应 DDQW 的边界
typedef ALIGN(64) struct syndrome_s
{
  uint64_t qw[3 * R_QW];
} syndrome_t;

// 利用 syndrome_t 的结构封装 ct
typedef struct dup_c_s
{
  syndrome_t val[N0];
} dup_c_t;

// 为了对每行 H0 和 H1 旋转，我们类似复制一次 h
typedef ALIGN(64) struct single_h_s
{
  uint64_t qw[2 * R_QW];
} single_h_t;

// 利用 h_t 的结构封装 h
typedef struct h_s
{
  single_h_t val[N0];
} h_t;

// modify: equations_t -> uint16_t equations
// // 存放方程索引的二维数组
// // 这里用 10 列基本可以包含所有未知数索引
// typedef struct single_equation_s
// {
//   uint16_t eq[R_BITS][10];
// } single_equation_t;

// // 封装两个二维数组
// typedef struct equation_s
// {
//   single_equation_t val[N0];
// } equations_t;

typedef struct upc_slice_s
{
  union
  {
    padded_r_t r;
    uint64_t   qw[sizeof(padded_r_t) / 8];
  } u;
} upc_slice_t;

typedef struct upc_s
{
  upc_slice_t slice[SLICES]; // SLICES=8
} upc_t;

#pragma pack(pop)
