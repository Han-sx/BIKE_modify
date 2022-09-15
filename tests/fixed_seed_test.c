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

#include "kem.h"
#include "measurements.h"
#include "utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

////////////////////////////////////////////////////////////////
//                 Main function for testing
////////////////////////////////////////////////////////////////
int
main()
{
#ifdef FIXED_SEED
  srand(0);
#else
  srand(time(NULL));
#endif

  uint8_t sk[sizeof(sk_t)]    = {0}; // private-key: (h0, h1)
  uint8_t pk[sizeof(pk_t)]    = {0}; // public-key:  (g0, g1)
  uint8_t ct[sizeof(ct_t)]    = {0}; // ciphertext:  (c0, c1)
  uint8_t k_enc[sizeof(ss_t)] = {0}; // shared secret after encapsulate
  uint8_t k_dec[sizeof(ss_t)] = {0}; // shared secret after decapsulate

  for(uint32_t i = 1; i <= NUM_OF_TESTS; ++i)
  {
    int res = 0;

    MSG("Code test: %d\n\n", i);

    // Key generation 密钥生成
    // h0 h1, σ0 σ1, g 在此步骤中均随机采样
    MEASURE("  keypair", res = crypto_kem_keypair(pk, sk););

    if(res != 0)
    {
      MSG("Keypair failed with error: %d\n", res);
      continue;
    }

    uint32_t dec_rc = 0;

    // Encapsulate 密钥封装，IN pk, OUT ct and k_enc
    // m 在此步骤中随机采样，过程中获取到 (e0,e1), (c0,c1)，协商密钥 k
    MEASURE("  encaps", res = crypto_kem_enc(ct, k_enc, pk););
    if(res != 0)
    {
      MSG("encapsulate failed with error: %d\n", res);
      continue;
    }

    // Decapsulate 解封装, IN ct and sk, OUT k_dec
    MEASURE("  decaps", dec_rc = crypto_kem_dec(k_dec, ct, sk););

    if(dec_rc != 0)
    {
      printf("Decoding failed after %d code tests!\n", i);
    }
    else
    {
      if(secure_cmp(k_enc, k_dec, sizeof(k_dec) / sizeof(uint64_t)))
      {
        MSG("Success! decapsulated key is the same as encapsulated "
            "key!\n");
      }
      else
      {
        MSG("Failure! decapsulated key is NOT the same as encapsulated "
            "key!\n");
      }
    }

    print("Initiator's generated key (K) of 256 bits = ", (uint64_t *)k_enc,
          SIZEOF_BITS(k_enc));
    print("Responder's computed key (K) of 256 bits  = ", (uint64_t *)k_dec,
          SIZEOF_BITS(k_enc));
  }

  return 0;
}
