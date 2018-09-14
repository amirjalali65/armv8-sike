/********************************************************************************************
* Abstract: Finite field arithmetic for ARM64 using code modified from the original x86_64
*           and generic implementations by microsoft.
*
* The multiplication is implemented using A64 + ASIMD mixed implementation 
* and 964-bit 2-level Karatsuba multiplication
* Written by Amir Jalali from the base generic implementation by Microsoft Research 
*
* ajalali2016@fau.edu    
* Florida Atlantic University
* All rights reserved        
*
*********************************************************************************************/

#include "../P964_internal.h"


// Global constants
extern const uint64_t p964[NWORDS_FIELD];
extern const uint64_t p964p1[NWORDS_FIELD]; 

#include <stdio.h>
#include <stdlib.h>

__inline void fpadd964(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular addition, c = a+b mod p964.
  // Inputs: a, b in [0, p964-1] 
  // Output: c in [0, p964-1]

    fpadd964_arm64_asm(a, b, c);
} 


__inline void fpsub964(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular subtraction, c = a-b mod p964.
  // Inputs: a, b in [0, p964-1] 
  // Output: c in [0, p964-1] 

    fpsub964_arm64_asm(a, b, c);
}


//  identical to x86_64 implementation
__inline void fpneg964(digit_t* a)
{ // Modular negation, a = -a mod p964.
  // Input/output: a in [0, p964-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, ((digit_t*)p964)[i], a[i], borrow, a[i]); 
    }
}


//  identical to x86_64 implementation
void fpdiv2_964(const digit_t* a, digit_t* c)
{ // Modular division by two, c = a/2 mod p964.
  // Input : a in [0, p964-1] 
  // Output: c in [0, p964-1] 
    unsigned int i, carry = 0;
    digit_t mask;
        
    mask = 0 - (digit_t)(a[0] & 1);    // If a is odd compute a+p521
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(carry, a[i], ((digit_t*)p964)[i] & mask, carry, c[i]); 
    }

    mp_shiftr1(c, NWORDS_FIELD);
} 

void fpcorrection964(digit_t* a)
{ // Modular correction to reduce field element a in [0, 2*p964-1] to [0, p964-1].
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], ((digit_t*)p964)[i], borrow, a[i]);
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, a[i], ((digit_t*)p964)[i] & mask, borrow, a[i]);
    }
}

void fpmul2x256_mixed_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{
    // Multiplication of a * b using mixed assembly
    // c[0...7] = a[0...3] * b[0...3] -> ASIMD assembly
    // c[8...15] = a[4...7] * b[4...7] -> A64 assembly
    asm(
        // load a[0...3] and b[0...3]
        "ldp x3, x4, [%0, #32]      \n\t"
        "ldp x5, x6, [%0, #48]      \n\t"
        "ldp x7, x8, [%1, #32]      \n\t"
        "ldp x9, x10, [%1, #48]     \n\t"

        // load a[4...7] and b[4...7]
        "ld1 {v0.2d, v1.2d}, [%0]  \n\t"
        "ld1 {v2.2d, v3.2d}, [%1]  \n\t"

        // Shuffle a[4...7] -> a[0|1|2|3|4|5|6|7] -> a[0|4|2|6|1|5|3|7]
        "trn1 v4.4s, v0.4s, v1.4s       \n\t"
        "trn2 v5.4s, v0.4s, v1.4s       \n\t"
        
        // Set x20 as zero register for A64 arithmetic
        "eor x20, x20, x20          \n\t"
        
        "mul x11, x3, x7            \n\t"
        "umulh x12, x3, x7          \n\t"

        // step 1
        "umull v6.2d, v4.2s, v2.s[0]    \n\t"   // v6 = a0b0 | a4b0   
        "umull v7.2d, v5.2s, v2.s[0]    \n\t"   // v7 = a1b0 | a5b0
        
        "mul x13, x3, x8            \n\t"
        "umulh x14, x3, x8          \n\t"

        "umull2 v8.2d, v4.4s, v2.s[0]   \n\t"   // v8 = a2b0 | a6b0
        "umull2 v9.2d, v5.4s, v2.s[0]   \n\t"   // v9 = a3b0 | a7b0

        "adds x25, x12, x13         \n\t"
        "adcs x26, x14, x20         \n\t"
        
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v11.16b, v11.16b, v11.16b  \n\t"
        "eor v28.16b, v28.16b, v28.16b  \n\t"   // for sepcial purpose
        "eor v27.16b, v27.16b, v27.16b  \n\t"   // for special purpose     

        "trn2 v0.4s, v6.4s, v0.4s       \n\t"   // v0 = h(a0b0 | a4b0)
        "trn1 v6.4s, v6.4s, v28.4s      \n\t"   // v6 = l(a0b0 | a4b0)
        "trn2 v1.4s, v7.4s, v1.4s       \n\t"   // v1 = h(a1b0 | a5b0)
        "trn1 v7.4s, v7.4s, v28.4s      \n\t"   // v7 = l(a1b0 | a5b0)
        "trn2 v10.4s, v8.4s, v10.4s     \n\t"   // v10 = h(a2b0 | a6b0)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a2b0 | a6b0)
        "trn2 v11.4s, v9.4s, v11.4s     \n\t"   // v11 = h(a3b0 | a7b0)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a3b0 | a7b0)

        "add v7.2d, v7.2d, v0.2d    \n\t"
        "add v8.2d, v8.2d, v1.2d    \n\t"
        "add v9.2d, v9.2d, v10.2d   \n\t"
        "mov v27.d[0], v6.d[1]      \n\t"
        "uqadd v11.2d, v11.2d, v27.2d \n\t"        

        //final results are in v20, v21, v22, v23
        "mov v20.s[0], v6.s[0]      \n\t"       // c0

        "mul x13, x4, x7            \n\t"
        "umulh x14, x4, x7          \n\t"

        //step 2
        "umlal v7.2d, v4.2s, v2.s[1]    \n\t"   // v7 += a0b1 | a4b1
        "umlal v8.2d, v5.2s, v2.s[1]    \n\t"   // v8 += a1b1 | a5b1


        "adds x25, x25, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"
        
        "umlal2 v9.2d, v4.4s, v2.s[1]   \n\t"   // v9 += a2b1 | a6b1
        "umlal2 v11.2d, v5.4s, v2.s[1]  \n\t"   // v11 += a3b1 | a7b1

        "stp x11, x25, [%2, #64]         \n\t"
        
        "mul x13, x3, x9            \n\t"
        "umulh x14, x3, x9          \n\t"
        
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v12.16b, v12.16b, v12.16b  \n\t"

        "trn2 v0.4s, v7.4s, v0.4s       \n\t"   // v0 = h(a0b1 | a4b1)
        "trn1 v7.4s, v7.4s, v28.4s      \n\t"   // v7 = l(a0b1 | a4b1)
        "trn2 v1.4s, v8.4s, v1.4s       \n\t"   // v1 = h(a1b1 | a5b1)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a1b1 | a5b1)
        "trn2 v10.4s, v9.4s, v10.4s     \n\t"   // v10 = h(a2b1 | a6b1)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a2b1 | a6b1)
        "trn2 v12.4s, v11.4s, v12.4s    \n\t"   // v12 = h(a3b1 | a7b1)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a3b1 | a7b1)

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "add v8.2d, v8.2d, v0.2d        \n\t"
        "add v9.2d, v9.2d, v1.2d        \n\t"
        "add v11.2d, v11.2d, v10.2d     \n\t"
        "mov v27.d[0], v7.d[1]          \n\t"
        "uqadd v12.2d, v12.2d, v27.2d   \n\t"        

        "mov v20.s[1], v7.s[0]      \n\t"       // c1

        "mul x13, x4, x8            \n\t"
        "umulh x14, x4, x8          \n\t"

        //step 3
        "umlal v8.2d, v4.2s, v2.s[2]    \n\t"   // v8 += a0b2 | a4b2
        "umlal v9.2d, v5.2s, v2.s[2]    \n\t"   // v9 += a1b2 | a5b2

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "umlal2 v11.2d, v4.4s, v2.s[2]   \n\t"   // v11 += a2b2 | a6b2
        "umlal2 v12.2d, v5.4s, v2.s[2]  \n\t"   // v12 += a3b2 | a7b2
       
        "mul x13, x5, x7            \n\t"
        "umulh x14, x5, x7          \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v6.16b, v6.16b, v6.16b     \n\t"

        "trn2 v0.4s, v8.4s, v0.4s       \n\t"   // v0 = h(a0b2 | a4b2)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a0b2 | a4b2)
        "trn2 v1.4s, v9.4s, v1.4s       \n\t"   // v1 = h(a1b2 | a5b2)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a1b2 | a5b2)
        "trn2 v10.4s, v11.4s, v10.4s    \n\t"   // v10 = h(a2b2 | a6b2)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a2b2 | a6b2)
        "trn2 v6.4s, v12.4s, v6.4s      \n\t"   // v6 = h(a3b2 | a7b2)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a3b2 | a7b2)

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "add v9.2d, v9.2d, v0.2d        \n\t"
        "add v11.2d, v11.2d, v1.2d      \n\t"
        "add v12.2d, v12.2d, v10.2d     \n\t"
        "mov v27.d[0], v8.d[1]          \n\t"
        "uqadd v6.2d, v6.2d, v27.2d     \n\t"        

        "mov v20.s[2], v8.s[0]      \n\t"       // c2
        
        "mul x13, x3, x10           \n\t"
        "umulh x14, x3, x10         \n\t"

        //step 4 
        "umlal v9.2d, v4.2s, v2.s[3]    \n\t"   // v9 += a0b3 | a4b3
        "umlal v11.2d, v5.2s, v2.s[3]   \n\t"   // v11 += a1b3 | a5b3

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x20, x20         \n\t"

        "umlal2 v12.2d, v4.4s, v2.s[3]  \n\t"   // v12 += a2b3 | a6b3
        "umlal2 v6.2d, v5.4s, v2.s[3]   \n\t"   // v6 += a3b3 | a7b3        

        "mul x13, x4, x9            \n\t"
        "umulh x14, x4, x9          \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v8.16b, v8.16b, v8.16b     \n\t"

        "trn2 v0.4s, v9.4s, v0.4s       \n\t"   // v0 = h(a0b3 | a4b3)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a0b3 | a4b3)
        "trn2 v1.4s, v11.4s, v1.4s      \n\t"   // v1 = h(a1b3 | a5b3)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a1b3 | a5b3)
        "trn2 v10.4s, v12.4s, v10.4s    \n\t"   // v10 = h(a2b3 | a6b3)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a2b3 | a6b3)
        "trn2 v8.4s, v6.4s, v8.4s       \n\t"   // v8 = h(a3b3 | a7b3)
        "trn1 v6.4s, v6.4s, v28.4s      \n\t"   // v6 = l(a3b3 | a7b3)
 
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"
        
        "add v11.2d, v11.2d, v0.2d      \n\t"
        "add v12.2d, v12.2d, v1.2d      \n\t"
        "add v6.2d, v6.2d, v10.2d       \n\t"
        "mov v27.d[0], v9.d[1]          \n\t"
        "uqadd v8.2d, v8.2d, v27.2d     \n\t"        

        "mov v20.s[3], v9.s[0]          \n\t"   // c3
                
        "mul x13, x5, x8            \n\t"
        "umulh x14, x5, x8          \n\t"

        //step 5
        "umlal v11.2d, v4.2s, v3.s[0]    \n\t"   // v11 += a0b4 | a4b4
        "umlal v12.2d, v5.2s, v3.s[0]    \n\t"   // v12 += a1b4 | a5b4

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "umlal2 v6.2d, v4.4s, v3.s[0]   \n\t"   // v6 += a2b4 | a6b4
        "umlal2 v8.2d, v5.4s, v3.s[0]   \n\t"   // v8 += a3b4 | a7b4        

        "mul x13, x6, x7            \n\t"
        "umulh x14, x6, x7          \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v7.16b, v7.16b, v7.16b     \n\t"

        "trn2 v0.4s, v11.4s, v0.4s      \n\t"   // v0 = h(a0b4 | a4b4)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a0b4 | a4b4)
        "trn2 v1.4s, v12.4s, v1.4s      \n\t"   // v1 = h(a1b4 | a5b4)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a1b4 | a5b4)
        "trn2 v10.4s, v6.4s, v10.4s     \n\t"   // v10 = h(a2b4 | a6b4)
        "trn1 v6.4s, v6.4s, v28.4s      \n\t"   // v6 = l(a2b4 | a6b4)
        "trn2 v7.4s, v8.4s, v7.4s       \n\t"   // v7 = h(a3b4 | a7b4)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a3b4 | a7b4)

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"
        
        "add v12.2d, v12.2d, v0.2d      \n\t"
        "add v6.2d, v6.2d, v1.2d        \n\t"
        "add v8.2d, v8.2d, v10.2d       \n\t"
        "mov v27.d[0], v11.d[1]         \n\t"
        "uqadd v7.2d, v7.2d, v27.2d     \n\t"        

        "mov v21.s[0], v11.s[0]         \n\t"   // c4

        "stp x26, x27, [%2, #80]     \n\t"
        
        "mul x13, x4, x10           \n\t"
        "umulh x14, x4, x10         \n\t"

        //step 6
        "umlal v12.2d, v4.2s, v3.s[1]    \n\t"   // v12 += a0b5 | a4b5
        "umlal v6.2d, v5.2s, v3.s[1]     \n\t"   // v6 += a1b5 | a5b5

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"
        
        "umlal2 v8.2d, v4.4s, v3.s[1]   \n\t"   // v8 += a2b5 | a6b5
        "umlal2 v7.2d, v5.4s, v3.s[1]   \n\t"   // v7 += a3b5 | a7b5  
        
        "mul x13, x5, x9            \n\t"
        "umulh x14, x5, x9          \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v11.16b, v11.16b, v11.16b  \n\t"

        "trn2 v0.4s, v12.4s, v0.4s      \n\t"   // v0 = h(a0b5 | a4b5)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a0b5 | a4b5)
        "trn2 v1.4s, v6.4s, v1.4s       \n\t"   // v1 = h(a1b5 | a5b5)
        "trn1 v6.4s, v6.4s, v28.4s      \n\t"   // v6 = l(a1b5 | a5b5)
        "trn2 v10.4s, v8.4s, v10.4s     \n\t"   // v10 = h(a2b5 | a6b5)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a2b5 | a6b5)
        "trn2 v11.4s, v7.4s, v11.4s     \n\t"   // v11 = h(a3b5 | a7b5)
        "trn1 v7.4s, v7.4s, v28.4s      \n\t"   // v7 = l(a3b5 | a7b5)

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "add v6.2d, v6.2d, v0.2d       \n\t"
        "add v8.2d, v8.2d, v1.2d        \n\t"
        "add v7.2d, v7.2d, v10.2d       \n\t"
        "mov v27.d[0], v12.d[1]         \n\t"
        "uqadd v11.2d, v11.2d, v27.2d   \n\t"        

        "mov v21.s[1], v12.s[0]         \n\t"   // c5
        
        "mul x13, x6, x8            \n\t"
        "umulh x14, x6, x8          \n\t"
        
        //step 7
        "umlal v6.2d, v4.2s, v3.s[2]    \n\t"   // v6 += a0b6 | a4b6
        "umlal v8.2d, v5.2s, v3.s[2]    \n\t"   // v8 += a1b6 | a5b6

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "umlal2 v7.2d, v4.4s, v3.s[2]   \n\t"   // v7 += a2b6 | a6b6
        "umlal2 v11.2d, v5.4s, v3.s[2]  \n\t"   // v11 += a3b6 | a7b6
        
        "mul x13, x5, x10           \n\t"
        "umulh x14, x5, x10         \n\t"
        
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v12.16b, v12.16b, v12.16b  \n\t"

        "trn2 v0.4s, v6.4s, v0.4s       \n\t"   // v0 = h(a0b6 | a4b6)
        "trn1 v6.4s, v6.4s, v28.4s      \n\t"   // v6 = l(a0b6 | a4b6)
        "trn2 v1.4s, v8.4s, v1.4s       \n\t"   // v1 = h(a1b6 | a5b6)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a1b6 | a5b6)
        "trn2 v10.4s, v7.4s, v10.4s     \n\t"   // v10 = h(a2b6 | a6b6)
        "trn1 v7.4s, v7.4s, v28.4s      \n\t"   // v7 = l(a2b6 | a6b6)
        "trn2 v12.4s, v11.4s, v12.4s    \n\t"   // v12 = h(a3b6 | a7b6)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a3b6 | a7b6)
      
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "add v8.2d, v8.2d, v0.2d       \n\t"
        "add v7.2d, v7.2d, v1.2d       \n\t"
        "add v11.2d, v11.2d, v10.2d    \n\t"
        "mov v27.d[0], v6.d[1]         \n\t"
        "uqadd v12.2d, v12.2d, v27.2d  \n\t"        

        "mov v21.s[2], v6.s[0]         \n\t"    // c6
        
        "mul x13, x6, x9            \n\t"
        "umulh x14, x6, x9          \n\t"

        // step 8
        "umlal v8.2d, v4.2s, v3.s[3]    \n\t"   // v8 += a0b7 | a4b7
        "umlal v7.2d, v5.2s, v3.s[3]    \n\t"   // v7 += a1b7 | a5b7

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "umlal2 v11.2d, v4.4s, v3.s[3]  \n\t"   // v11 += a2b7 | a6b7
        "umlal2 v12.2d, v5.4s, v3.s[3]  \n\t"   // v12 += a3b7 | a7b7

        "stp x25, x18, [%2, #96]     \n\t"
        
        "mul x13, x6, x10           \n\t"
        "umulh x14, x6, x10         \n\t"
       
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v6.16b, v6.16b, v6.16b     \n\t"

        "trn2 v0.4s, v8.4s, v0.4s       \n\t"   // v0 = h(a0b7 | a4b7)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a0b7 | a4b7)
        "trn2 v1.4s, v7.4s, v1.4s       \n\t"   // v1 = h(a1b7 | a5b7)
        "trn1 v7.4s, v7.4s, v28.4s      \n\t"   // v7 = l(a1b7 | a5b7)
        "trn2 v10.4s, v11.4s, v10.4s    \n\t"   // v10 = h(a2b7 | a6b7)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a2b7 | a6b7)
        "trn2 v6.4s, v12.4s, v6.4s      \n\t"   // v6 = h(a3b7 | a7b7)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a3b7 | a7b7)

        "adds x26, x26, x13         \n\t"
        "adc x27, x27, x14          \n\t"

        "add v7.2d, v7.2d, v0.2d       \n\t"
        "add v11.2d, v11.2d, v1.2d     \n\t"
        "add v12.2d, v12.2d, v10.2d    \n\t"
        "mov v27.d[0], v8.d[1]         \n\t"
        "uqadd v6.2d, v6.2d, v27.2d    \n\t"        

        "mov v21.s[3], v8.s[0]         \n\t"    // c7

        "stp x26, x27, [%2, #112]    \n\t"

        // final additions and carry propagations
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v8.16b, v8.16b, v8.16b     \n\t"

        "trn2 v0.4s, v7.4s, v0.4s       \n\t"   // v0 = h(c8 | c12)
        "trn1 v7.4s, v7.4s, v28.4s      \n\t"   // v7 = l(c8 | c12)
        "uqadd v11.2d, v11.2d, v0.2d   \n\t"        

        "mov v22.s[0], v7.s[0]         \n\t"    // c8

        "trn2 v1.4s, v11.4s, v1.4s      \n\t"   // v1 = h(c9 | c13)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(c9 | c13)
        "uqadd v12.2d, v12.2d, v1.2d   \n\t"        

        "mov v22.s[1], v11.s[0]         \n\t"    // c9


        "trn2 v10.4s, v12.4s, v10.4s    \n\t"   // v10 = h(c10 | c14)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(c10 | c14)
        "uqadd v6.2d, v6.2d, v10.2d     \n\t"        

        "mov v22.s[2], v12.s[0]         \n\t"    // c10


        "trn2 v8.4s, v6.4s, v8.4s       \n\t"   // v8 = h(c11 | c15)
        "trn1 v6.4s, v6.4s, v28.4s      \n\t"   // v6 = l(c11 | c15)
        "mov v27.d[1], v8.d[0]          \n\t"
        "uqadd v7.2d, v7.2d, v27.2d     \n\t"        

        "mov v22.s[3], v6.s[0]         \n\t"    // c11

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"
        "eor v8.16b, v8.16b, v8.16b     \n\t"

        "trn2 v0.4s, v7.4s, v0.4s       \n\t"   // v0 = h(c8 | c12)
        "trn1 v7.4s, v7.4s, v28.4s      \n\t"   // v7 = l(c8 | c12)
        "uqadd v11.2d, v11.2d, v0.2d   \n\t"        

        "mov v23.s[0], v7.s[2]         \n\t"    // c12
        
        "trn2 v1.4s, v11.4s, v1.4s      \n\t"   // v1 = h(c9 | c13)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(c9 | c13)
        "uqadd v12.2d, v12.2d, v1.2d   \n\t"        

        "mov v23.s[1], v11.s[2]         \n\t"    // c13

        "trn2 v10.4s, v12.4s, v10.4s    \n\t"   // v10 = h(c10 | c14)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(c10 | c14)
        "uqadd v6.2d, v6.2d, v10.2d     \n\t"        

        "mov v23.s[2], v12.s[2]         \n\t"    // c14

        "mov v23.s[3], v6.s[2]          \n\t"    // c15

        "st1 {v20.2d, v21.2d, v22.2d, v23.2d}, [%2]        \n\t"

        :
		:"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x18", "x20", "x21", "x22", "x25", "x26", "x27"
    ); 
}

void fpadd256_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{    //Add two 256-bit values
    asm(

	    "ldp x3, x4, [%0]           \n\t"
	    "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x7, x8, [%1]           \n\t"
        "ldp x9, x10, [%1, #16]      \n\t"
        "adds x3, x3, x7            \n\t"
        "adcs x4, x4, x8            \n\t"
        "adcs x5, x5, x9            \n\t"
        "adcs x6, x6, x10           \n\t"
        "eor x11, x11, x11         \n\t"
        "adc x12, x11, x11          \n\t"
        "stp x3, x4, [%2]           \n\t"
        "stp x5, x6, [%2, #16]      \n\t"
        "str x12, [%2, #32]         \n\t"
		:
		:"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
		:"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12"
    );
}

void fpadd576_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{   // Add two 576-bit values
    asm(

        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x7, x8, [%0, #32]      \n\t"
        "ldp x9, x10, [%0, #48]     \n\t"
        "ldr x20, [%0, #64]         \n\t"
        "ldp x11, x12, [%1]         \n\t"
        "ldp x13, x14, [%1, #16]    \n\t"
        "ldp x15, x16, [%1, #32]    \n\t"
        "ldp x17, x18, [%1, #48]    \n\t"
        "ldr x19, [%1, #64]         \n\t"
        "adds x3, x3, x11           \n\t"
        "adcs x4, x4, x12           \n\t"
        "adcs x5, x5, x13           \n\t"
        "adcs x6, x6, x14           \n\t"
        "adcs x7, x7, x15           \n\t"
        "adcs x8, x8, x16           \n\t"
        "adcs x9, x9, x17           \n\t"
        "adcs x10, x10, x18         \n\t"
        "adc x20, x20, x19          \n\t"
        "stp x3, x4, [%2]           \n\t"
        "stp x5, x6, [%2, #16]      \n\t"
        "stp x7, x8, [%2, #32]      \n\t"
        "stp x9, x10, [%2, #48]     \n\t"
        "str x20, [%2, #64]         \n\t"
        :
		:"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x15", "x16", "x17", "x18", "x19", "x20"
    );
}

void fpsub576_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{   // Sub two 576-bit values
    asm(

        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x7, x8, [%0, #32]      \n\t"
        "ldp x9, x10, [%0, #48]     \n\t"
        "ldr x19, [%0, #64]         \n\t"
    
        "ldp x11, x12, [%1]         \n\t"
        "ldp x13, x14, [%1, #16]    \n\t"
        "ldp x15, x16, [%1, #32]    \n\t"
        "ldp x17, x18, [%1, #48]    \n\t"
    
        "subs x3, x3, x11           \n\t"
        "sbcs x4, x4, x12           \n\t"
        "sbcs x5, x5, x13           \n\t"
        "sbcs x6, x6, x14           \n\t"
        "sbcs x7, x7, x15           \n\t"
        "sbcs x8, x8, x16           \n\t"
        "sbcs x9, x9, x17           \n\t"
        "sbcs x10, x10, x18         \n\t"
        "eor x18, x18, x18          \n\t"
        "sbc x19, x19, x18          \n\t"
    
        "stp x3, x4, [%2]           \n\t"
        "stp x5, x6, [%2, #16]      \n\t"
        "stp x7, x8, [%2, #32]      \n\t"
        "stp x9, x10, [%2, #48]     \n\t"
        "str x19, [%2, #64]         \n\t"

        :
		:"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x15", "x16", "x17", "x18", "x19", "x20"
    );
}

void fpmul320_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{   // multiplication of two 320-bit values
    asm(
        
        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldr x20, [%0, #32]         \n\t"
        "ldp x7, x8, [%1]           \n\t"
        "ldp x9, x10, [%1, #16]     \n\t"
        "ldr x21, [%1, #32]         \n\t"
        "eor x22, x22, x22          \n\t"
        "mul x11, x3, x7            \n\t"
        "umulh x12, x3, x7          \n\t"
        "mul x13, x3, x8            \n\t"
        "umulh x14, x3, x8          \n\t"
        "adds x25, x12, x13         \n\t"
        "adcs x26, x14, x22         \n\t"
        "mul x13, x4, x7            \n\t"
        "umulh x14, x4, x7          \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x22, x22         \n\t"
        "stp x11, x25, [%2]         \n\t"
        "mul x13, x3, x9            \n\t"
        "umulh x14, x3, x9          \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x22, x22         \n\t"
        "mul x13, x4, x8            \n\t"
        "umulh x14, x4, x8          \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x22         \n\t"
        "mul x13, x5, x7            \n\t"
        "umulh x14, x5, x7          \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x22         \n\t"
        "mul x13, x3, x10           \n\t"
        "umulh x14, x3, x10         \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x22, x22         \n\t"
        "mul x13, x4, x9            \n\t"
        "umulh x14, x4, x9          \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x22         \n\t"
        "mul x13, x5, x8            \n\t"
        "umulh x14, x5, x8          \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x22         \n\t"
        "mul x13, x6, x7            \n\t"
        "umulh x14, x6, x7          \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x22         \n\t"
        "stp x26, x27, [%2, #16]    \n\t"
        "mul x13, x4, x10           \n\t"
        "umulh x14, x4, x10         \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x22, x22         \n\t"
        "mul x13, x5, x9            \n\t"
        "umulh x14, x5, x9          \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x22         \n\t"
        "mul x13, x6, x8            \n\t"
        "umulh x14, x6, x8          \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x22         \n\t"
        "mul x13, x3, x21           \n\t"
        "umulh x14, x3, x21         \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x22         \n\t"
        "mul x13, x20, x7           \n\t"
        "umulh x14, x20, x7         \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x22         \n\t"
        "mul x13, x5, x10           \n\t"
        "umulh x14, x5, x10         \n\t"
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x22, x22         \n\t"
        "mul x13, x6, x9            \n\t"
        "umulh x14, x6, x9          \n\t"
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x22         \n\t"
        "mul x13, x4, x21           \n\t"
        "umulh x14, x4, x21         \n\t"
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x22         \n\t"
        "mul x13, x20, x8           \n\t"
        "umulh x14, x20, x8         \n\t"
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x22         \n\t"
        "stp x25, x18, [%2, #32]    \n\t"
        "mul x13, x6, x10           \n\t"
        "umulh x14, x6, x10         \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x22, x22         \n\t"
        "mul x13, x5, x21           \n\t"
        "umulh x14, x5, x21         \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x22         \n\t"
        "mul x13, x20, x9           \n\t"
        "umulh x14, x20, x9         \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x22         \n\t"
        "mul x13, x6, x21           \n\t"
        "umulh x14, x6, x21         \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x22, x22         \n\t"
        "mul x13, x20, x10          \n\t"
        "umulh x14, x20, x10        \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x22         \n\t"
        "stp x26, x27, [%2, #48]    \n\t"
        "mul x13, x20, x21          \n\t"
        "umulh x14, x20, x21        \n\t"
        "adds x25, x25, x13         \n\t"
        "adc x18, x18, x14          \n\t"
        "stp x25, x18, [%2, #64]    \n\t"
        :
		:"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x18", "x19", "x20", "x21", "x22", "x25", "x26", "x27"
    );
}

void fpmul256_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{   // Multiplication of two 256-bit values
    asm(
        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x7, x8, [%1]           \n\t"
        "ldp x9, x10, [%1, #16]     \n\t"
        "eor x20, x20, x20          \n\t"
        "mul x11, x3, x7            \n\t"
        "umulh x12, x3, x7          \n\t"
        "mul x13, x3, x8            \n\t"
        "umulh x14, x3, x8          \n\t"
        "adds x25, x12, x13         \n\t"
        "adcs x26, x14, x20         \n\t"
        "mul x13, x4, x7            \n\t"
        "umulh x14, x4, x7          \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"
        "stp x11, x25, [%2]         \n\t"
        "mul x13, x3, x9            \n\t"
        "umulh x14, x3, x9          \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"
        "mul x13, x4, x8            \n\t"
        "umulh x14, x4, x8          \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"
        "mul x13, x5, x7            \n\t"
        "umulh x14, x5, x7          \n\t"
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"
        "mul x13, x3, x10           \n\t"
        "umulh x14, x3, x10         \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x20, x20         \n\t"
        "mul x13, x4, x9            \n\t"
        "umulh x14, x4, x9          \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"
        "mul x13, x5, x8            \n\t"
        "umulh x14, x5, x8          \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"
        "mul x13, x6, x7            \n\t"
        "umulh x14, x6, x7          \n\t"
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"
        "stp x26, x27, [%2, #16]     \n\t"
        "mul x13, x4, x10           \n\t"
        "umulh x14, x4, x10         \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"
        "mul x13, x5, x9            \n\t"
        "umulh x14, x5, x9          \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        "mul x13, x6, x8            \n\t"
        "umulh x14, x6, x8          \n\t"
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        "mul x13, x5, x10           \n\t"
        "umulh x14, x5, x10         \n\t"
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"
        "mul x13, x6, x9            \n\t"
        "umulh x14, x6, x9          \n\t"
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"
        "stp x25, x18, [%2, #32]     \n\t"
        "mul x13, x6, x10           \n\t"
        "umulh x14, x6, x10         \n\t"
        "adds x26, x26, x13         \n\t"
        "adc x27, x27, x14          \n\t"
        "stp x26, x27, [%2, #48]    \n\t"
        :
		:"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x18", "x20", "x21", "x22", "x25", "x26", "x27"
    );
}

void fpadd512_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{   // Add two 512-bit values
    asm(

        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x7, x8, [%0, #32]      \n\t"
        "ldp x9, x10, [%0, #48]     \n\t"
        "ldp x11, x12, [%1]         \n\t"
        "ldp x13, x14, [%1, #16]    \n\t"
        "ldp x15, x16, [%1, #32]    \n\t"
        "ldp x17, x18, [%1, #48]    \n\t"
        "adds x3, x3, x11           \n\t"
        "adcs x4, x4, x12           \n\t"
        "adcs x5, x5, x13           \n\t"
        "adcs x6, x6, x14           \n\t"
        "adcs x7, x7, x15           \n\t"
        "adcs x8, x8, x16           \n\t"
        "adcs x9, x9, x17           \n\t"
        "adc x10, x10, x18         \n\t"
        "stp x3, x4, [%2]           \n\t"
        "stp x5, x6, [%2, #16]      \n\t"
        "stp x7, x8, [%2, #32]      \n\t"
        "stp x9, x10, [%2, #48]     \n\t"

        :
		:"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x15", "x16", "x17", "x18"
    );
}

void fpmul512_karatsuba(uint64_t *a, uint64_t *b, uint64_t *c)
{
    // level-two Karatsuba
    uint64_t aplus512[5];
    uint64_t bplus512[5];
    uint64_t rplus512[10];
    fpmul2x256_mixed_arm64(a, b, c);
    fpadd256_arm64(a, a+4, aplus512);
    fpadd256_arm64(b, b+4, bplus512);
    fpmul320_arm64(aplus512, bplus512, rplus512);
    fpsub576_arm64(rplus512, c, rplus512);
	fpsub576_arm64(rplus512, &c[8], rplus512);
    fpadd576_arm64(&c[4], rplus512, &c[4]); 
}

void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{
    // level-one Karatsuba
    uint64_t aplus[8];
    uint64_t bplus[8];
    uint64_t rplus[16];

    fpmul512_karatsuba(a, b, c);
    fpmul512_karatsuba(&a[8], &b[8], &c[16]);
	fpadd512_arm64(a, &a[8], aplus);
	fpadd512_arm64(b, &b[8], bplus);
	fpmul512_karatsuba(aplus, bplus, rplus);
	fpsub1024_arm64(rplus, c, rplus);
	fpsub1024_arm64(rplus, &c[16], rplus);
	fpadd1024_arm64(&c[8], rplus, &c[8]);
}

void rdc_mont(const uint64_t *ma, uint64_t *mc)
{ // Optimized Montgomery reduction using comba and exploiting the special form of the prime p964.
  // mc = ma*mb*R^-1 mod p964, where ma,mb,mc in [0, p964-1] and R = 2^1024.
  // ma and mb are assumed to be in Montgomery representation.
    rdc964_arm64_asm(ma, mc);
}
