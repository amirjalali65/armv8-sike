/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: modular arithmetic optimized for ARM64 platforms for P751
*********************************************************************************************/


#include "../P751_internal.h"

// Global constants
extern const uint64_t p751[NWORDS_FIELD];
extern const uint64_t p751x2[NWORDS_FIELD]; 

void fpmul2x384_mixed_arm64(uint64_t *a, uint64_t *b, uint64_t *c, uint64_t *aplusbplus)
{   
    // This function computes al * bl and ah * bh where a and b are 768-bit integers
    // The function computes the result of multiplication using mixed AARCH64 and ASIMD assembly 
    // operation:
    //      c[0...11] = a[0...5] * b[0...5] -> ASIMD assembly instructions
    //      c[12...23] = a[6...11] * b[6...11] -> AARCH64 assembly instructions
   asm volatile(
        "ldp x3, x4, [%0, #48]      \n\t"
        "ldp x5, x6, [%0, #64]      \n\t"
        "ldp x21, x22, [%0, #80]    \n\t"

        "ld1 {v0.2d, v1.2d, v2.2d}, [%0]        \n\t"
        "ld1 {v3.2d, v4.2d, v5.2d}, [%1]        \n\t"

        "ldp x7, x8, [%1, #48]      \n\t"
        "ldp x9, x10, [%1, #64]     \n\t"
        "ldp x23, x24, [%1, #80]    \n\t"

        // Shuffle a[0...11] -> a[0|4|2|6|1|5|3|7|8|10|9|11]
        "trn1 v6.4s, v0.4s, v1.4s       \n\t"
        "trn2 v7.4s, v0.4s, v1.4s       \n\t"
        "mov v8.s[0], v2.s[1]           \n\t"
        "mov v8.s[1], v2.s[2]           \n\t"
        "mov v2.s[1], v8.s[1]           \n\t"
        "mov v2.s[2], v8.s[0]           \n\t"
        // a in v6, v7, v2
        
        "eor x20, x20, x20          \n\t"   // use as zero register
        //step1
        "umull v8.2d, v6.2s, v3.s[0]    \n\t"   // v8 = a0b0 | a4b0
        "umull v9.2d, v7.2s, v3.s[0]    \n\t"   // v9 = a1b0 | a5b0

        "mul x11, x3, x7            \n\t"
        "umulh x12, x3, x7          \n\t"

        "umull2 v11.2d, v6.4s, v3.s[0]  \n\t"   // v11 = a2b0 | a6b0
        "umull2 v12.2d, v7.4s, v3.s[0]  \n\t"   // v12 = a3b0 | a7b0
       
        "mul x13, x3, x8            \n\t"
        "umulh x14, x3, x8          \n\t"

        "umull v10.2d, v2.2s, v3.s[0]   \n\t"   // v10 = a8b0 | a10b0
        "umull2 v13.2d, v2.4s, v3.s[0]  \n\t"   // v13 = a9b0 | a11b0

        "adds x25, x12, x13         \n\t"
        "adcs x26, x14, x20         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v16.16b, v16.16b, v16.16b  \n\t"
        
        
        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v17.16b, v17.16b, v17.16b  \n\t"

        "eor v28.16b, v28.16b, v28.16b  \n\t"   // for sepcial purpose
        "eor v27.16b, v27.16b, v27.16b  \n\t"   // for special purpose     
        "eor v26.16b, v26.16b, v26.16b  \n\t"   // for special purpose

        "trn2 v0.4s, v8.4s, v0.4s       \n\t"   // v0 = h(a0b0 | a4b0)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a0b0 | a4b0)
        "trn2 v1.4s, v9.4s, v1.4s       \n\t"   // v1 = h(a1b0 | a5b0)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a1b0 | a5b0)
        "trn2 v15.4s, v11.4s, v15.4s    \n\t"   // v15 = h(a2b0 | a6b0)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a2b0 | a6b0)
        "trn2 v16.4s, v12.4s, v16.4s    \n\t"   // v16 = h(a3b0 | a7b0)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a3b0 | a7b0)
        
        "trn2 v14.4s, v10.4s, v14.4s    \n\t"   // v14 = h(a8b0 | a10b0)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a8b0 | a10b0)
        "trn2 v17.4s, v13.4s, v17.4s    \n\t"   // v17 = h(a9b0 | a11b0)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a9b0 | a11b0)
        
        "mul x13, x4, x7            \n\t"
        "umulh x14, x4, x7          \n\t"

        "add v9.2d, v9.2d, v0.2d        \n\t"
        "add v11.2d, v11.2d, v1.2d      \n\t"
        "add v12.2d, v12.2d, v15.2d     \n\t"
        "add v13.2d, v13.2d, v14.2d     \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "mov v27.d[0], v8.d[1]          \n\t"
        "mov v27.d[1], v10.d[0]         \n\t"
        "uqadd v16.2d, v16.2d, v27.2d   \n\t"        

        "mov v26.d[0], v10.d[1]         \n\t"
        "uqadd v17.2d, v17.2d, v26.2d   \n\t"

        // final result is in v20, v21, v22, v23, v24, v25
        "mov v20.s[0], v8.s[0]          \n\t"   // c0
        
        "stp x11, x25, [%2, #96]    \n\t"   // c0, c1
        
        //step2
        "umlal v9.2d, v6.2s, v3.s[1]    \n\t"   // v9 += a0b1 | a4b1
        "umlal v11.2d, v7.2s, v3.s[1]   \n\t"   // v11 += a1b1 | a5b1

        "mul x13, x3, x9            \n\t"
        "umulh x14, x3, x9          \n\t"
        
        "umlal2 v12.2d, v6.4s, v3.s[1]  \n\t"   // v12 += a2b1 | a6b1
        "umlal2 v16.2d, v7.4s, v3.s[1]  \n\t"   // v16 += a3b1 | a7b1

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "umlal v13.2d, v2.2s, v3.s[1]   \n\t"   // v13 += a8b1 | a10b1
        "umlal2 v17.2d, v2.4s, v3.s[1]  \n\t"   // v17 += a9b1 | a11b1
        
        "mul x13, x4, x8            \n\t"
        "umulh x14, x4, x8          \n\t"
        
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v8.16b, v8.16b, v8.16b     \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"

        "trn2 v0.4s, v9.4s, v0.4s       \n\t"   // v0 = h(a0b1 | a4b1)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a0b1 | a4b1)
        "trn2 v1.4s, v11.4s, v1.4s      \n\t"   // v1 = h(a1b1 | a5b1)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a1b1 | a5b1)
        "trn2 v15.4s, v12.4s, v15.4s    \n\t"   // v15 = h(a2b1 | a6b1)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a2b1 | a6b1)
        "trn2 v8.4s, v16.4s, v8.4s      \n\t"   // v8 = h(a3b1 | a7b1)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a3b1 | a7b1)
        
        "trn2 v14.4s, v13.4s, v14.4s    \n\t"   // v14 = h(a8b1 | a10b1)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a8b1 | a10b1)
        "trn2 v10.4s, v17.4s, v10.4s    \n\t"   // v10 = h(a9b1 | a11b1)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a9b1 | a11b1)

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"
        
        "add v11.2d, v11.2d, v0.2d      \n\t"
        "add v12.2d, v12.2d, v1.2d      \n\t"
        "add v16.2d, v16.2d, v15.2d     \n\t"
        "add v17.2d, v17.2d, v14.2d     \n\t"

        "mul x13, x5, x7            \n\t"
        "umulh x14, x5, x7          \n\t"
        
        "mov v27.d[0], v9.d[1]          \n\t"
        "mov v27.d[1], v13.d[0]         \n\t"
        "uqadd v8.2d, v8.2d, v27.2d     \n\t"        

        "mov v26.d[0], v13.d[1]         \n\t"
        "uqadd v10.2d, v10.2d, v26.2d   \n\t"

        "mov v20.s[1], v9.s[0]          \n\t"   // c1

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"
        
        //step3
        "umlal v11.2d, v6.2s, v3.s[2]   \n\t"   // v9 += a0b2 | a4b2
        "umlal v12.2d, v7.2s, v3.s[2]   \n\t"   // v11 += a1b2 | a5b2

        "mul x13, x3, x10           \n\t"
        "umulh x14, x3, x10         \n\t"

        "umlal2 v16.2d, v6.4s, v3.s[2]  \n\t"   // v12 += a2b2 | a6b2
        "umlal2 v8.2d, v7.4s, v3.s[2]   \n\t"   // v16 += a3b2 | a7b2
        
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x20, x20         \n\t"
        
        "umlal v17.2d, v2.2s, v3.s[2]   \n\t"   // v13 += a8b2 | a10b2
        "umlal2 v10.2d, v2.4s, v3.s[2]  \n\t"   // v17 += a9b2 | a11b2

        "mul x13, x4, x9            \n\t"
        "umulh x14, x4, x9          \n\t"
        
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v9.16b, v9.16b, v9.16b     \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v13.16b, v13.16b, v13.16b  \n\t"

        "trn2 v0.4s, v11.4s, v0.4s      \n\t"   // v0 = h(a0b2 | a4b2)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a0b2 | a4b2)
        "trn2 v1.4s, v12.4s, v1.4s      \n\t"   // v1 = h(a1b2 | a5b2)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a1b2 | a5b2)
        "trn2 v15.4s, v16.4s, v15.4s    \n\t"   // v15 = h(a2b2 | a6b2)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a2b2 | a6b2)
        "trn2 v9.4s, v8.4s, v9.4s       \n\t"   // v9 = h(a3b2 | a7b2)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a3b2 | a7b2)
        
        "trn2 v14.4s, v17.4s, v14.4s    \n\t"   // v14 = h(a8b2 | a10b2)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a8b2 | a10b2)
        "trn2 v13.4s, v10.4s, v13.4s    \n\t"   // v13 = h(a9b2 | a11b2)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a9b2 | a11b2)

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"
        
        "add v12.2d, v12.2d, v0.2d      \n\t"
        "add v16.2d, v16.2d, v1.2d      \n\t"
        "add v8.2d, v8.2d, v15.2d       \n\t"
        "add v10.2d, v10.2d, v14.2d     \n\t"

        "mul x13, x5, x8            \n\t"
        "umulh x14, x5, x8          \n\t"
        
        "mov v27.d[0], v11.d[1]          \n\t"
        "mov v27.d[1], v17.d[0]         \n\t"
        "uqadd v9.2d, v9.2d, v27.2d     \n\t"        

        "mov v26.d[0], v17.d[1]         \n\t"
        "uqadd v13.2d, v13.2d, v26.2d   \n\t"

        "mov v20.s[2], v11.s[0]         \n\t"   // c2

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"
        
        //step4
        "umlal v12.2d, v6.2s, v3.s[3]   \n\t"   // v12 += a0b3 | a4b3
        "umlal v16.2d, v7.2s, v3.s[3]   \n\t"   // v16 += a1b3 | a5b3

        "mul x13, x6, x7            \n\t"
        "umulh x14, x6, x7          \n\t"

        "umlal2 v8.2d, v6.4s, v3.s[3]   \n\t"   // v8 += a2b3 | a6b3
        "umlal2 v9.2d, v7.4s, v3.s[3]   \n\t"   // v9 += a3b3 | a7b3        

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "umlal v10.2d, v2.2s, v3.s[3]   \n\t"   // v10 += a8b3 | a10b3
        "umlal2 v13.2d, v2.4s, v3.s[3]  \n\t"   // v13 += a9b3 | a11b3

        "stp x26, x27, [%2, #112]   \n\t"  // c2, c3

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v11.16b, v11.16b, v11.16b     \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v17.16b, v17.16b, v17.16b  \n\t"

        "trn2 v0.4s, v12.4s, v0.4s      \n\t"   // v0 = h(a0b3 | a4b3)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a0b3 | a4b3)
        "trn2 v1.4s, v16.4s, v1.4s      \n\t"   // v1 = h(a1b3 | a5b3)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a1b3 | a5b3)
        "trn2 v15.4s, v8.4s, v15.4s     \n\t"   // v15 = h(a2b3 | a6b3)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a2b3 | a6b3)
        "trn2 v11.4s, v9.4s, v11.4s     \n\t"   // v11 = h(a3b3 | a7b3)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a3b3 | a7b3)
        
        "trn2 v14.4s, v10.4s, v14.4s    \n\t"   // v14 = h(a8b3 | a10b3)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a8b3 | a10b3)
        "trn2 v17.4s, v13.4s, v17.4s    \n\t"   // v17 = h(a9b3 | a11b3)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a9b3 | a11b3)

        "mul x13, x3, x23           \n\t"
        "umulh x14, x3, x23         \n\t"

        "add v16.2d, v16.2d, v0.2d      \n\t"
        "add v8.2d, v8.2d, v1.2d        \n\t"
        "add v9.2d, v9.2d, v15.2d       \n\t"
        "add v13.2d, v13.2d, v14.2d     \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"

        "mov v27.d[0], v12.d[1]         \n\t"
        "mov v27.d[1], v10.d[0]         \n\t"
        "uqadd v11.2d, v11.2d, v27.2d   \n\t"        

        "mov v26.d[0], v10.d[1]         \n\t"
        "uqadd v17.2d, v17.2d, v26.2d   \n\t"

        "mov v20.s[3], v12.s[0]         \n\t"   // c3
        
        //step5
        "umlal v16.2d, v6.2s, v4.s[0]   \n\t"   // v16 += a0b4 | a4b4
        "umlal v8.2d, v7.2s, v4.s[0]    \n\t"   // v8 += a1b4 | a5b4

        "mul x13, x7, x21           \n\t"
        "umulh x14, x7, x21         \n\t"

        "umlal2 v9.2d, v6.4s, v4.s[0]   \n\t"   // v9 += a2b4 | a6b4
        "umlal2 v11.2d, v7.4s, v4.s[0]  \n\t"   // v11 += a3b4 | a7b4

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "umlal v13.2d, v2.2s, v4.s[0]   \n\t"   // v13 += a8b4 | a10b4
        "umlal2 v17.2d, v2.4s, v4.s[0]  \n\t"   // v17 += a9b4 | a11b4

        "mul x13, x4, x10           \n\t"
        "umulh x14, x4, x10         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v12.16b, v12.16b, v12.16b  \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"

        "trn2 v0.4s, v16.4s, v0.4s      \n\t"   // v0 = h(a0b4 | a4b4)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a0b4 | a4b4)
        "trn2 v1.4s, v8.4s, v1.4s       \n\t"   // v1 = h(a1b4 | a5b4)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a1b4 | a5b4)
        "trn2 v15.4s, v9.4s, v15.4s     \n\t"   // v15 = h(a2b4 | a6b4)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a2b4 | a6b4)
        "trn2 v12.4s, v11.4s, v12.4s    \n\t"   // v12 = h(a3b4 | a7b4)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a3b4 | a7b4)
        
        "trn2 v14.4s, v13.4s, v14.4s    \n\t"   // v14 = h(a8b4 | a10b4)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a8b4 | a10b4)
        "trn2 v10.4s, v17.4s, v10.4s    \n\t"   // v10 = h(a9b4 | a11b4)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a9b4 | a11b4)
        
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        
        "add v8.2d, v8.2d, v0.2d        \n\t"
        "add v9.2d, v9.2d, v1.2d        \n\t"
        "add v11.2d, v11.2d, v15.2d     \n\t"
        "add v17.2d, v17.2d, v14.2d     \n\t"

        "mul x13, x5, x9            \n\t"
        "umulh x14, x5, x9          \n\t"

        "mov v27.d[0], v16.d[1]         \n\t"
        "mov v27.d[1], v13.d[0]         \n\t"
        "uqadd v12.2d, v12.2d, v27.2d   \n\t"        

        "mov v26.d[0], v13.d[1]         \n\t"
        "uqadd v10.2d, v10.2d, v26.2d   \n\t"

        "mov v21.s[0], v16.s[0]         \n\t"   // c4
        
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        //step6
        "umlal v8.2d, v6.2s, v4.s[1]    \n\t"   // v8 += a0b5 | a4b5
        "umlal v9.2d, v7.2s, v4.s[1]    \n\t"   // v9 += a1b5 | a5b5
        
        "mul x13, x6, x8            \n\t"
        "umulh x14, x6, x8          \n\t"

        "umlal2 v11.2d, v6.4s, v4.s[1]  \n\t"   // v11 += a2b5 | a6b5
        "umlal2 v12.2d, v7.4s, v4.s[1]  \n\t"   // v12 += a3b5 | a7b5
        
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "umlal v17.2d, v2.2s, v4.s[1]   \n\t"   // v17 += a8b5 | a10b5
        "umlal2 v10.2d, v2.4s, v4.s[1]  \n\t"   // v10 += a9b5 | a11b5

        "mul x13, x3, x24           \n\t"
        "umulh x14, x3, x24         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v16.16b, v16.16b, v16.16b  \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v13.16b, v13.16b, v13.16b  \n\t"

        "trn2 v0.4s, v8.4s, v0.4s       \n\t"   // v0 = h(a0b5 | a4b5)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a0b5 | a4b5)
        "trn2 v1.4s, v9.4s, v1.4s       \n\t"   // v1 = h(a1b5 | a5b5)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a1b5 | a5b5)
        "trn2 v15.4s, v11.4s, v15.4s    \n\t"   // v15 = h(a2b5 | a6b5)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a2b5 | a6b5)
        "trn2 v16.4s, v12.4s, v16.4s    \n\t"   // v16 = h(a3b5 | a7b5)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a3b5 | a7b5)
        
        "trn2 v14.4s, v17.4s, v14.4s    \n\t"   // v14 = h(a8b5 | a10b5)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a8b5 | a10b5)
        "trn2 v13.4s, v10.4s, v13.4s    \n\t"   // v13 = h(a9b5 | a11b5)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a9b5 | a11b5)

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "add v9.2d, v9.2d, v0.2d        \n\t"
        "add v11.2d, v11.2d, v1.2d      \n\t"
        "add v12.2d, v12.2d, v15.2d     \n\t"
        "add v10.2d, v10.2d, v14.2d     \n\t"

        "mul x13, x7, x22           \n\t"
        "umulh x14, x7, x22         \n\t"

        "mov v27.d[0], v8.d[1]          \n\t"
        "mov v27.d[1], v17.d[0]         \n\t"
        "uqadd v16.2d, v16.2d, v27.2d   \n\t"        

        "mov v26.d[0], v17.d[1]         \n\t"
        "uqadd v13.2d, v13.2d, v26.2d   \n\t"

        "mov v21.s[1], v8.s[0]          \n\t"   // c5
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        //step7
        "umlal v9.2d, v6.2s, v4.s[2]    \n\t"   // v9 += a0b6 | a4b6
        "umlal v11.2d, v7.2s, v4.s[2]   \n\t"   // v11 += a1b6 | a5b6

        "mul x13, x4, x23           \n\t"
        "umulh x14, x4, x23         \n\t"

        "umlal2 v12.2d, v6.4s, v4.s[2]  \n\t"   // v12 += a2b6 | a6b6
        "umlal2 v16.2d, v7.4s, v4.s[2]  \n\t"   // v16 += a3b6 | a7b6
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "umlal v10.2d, v2.2s, v4.s[2]   \n\t"   // v10 += a8b6 | a10b6
        "umlal2 v13.2d, v2.4s, v4.s[2]  \n\t"   // v13 += a9b6 | a11b6

        "mul x13, x8, x21           \n\t"
        "umulh x14, x8, x21         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v8.16b, v8.16b, v8.16b  \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v17.16b, v17.16b, v17.16b  \n\t"

        "trn2 v0.4s, v9.4s, v0.4s       \n\t"   // v0 = h(a0b6 | a4b6)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a0b6 | a4b6)
        "trn2 v1.4s, v11.4s, v1.4s      \n\t"   // v1 = h(a1b6 | a5b6)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a1b6 | a5b6)
        "trn2 v15.4s, v12.4s, v15.4s    \n\t"   // v15 = h(a2b6 | a6b6)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a2b6 | a6b6)
        "trn2 v8.4s, v16.4s, v8.4s      \n\t"   // v8 = h(a3b6 | a7b6)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a3b6 | a7b6)
        
        "trn2 v14.4s, v10.4s, v14.4s    \n\t"   // v14 = h(a8b6 | a10b6)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a8b6 | a10b6)
        "trn2 v17.4s, v13.4s, v17.4s    \n\t"   // v17 = h(a9b6 | a11b6)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a9b6 | a11b6)
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "add v11.2d, v11.2d, v0.2d      \n\t"
        "add v12.2d, v12.2d, v1.2d      \n\t"
        "add v16.2d, v16.2d, v15.2d     \n\t"
        "add v13.2d, v13.2d, v14.2d     \n\t"

        "mul x13, x5, x10           \n\t"
        "umulh x14, x5, x10         \n\t"

        "mov v27.d[0], v9.d[1]          \n\t"
        "mov v27.d[1], v10.d[0]         \n\t"
        "uqadd v8.2d, v8.2d, v27.2d     \n\t"        

        "mov v26.d[0], v10.d[1]         \n\t"
        "uqadd v17.2d, v17.2d, v26.2d   \n\t"

        "mov v21.s[2], v9.s[0]          \n\t"   // c6
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        //step8
        "umlal v11.2d, v6.2s, v4.s[3]   \n\t"   // v11 += a0b7 | a4b7
        "umlal v12.2d, v7.2s, v4.s[3]   \n\t"   // v12 += a1b7 | a5b7
        
        "mul x13, x6, x9            \n\t"
        "umulh x14, x6, x9          \n\t"

        "umlal2 v16.2d, v6.4s, v4.s[3]  \n\t"   // v16 += a2b7 | a6b7
        "umlal2 v8.2d, v7.4s, v4.s[3]   \n\t"   // v8 += a3b7 | a7b7
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "umlal v13.2d, v2.2s, v4.s[3]   \n\t"   // v13 += a8b7 | a10b7
        "umlal2 v17.2d, v2.4s, v4.s[3]  \n\t"   // v17 += a9b7 | a11b7
        
        "stp x25, x18, [%2, #128]   \n\t"  // c4, c5

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v9.16b, v9.16b, v9.16b  \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"

        "trn2 v0.4s, v11.4s, v0.4s      \n\t"   // v0 = h(a0b7 | a4b7)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a0b7 | a4b7)
        "trn2 v1.4s, v12.4s, v1.4s      \n\t"   // v1 = h(a1b7 | a5b7)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a1b7 | a5b7)
        "trn2 v15.4s, v16.4s, v15.4s    \n\t"   // v15 = h(a2b7 | a6b7)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a2b7 | a6b7)
        "trn2 v9.4s, v8.4s, v9.4s       \n\t"   // v9 = h(a3b7 | a7b7)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a3b7 | a7b7)
        
        "trn2 v14.4s, v13.4s, v14.4s    \n\t"   // v14 = h(a8b7 | a10b7)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a8b7 | a10b7)
        "trn2 v10.4s, v17.4s, v10.4s    \n\t"   // v10 = h(a9b7 | a11b7)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a9b7 | a11b7)

        "mul x13, x4, x24           \n\t"
        "umulh x14, x4, x24         \n\t"

        "add v12.2d, v12.2d, v0.2d      \n\t"
        "add v16.2d, v16.2d, v1.2d      \n\t"
        "add v8.2d, v8.2d, v15.2d       \n\t"
        "add v17.2d, v17.2d, v14.2d     \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "mov v27.d[0], v11.d[1]         \n\t"
        "mov v27.d[1], v13.d[0]         \n\t"
        "uqadd v9.2d, v9.2d, v27.2d    \n\t"        

        "mov v26.d[0], v13.d[1]         \n\t"
        "uqadd v10.2d, v10.2d, v26.2d   \n\t"

        "mov v21.s[3], v11.s[0]         \n\t"   // c7

        "mul x13, x5, x23           \n\t"
        "umulh x14, x5, x23         \n\t"

        //step9
        "umlal v12.2d, v6.2s, v5.s[0]   \n\t"   // v12 += a0b8 | a4b8
        "umlal v16.2d, v7.2s, v5.s[0]   \n\t"   // v16 += a1b8 | a5b8

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "umlal2 v8.2d, v6.4s, v5.s[0]   \n\t"   // v8 += a2b8 | a6b8
        "umlal2 v9.2d, v7.4s, v5.s[0]   \n\t"   // v9 += a3b8 | a7b8

        "mul x13, x6, x10           \n\t"
        "umulh x14, x6, x10         \n\t"

        "umlal v17.2d, v2.2s, v5.s[0]   \n\t"   // v17 += a8b8 | a10b8
        "umlal2 v10.2d, v2.4s, v5.s[0]  \n\t"   // v10 += a9b8 | a11b8

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v11.16b, v11.16b, v11.16b  \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v13.16b, v13.16b, v13.16b  \n\t"

        "trn2 v0.4s, v12.4s, v0.4s      \n\t"   // v0 = h(a0b8 | a4b8)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a0b8 | a4b8)
        "trn2 v1.4s, v16.4s, v1.4s      \n\t"   // v1 = h(a1b8 | a5b8)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a1b8 | a5b8)
        "trn2 v15.4s, v8.4s, v15.4s     \n\t"   // v15 = h(a2b8 | a6b8)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a2b8 | a6b8)
        "trn2 v11.4s, v9.4s, v11.4s     \n\t"   // v11 = h(a3b8 | a7b8)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a3b8 | a7b8)
        
        "trn2 v14.4s, v17.4s, v14.4s    \n\t"   // v14 = h(a8b8 | a10b8)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a8b8 | a10b8)
        "trn2 v13.4s, v10.4s, v13.4s    \n\t"   // v13 = h(a9b8 | a11b8)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a9b8 | a11b8)

        "mul x13, x21, x9           \n\t"
        "umulh x14, x21, x9         \n\t"

        "add v16.2d, v16.2d, v0.2d      \n\t"
        "add v8.2d, v8.2d, v1.2d        \n\t"
        "add v9.2d, v9.2d, v15.2d       \n\t"
        "add v10.2d, v10.2d, v14.2d     \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mov v27.d[0], v12.d[1]         \n\t"
        "mov v27.d[1], v17.d[0]         \n\t"
        "uqadd v11.2d, v11.2d, v27.2d   \n\t"        

        "mov v26.d[0], v17.d[1]         \n\t"
        "uqadd v13.2d, v13.2d, v26.2d   \n\t"

        "mov v22.s[0], v12.s[0]         \n\t"   // c8

        "mul x13, x22, x8           \n\t"
        "umulh x14, x22, x8         \n\t"

        //step10
        "umlal v16.2d, v6.2s, v5.s[1]   \n\t"   // v16 += a0b9 | a4b9
        "umlal v8.2d, v7.2s, v5.s[1]    \n\t"   // v8 += a1b9 | a5b9

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "umlal2 v9.2d, v6.4s, v5.s[1]   \n\t"   // v9 += a2b9 | a6b9
        "umlal2 v11.2d, v7.4s, v5.s[1]  \n\t"   // v11 += a3b9 | a7b9

        "mul x13, x5, x24           \n\t"
        "umulh x14, x5, x24         \n\t"

        "umlal v10.2d, v2.2s, v5.s[1]   \n\t"   // v10 += a8b9 | a10b9
        "umlal2 v13.2d, v2.4s, v5.s[1]  \n\t"   // v13 += a9b9 | a11b9

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x20, x20         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v12.16b, v12.16b, v12.16b  \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v17.16b, v17.16b, v17.16b  \n\t"

        "trn2 v0.4s, v16.4s, v0.4s      \n\t"   // v0 = h(a0b9 | a4b9)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a0b9 | a4b9)
        "trn2 v1.4s, v8.4s, v1.4s       \n\t"   // v1 = h(a1b9 | a5b9)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a1b9 | a5b9)
        "trn2 v15.4s, v9.4s, v15.4s     \n\t"   // v15 = h(a2b9 | a6b9)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a2b9 | a6b9)
        "trn2 v12.4s, v11.4s, v12.4s    \n\t"   // v12 = h(a3b9 | a7b9)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a3b9 | a7b9)
        
        "trn2 v14.4s, v10.4s, v14.4s    \n\t"   // v14 = h(a8b9 | a10b9)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a8b9 | a10b9)
        "trn2 v17.4s, v13.4s, v17.4s    \n\t"   // v17 = h(a9b9 | a11b9)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a9b9 | a11b9)

        "mul x13, x6, x23           \n\t"
        "umulh x14, x6, x23         \n\t"

        "add v8.2d, v8.2d, v0.2d        \n\t"
        "add v9.2d, v9.2d, v1.2d        \n\t"
        "add v11.2d, v11.2d, v15.2d     \n\t"
        "add v13.2d, v13.2d, v14.2d     \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mov v27.d[0], v16.d[1]         \n\t"
        "mov v27.d[1], v10.d[0]         \n\t"
        "uqadd v12.2d, v12.2d, v27.2d   \n\t"        

        "mov v26.d[0], v10.d[1]         \n\t"
        "uqadd v17.2d, v17.2d, v26.2d   \n\t"

        "mov v22.s[1], v16.s[0]         \n\t"   // c9

        "mul x13, x21, x10           \n\t"
        "umulh x14, x21, x10         \n\t"

        //step11
        "umlal v8.2d, v6.2s, v5.s[2]    \n\t"   // v8 += a0b10 | a4b10
        "umlal v9.2d, v7.2s, v5.s[2]    \n\t"   // v9 += a1b10 | a5b10

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "umlal2 v11.2d, v6.4s, v5.s[2]  \n\t"   // v11 += a2b10 | a6b10
        "umlal2 v12.2d, v7.4s, v5.s[2]  \n\t"   // v12 += a3b10 | a7b10

        "mul x13, x22, x9           \n\t"
        "umulh x14, x22, x9         \n\t"

        "umlal v13.2d, v2.2s, v5.s[2]   \n\t"   // v13 += a8b10 | a10b10
        "umlal2 v17.2d, v2.4s, v5.s[2]  \n\t"   // v17 += a9b10 | a11b10

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v16.16b, v16.16b, v16.16b  \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v10.16b, v10.16b, v10.16b  \n\t"

        "trn2 v0.4s, v8.4s, v0.4s       \n\t"   // v0 = h(a0b10 | a4b10)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a0b10 | a4b10)
        "trn2 v1.4s, v9.4s, v1.4s       \n\t"   // v1 = h(a1b10 | a5b10)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a1b10 | a5b10)
        "trn2 v15.4s, v11.4s, v15.4s    \n\t"   // v15 = h(a2b10 | a6b10)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a2b10 | a6b10)
        "trn2 v16.4s, v12.4s, v16.4s    \n\t"   // v16 = h(a3b10 | a7b10)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a3b10 | a7b10)
        
        "trn2 v14.4s, v13.4s, v14.4s    \n\t"   // v14 = h(a8b10 | a10b10)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a8b10 | a10b10)
        "trn2 v10.4s, v17.4s, v10.4s    \n\t"   // v10 = h(a9b10 | a11b10)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a9b10 | a11b10)

        "stp x26, x27, [%2, #144]   \n\t"   // c6, c7

        "add v9.2d, v9.2d, v0.2d        \n\t"
        "add v11.2d, v11.2d, v1.2d      \n\t"
        "add v12.2d, v12.2d, v15.2d     \n\t"
        "add v17.2d, v17.2d, v14.2d     \n\t"

        "mul x13, x6, x24           \n\t"
        "umulh x14, x6, x24         \n\t"

        "mov v27.d[0], v8.d[1]          \n\t"
        "mov v27.d[1], v13.d[0]         \n\t"
        "uqadd v16.2d, v16.2d, v27.2d   \n\t"        

        "mov v26.d[0], v13.d[1]         \n\t"
        "uqadd v10.2d, v10.2d, v26.2d   \n\t"

        "mov v22.s[2], v8.s[0]          \n\t"   // c10

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"

        //step12
        "umlal v9.2d, v6.2s, v5.s[3]    \n\t"   // v9 += a0b11 | a4b11
        "umlal v11.2d, v7.2s, v5.s[3]   \n\t"   // v11 += a1b11 | a5b11

        "mul x13, x21, x23           \n\t"
        "umulh x14, x21, x23         \n\t"

        "umlal2 v12.2d, v6.4s, v5.s[3]  \n\t"   // v12 += a2b11 | a6b11
        "umlal2 v16.2d, v7.4s, v5.s[3]  \n\t"   // v16 += a3b11 | a7b11

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "umlal v17.2d, v2.2s, v5.s[3]   \n\t"   // v17 += a8b11 | a10b11
        "umlal2 v10.2d, v2.4s, v5.s[3]  \n\t"   // v10 += a9b11 | a11b11

        "mul x13, x22, x10           \n\t"
        "umulh x14, x22, x10         \n\t"

        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v8.16b, v8.16b, v8.16b     \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v13.16b, v13.16b, v13.16b  \n\t"

        "trn2 v0.4s, v9.4s, v0.4s       \n\t"   // v0 = h(a0b11 | a4b11)
        "trn1 v9.4s, v9.4s, v28.4s      \n\t"   // v9 = l(a0b11 | a4b11)
        "trn2 v1.4s, v11.4s, v1.4s      \n\t"   // v1 = h(a1b11 | a5b11)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a1b11 | a5b11)
        "trn2 v15.4s, v12.4s, v15.4s    \n\t"   // v15 = h(a2b11 | a6b11)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a2b11 | a6b11)
        "trn2 v8.4s, v16.4s, v8.4s      \n\t"   // v8 = h(a3b11 | a7b11)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a3b11 | a7b11)
        
        "trn2 v14.4s, v17.4s, v14.4s    \n\t"   // v14 = h(a8b11 | a10b11)
        "trn1 v17.4s, v17.4s, v28.4s    \n\t"   // v17 = l(a8b11 | a10b11)
        "trn2 v13.4s, v10.4s, v13.4s    \n\t"   // v13 = h(a9b11 | a11b11)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a9b11 | a11b11)

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "add v11.2d, v11.2d, v0.2d      \n\t"
        "add v12.2d, v12.2d, v1.2d      \n\t"
        "add v16.2d, v16.2d, v15.2d     \n\t"
        "add v10.2d, v10.2d, v14.2d     \n\t"

        "mul x13, x21, x24           \n\t"
        "umulh x14, x21, x24         \n\t"

        "mov v27.d[0], v9.d[1]          \n\t"
        "mov v27.d[1], v17.d[0]         \n\t"
        "uqadd v8.2d, v8.2d, v27.2d     \n\t"        

        "mov v26.d[0], v17.d[1]         \n\t"
        "uqadd v13.2d, v13.2d, v26.2d   \n\t"

        "mov v22.s[3], v9.s[0]         \n\t"   // c11

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        // final additions and carry propagations
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v9.16b, v9.16b, v9.16b     \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v17.16b, v17.16b, v17.16b  \n\t"

        "trn2 v0.4s, v11.4s, v0.4s      \n\t"   // v0 = h(a0b11 | a4b11)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a0b11 | a4b11)
        "uqadd v12.2d, v12.2d, v0.2d    \n\t"

        "mov v23.s[0], v11.s[0]         \n\t"  // c12

        "mul x13, x22, x23           \n\t"
        "umulh x14, x22, x23         \n\t"

        "trn2 v1.4s, v12.4s, v1.4s      \n\t"   // v1 = h(a1b11 | a5b11)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"   // v12 = l(a1b11 | a5b11)
        "uqadd v16.2d, v16.2d, v1.2d    \n\t"

        "mov v23.s[1], v12.s[0]         \n\t"   // c13

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "trn2 v15.4s, v16.4s, v15.4s    \n\t"   // v15 = h(a2b11 | a6b11)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a2b11 | a6b11)
        "uqadd v8.2d, v8.2d, v15.2d     \n\t"
        
        "mov v23.s[2], v16.s[0]         \n\t"   // c14
        
        "trn2 v9.4s, v8.4s, v9.4s       \n\t"   // v9 = h(a3b11 | a7b11)
        "trn1 v8.4s, v8.4s, v28.4s      \n\t"   // v8 = l(a3b11 | a7b11)
        "mov v26.d[1], v9.d[0]          \n\t"
        
        "uqadd v11.2d, v11.2d, v26.2d   \n\t"

        "stp x25, x18, [%2, #160]       \n\t"   // c8, c9

        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "mov v15.d[0], v9.d[1]          \n\t"
        "uqadd v10.2d, v10.2d, v15.2d   \n\t"

        "mov v23.s[3], v8.s[0]          \n\t"    // c15
        
        "trn2 v14.4s, v11.4s, v14.4s    \n\t"   // v14 = h(a0b11 | a4b11)
        "trn1 v11.4s, v11.4s, v28.4s    \n\t"   // v11 = l(a0b11 | a4b11)
        "uqadd v12.2d, v12.2d, v14.2d   \n\t"

        "mov v24.s[0], v11.s[2]         \n\t"   //c16

        "mul x13, x22, x24              \n\t"
        "umulh x14, x22, x24            \n\t"

        "trn2 v17.4s, v12.4s, v17.4s    \n\t"   // v17 = h(a1b11 | a5b11)
        "trn1 v12.4s, v12.4s, v28.4s    \n\t"     // v12 = l(a1b11 | a5b11)
        "uqadd v16.2d, v16.2d, v17.2d   \n\t"

        "mov v24.s[1], v12.s[2]         \n\t"   // c17
        
        "eor v0.16b, v0.16b, v0.16b     \n\t"
        "eor v1.16b, v1.16b, v1.16b     \n\t"
        "eor v15.16b, v15.16b, v15.16b  \n\t"
        "eor v9.16b, v9.16b, v9.16b     \n\t"

        "eor v14.16b, v14.16b, v14.16b  \n\t"
        "eor v17.16b, v17.16b, v17.16b  \n\t"
        
        "trn2 v0.4s, v16.4s, v0.4s      \n\t"   // v0 = h(a2b11 | a6b11)
        "trn1 v16.4s, v16.4s, v28.4s    \n\t"   // v16 = l(a2b11 | a6b11)
        "uqadd v8.2d, v8.2d, v0.2d      \n\t"
        
        "mov v24.s[2], v16.s[2]         \n\t"   // c18

        "adds x26, x26, x13             \n\t"
        "adc x27, x27, x14              \n\t"

        "trn2 v1.4s, v8.4s, v1.4s      \n\t"   // v1 = h(a3b11 | a7b11)
        "trn1 v8.4s, v8.4s, v28.4s     \n\t"   // v8 = l(a3b11 | a7b11)
        
        "eor v26.16b, v26.16b, v26.16b  \n\t"
        "mov v26.d[0], v1.d[1]         \n\t"
        "uqadd v10.2d, v10.2d, v26.2d  \n\t"

        "mov v24.s[3], v8.s[2]         \n\t"    // c19
        
        "trn2 v15.4s, v10.4s, v15.4s    \n\t"   // v15 = h(a8b11 | a10b11)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a8b11 | a10b11)
        "uqadd v13.2d, v13.2d, v15.2d   \n\t"

        "mov v25.s[0], v10.s[0]         \n\t"   // c20 *

        "stp x26, x27, [%2, #176]    \n\t"   // c8, c9

        "trn2 v9.4s, v13.4s, v9.4s    \n\t"   // v9 = h(a9b11 | a11b11)
        "trn1 v13.4s, v13.4s, v28.4s    \n\t"   // v13 = l(a9b11 | a11b11)
        
        "eor v26.16b, v26.16b, v26.16b  \n\t"
        "mov v26.d[1], v9.d[0]          \n\t"
        "uqadd v10.2d, v10.2d, v26.2d   \n\t"

        "mov v25.s[1], v13.s[0]         \n\t"   // c21

        "trn2 v14.4s, v10.4s, v14.4s    \n\t"   // v14 = h(a8b11 | a10b11)
        "trn1 v10.4s, v10.4s, v28.4s    \n\t"   // v10 = l(a8b11 | a10b11)
        "uqadd v13.2d, v13.2d, v14.2d   \n\t"

        "mov v25.s[2], v10.s[2]         \n\t"   // c22

        "mov v25.s[3], v13.s[2]         \n\t"    // c23
        
        "st1 {v20.2d, v21.2d}, [%2], #32    \n\t"
        "st1 {v22.2d, v23.2d, v24.2d, v25.2d}, [%2]       \n\t"
        // computing additions 
        "ldp x11, x12, [%0]         \n\t"
        "ldp x13, x14, [%0, #16]    \n\t"
        "ldp x25, x26, [%0, #32]    \n\t"

        "adds x3, x3, x11           \n\t"
        "adcs x4, x4, x12           \n\t"
        "adcs x5, x5, x13           \n\t"
        "adcs x6, x6, x14           \n\t"
        "adcs x21, x21, x25         \n\t"
        "adcs x22, x22, x26         \n\t"
        "eor x15, x15, x15          \n\t"
        "adc x15, x15, x15          \n\t"

        "ldp x11, x12, [%1]         \n\t"
        "ldp x13, x14, [%1, #16]    \n\t"
        "ldp x25, x26, [%1, #32]    \n\t"

        "adds x7, x7, x11           \n\t"
        "adcs x8, x8, x12           \n\t"
        "adcs x9, x9, x13           \n\t"
        "adcs x10, x10, x14         \n\t"
        "adcs x23, x23, x25         \n\t"
        "adcs x24, x24, x26         \n\t"
        "eor x16, x16, x16          \n\t"
        "adc x16, x16, x16          \n\t"

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

        "stp x11, x25, [%3]         \n\t"   // c0, c1
        
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

        "stp x26, x27, [%3, #16]     \n\t"  // c2, c3

        "mul x13, x3, x23           \n\t"
        "umulh x14, x3, x23         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"
        
        "mul x13, x7, x21           \n\t"
        "umulh x14, x7, x21         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x4, x10           \n\t"
        "umulh x14, x4, x10         \n\t"
        
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        
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

        "mul x13, x3, x24           \n\t"
        "umulh x14, x3, x24         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "mul x13, x7, x22           \n\t"
        "umulh x14, x7, x22         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x4, x23           \n\t"
        "umulh x14, x4, x23         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x8, x21           \n\t"
        "umulh x14, x8, x21         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x5, x10           \n\t"
        "umulh x14, x5, x10         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"
        
        "mul x13, x6, x9            \n\t"
        "umulh x14, x6, x9          \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"
        
        "stp x25, x18, [%3, #32]     \n\t"  // c4, c5

        "mul x13, x4, x24           \n\t"
        "umulh x14, x4, x24         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "mul x13, x5, x23           \n\t"
        "umulh x14, x5, x23         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x6, x10           \n\t"
        "umulh x14, x6, x10         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x21, x9           \n\t"
        "umulh x14, x21, x9         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x22, x8           \n\t"
        "umulh x14, x22, x8         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x3, x16           \n\t"
        "umulh x14, x3, x16          \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x7, x15           \n\t"
        "umulh x14, x7, x15         \n\t"
        
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"
        //-----------------------------------
        "mul x13, x5, x24           \n\t"
        "umulh x14, x5, x24         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x20, x20         \n\t"

        "mul x13, x6, x23           \n\t"
        "umulh x14, x6, x23         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x21, x10           \n\t"
        "umulh x14, x21, x10         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x22, x9           \n\t"
        "umulh x14, x22, x9         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x4, x16           \n\t"
        "umulh x14, x4, x16         \n\t"
        
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x8, x15           \n\t"
        "umulh x14, x8, x15          \n\t"
        
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "stp x26, x27, [%3, #48]    \n\t"   // c6, c7

        "mul x13, x6, x24           \n\t"
        "umulh x14, x6, x24         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"

        "mul x13, x21, x23           \n\t"
        "umulh x14, x21, x23         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x22, x10           \n\t"
        "umulh x14, x22, x10         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x5, x16           \n\t"
        "umulh x14, x5, x16          \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x9, x15           \n\t"
        "umulh x14, x9, x15          \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        //----------------------------------

        "mul x13, x21, x24           \n\t"
        "umulh x14, x21, x24         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "mul x13, x22, x23           \n\t"
        "umulh x14, x22, x23         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x6, x16           \n\t"
        "umulh x14, x6, x16         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x10, x15           \n\t"
        "umulh x14, x10, x15         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "stp x25, x18, [%3, #64]    \n\t"   // c8, c9

        "mul x13, x22, x24           \n\t"
        "umulh x14, x22, x24         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "mul x13, x21, x16          \n\t"
        "umulh x14, x21, x16        \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14          \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x23, x15          \n\t"
        "umulh x14, x23, x15        \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14          \n\t"
        "adcs x25, x25, x20         \n\t"
        //--------------------------------------
        "mul x13, x22, x16          \n\t"
        "umulh x14, x22, x16        \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14          \n\t"
        "adcs x18, x20, x20         \n\t"

        "mul x13, x24, x15          \n\t"
        "umulh x14, x24, x15        \n\t"
        
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14          \n\t"
        "adcs x18, x18, x20         \n\t"

        "stp x26, x27, [%3, #80]    \n\t"   // c10, c11

        "mul x13, x16, x15          \n\t"
        "umulh x14, x16, x15        \n\t"
        
        "adds x25, x25, x13         \n\t"
        "adc x18, x18, x14          \n\t"

        "stp x25, x18, [%3, #96]    \n\t"   // c12, c13




        :
        :"r"(&a[0]),"r"(&b[0]),"r"(&c[0]), "r"(&aplusbplus[0])
        :"memory", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16",
        "x18", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27"
    
    );    
}

void fpmul448_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{
    // Multiplication of two 384-bit integers using aarch64 assembly
    // Operation:
    // c[0...13] = a[0...6] * b[0...6] while a[6] and b[6] are either 0 or 1
   asm volatile(
        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x21, x22, [%0, #32]    \n\t"

        "ldp x7, x8, [%1]           \n\t"
        "ldp x9, x10, [%1, #16]     \n\t"
        "ldp x23, x24, [%1, #32]    \n\t"
        
        "eor x20, x20, x20          \n\t"   // use as zero register
        
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

        "stp x11, x25, [%2]         \n\t"   // c0, c1
        
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

        "stp x26, x27, [%2, #16]     \n\t"  // c2, c3

        "mul x13, x3, x23           \n\t"
        "umulh x14, x3, x23         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"
        
        "mul x13, x7, x21           \n\t"
        "umulh x14, x7, x21         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x4, x10           \n\t"
        "umulh x14, x4, x10         \n\t"
        
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        
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

        "mul x13, x3, x24           \n\t"
        "umulh x14, x3, x24         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "mul x13, x7, x22           \n\t"
        "umulh x14, x7, x22         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x4, x23           \n\t"
        "umulh x14, x4, x23         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x8, x21           \n\t"
        "umulh x14, x8, x21         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x5, x10           \n\t"
        "umulh x14, x5, x10         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"
        
        "mul x13, x6, x9            \n\t"
        "umulh x14, x6, x9          \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"
        
        "stp x25, x18, [%2, #32]     \n\t"  // c4, c5

        "mul x13, x4, x24           \n\t"
        "umulh x14, x4, x24         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "mul x13, x5, x23           \n\t"
        "umulh x14, x5, x23         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x6, x10           \n\t"
        "umulh x14, x6, x10         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x21, x9           \n\t"
        "umulh x14, x21, x9         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x22, x8           \n\t"
        "umulh x14, x22, x8         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "ldr x15, [%0, #48]         \n\t"   // a[6]
        "ldr x16, [%1, #48]         \n\t"   // b[6]

        "mul x13, x3, x16           \n\t"
        "umulh x14, x3, x16          \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x7, x15           \n\t"
        "umulh x14, x7, x15         \n\t"
        
        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"
        //-----------------------------------
        "mul x13, x5, x24           \n\t"
        "umulh x14, x5, x24         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x20, x20         \n\t"

        "mul x13, x6, x23           \n\t"
        "umulh x14, x6, x23         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x21, x10           \n\t"
        "umulh x14, x21, x10         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x22, x9           \n\t"
        "umulh x14, x22, x9         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x4, x16           \n\t"
        "umulh x14, x4, x16         \n\t"
        
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x8, x15           \n\t"
        "umulh x14, x8, x15          \n\t"
        
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "stp x26, x27, [%2, #48]    \n\t"   // c6, c7

        "mul x13, x6, x24           \n\t"
        "umulh x14, x6, x24         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"

        "mul x13, x21, x23           \n\t"
        "umulh x14, x21, x23         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x22, x10           \n\t"
        "umulh x14, x22, x10         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x5, x16           \n\t"
        "umulh x14, x5, x16          \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x9, x15           \n\t"
        "umulh x14, x9, x15          \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        //----------------------------------

        "mul x13, x21, x24           \n\t"
        "umulh x14, x21, x24         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "mul x13, x22, x23           \n\t"
        "umulh x14, x22, x23         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x6, x16           \n\t"
        "umulh x14, x6, x16         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x10, x15           \n\t"
        "umulh x14, x10, x15         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "stp x25, x18, [%2, #64]    \n\t"   // c8, c9

        "mul x13, x22, x24           \n\t"
        "umulh x14, x22, x24         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "mul x13, x21, x16          \n\t"
        "umulh x14, x21, x16        \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14          \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x23, x15          \n\t"
        "umulh x14, x23, x15        \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14          \n\t"
        "adcs x25, x25, x20         \n\t"
        //--------------------------------------
        "mul x13, x22, x16          \n\t"
        "umulh x14, x22, x16        \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14          \n\t"
        "adcs x18, x20, x20         \n\t"

        "mul x13, x24, x15          \n\t"
        "umulh x14, x24, x15        \n\t"
        
        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14          \n\t"
        "adcs x18, x18, x20         \n\t"

        "stp x26, x27, [%2, #80]    \n\t"   // c10, c11

        "mul x13, x16, x15          \n\t"
        "umulh x14, x16, x15        \n\t"
        
        "adds x25, x25, x13         \n\t"
        "adc x18, x18, x14          \n\t"

        "stp x25, x18, [%2, #96]    \n\t"   // c12, c13

        :
        :"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16",
        "x18", "x20", "x21", "x22", "x25", "x26", "x27"

    );
}


void fpadd384_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{    //Add two 384-bit values when there is no carry overflow
     // Operation:
     // c[0...6] = a[0...6] + b[0...6]
   asm volatile(

        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x13, x14, [%0, #32]    \n\t"

        "ldp x7, x8, [%1]           \n\t"
        "ldp x9, x10, [%1, #16]     \n\t"
        "ldp x15, x16, [%1, #32]    \n\t"

        "adds x3, x3, x7            \n\t"
        "adcs x4, x4, x8            \n\t"
        "adcs x5, x5, x9            \n\t"
        "adcs x6, x6, x10           \n\t"
        "adcs x13, x13, x15         \n\t"
        "adcs x14, x14, x16         \n\t"
        "eor x7, x7, x7             \n\t"
        "adc x7, x7, x7             \n\t"
        
        "stp x3, x4, [%2]           \n\t"
        "stp x5, x6, [%2, #16]      \n\t"
        "stp x13, x14, [%2, #32]    \n\t"
        "str x7, [%2, #48]          \n\t"

        :
        :"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12",
        "x13", "x14", "x15", "x16"
    );
}

void fpsub768_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{   // Sub two 768-bit values
   asm volatile(

        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x7, x8, [%0, #32]      \n\t"
        "ldp x9, x10, [%0, #48]     \n\t"
        "ldp x19, x20, [%0, #64]    \n\t"
        "ldp x21, x22, [%0, #80]    \n\t"

        "ldp x11, x12, [%1]         \n\t"
        "ldp x13, x14, [%1, #16]    \n\t"
        "ldp x15, x16, [%1, #32]    \n\t"
        "ldp x17, x18, [%1, #48]    \n\t"
        "ldp x23, x24, [%1, #64]    \n\t"
        "ldp x25, x26, [%1, #80]    \n\t"
    
        "subs x3, x3, x11           \n\t"
        "sbcs x4, x4, x12           \n\t"
        "sbcs x5, x5, x13           \n\t"
        "sbcs x6, x6, x14           \n\t"
        "sbcs x7, x7, x15           \n\t"
        "sbcs x8, x8, x16           \n\t"
        "sbcs x9, x9, x17           \n\t"
        "sbcs x10, x10, x18         \n\t"
        "sbcs x19, x19, x23         \n\t"
        "sbcs x20, x20, x24         \n\t"
        "sbcs x21, x21, x25         \n\t"
        "sbcs x22, x22, x26         \n\t"
        
        "ldr x11, [%0, #96]         \n\t"
        "eor x12, x12, x12          \n\t"
        "sbc x11, x11, x12          \n\t"
    
        "stp x3, x4, [%2]           \n\t"
        "stp x5, x6, [%2, #16]      \n\t"
        "stp x7, x8, [%2, #32]      \n\t"
        "stp x9, x10, [%2, #48]     \n\t"
        "stp x19, x20, [%2, #64]    \n\t"
        "stp x21, x22, [%2, #80]    \n\t"
        "str x11, [%2, #96]         \n\t"

        :
        :"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26"
    );
}

void fpadd832_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{   // add two 768-bit values and propagate the carry overflow
    // operation 
    // c[0...12] = a[0...12] + b[0...12]
   asm volatile(

        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x7, x8, [%0, #32]      \n\t"
        "ldp x9, x10, [%0, #48]     \n\t"
        "ldp x19, x20, [%0, #64]    \n\t"
        "ldp x21, x22, [%0, #80]    \n\t"


        "ldp x11, x12, [%1]         \n\t"
        "ldp x13, x14, [%1, #16]    \n\t"
        "ldp x15, x16, [%1, #32]    \n\t"
        "ldp x17, x18, [%1, #48]    \n\t"
        "ldp x23, x24, [%1, #64]    \n\t"
        "ldp x25, x26, [%1, #80]    \n\t"
    
        "adds x3, x3, x11           \n\t"
        "adcs x4, x4, x12           \n\t"
        "adcs x5, x5, x13           \n\t"
        "adcs x6, x6, x14           \n\t"
        "adcs x7, x7, x15           \n\t"
        "adcs x8, x8, x16           \n\t"
        "adcs x9, x9, x17           \n\t"
        "adcs x10, x10, x18         \n\t"
        "adcs x19, x19, x23         \n\t"
        "adcs x20, x20, x24         \n\t"
        "adcs x21, x21, x25         \n\t"
        "adcs x22, x22, x26         \n\t"

        "ldr x27, [%0, #96]         \n\t"
        "ldr x11, [%1, #96]         \n\t"

        "adc  x27, x27, x11         \n\t"
    
        "stp x3, x4, [%2]           \n\t"
        "stp x5, x6, [%2, #16]      \n\t"
        "stp x7, x8, [%2, #32]      \n\t"
        "stp x9, x10, [%2, #48]     \n\t"
        "stp x19, x20, [%2, #64]    \n\t"
        "stp x21, x22, [%2, #80]    \n\t"
        "str x27, [%2, #96]         \n\t"

        :
        :"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27"
    );
}

void fpmul384_arm64(uint64_t *a, uint64_t *b, uint64_t *c)
{
    // Multiplication of two 384-bit integers using aarch64 assembly
    // Operation:
    // c[0...11] = a[0...5] * b[0...5]
   asm volatile(
        "ldp x3, x4, [%0]           \n\t"
        "ldp x5, x6, [%0, #16]      \n\t"
        "ldp x21, x22, [%0, #32]    \n\t"

        "ldp x7, x8, [%1]           \n\t"
        "ldp x9, x10, [%1, #16]     \n\t"
        "ldp x23, x24, [%1, #32]    \n\t"
        
        "eor x20, x20, x20          \n\t"   // use as zero register
        
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

        "stp x11, x25, [%2]         \n\t"   // c0, c1
        
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

        "stp x26, x27, [%2, #16]     \n\t"  // c2, c3

        "mul x13, x3, x23           \n\t"
        "umulh x14, x3, x23         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"
        
        "mul x13, x7, x21           \n\t"
        "umulh x14, x7, x21         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x4, x10           \n\t"
        "umulh x14, x4, x10         \n\t"
        
        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"
        
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
        //-----------------------------------
        "mul x13, x3, x24           \n\t"
        "umulh x14, x3, x24         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "mul x13, x7, x22           \n\t"
        "umulh x14, x7, x22         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x4, x23           \n\t"
        "umulh x14, x4, x23         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x8, x21           \n\t"
        "umulh x14, x8, x21         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "mul x13, x5, x10           \n\t"
        "umulh x14, x5, x10         \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"
        
        "mul x13, x6, x9            \n\t"
        "umulh x14, x6, x9          \n\t"
        
        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"
        
        "stp x25, x18, [%2, #32]     \n\t"  // c4, c5

        "mul x13, x4, x24           \n\t"
        "umulh x14, x4, x24         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x20, x20         \n\t"

        "mul x13, x5, x23           \n\t"
        "umulh x14, x5, x23         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x6, x10           \n\t"
        "umulh x14, x6, x10         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x21, x9           \n\t"
        "umulh x14, x21, x9         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x22, x8           \n\t"
        "umulh x14, x22, x8         \n\t"

        "adds x26, x26, x13         \n\t"
        "adcs x27, x27, x14         \n\t"
        "adcs x25, x25, x20         \n\t"

        "mul x13, x5, x24           \n\t"
        "umulh x14, x5, x24         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x20, x20         \n\t"

        "mul x13, x6, x23           \n\t"
        "umulh x14, x6, x23         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x21, x10           \n\t"
        "umulh x14, x21, x10         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "mul x13, x22, x9           \n\t"
        "umulh x14, x22, x9         \n\t"

        "adds x27, x27, x13         \n\t"
        "adcs x25, x25, x14         \n\t"
        "adcs x18, x18, x20         \n\t"

        "stp x26, x27, [%2, #48]    \n\t"   // c6, c7

        "mul x13, x6, x24           \n\t"
        "umulh x14, x6, x24         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x20, x20         \n\t"

        "mul x13, x21, x23           \n\t"
        "umulh x14, x21, x23         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x22, x10           \n\t"
        "umulh x14, x22, x10         \n\t"

        "adds x25, x25, x13         \n\t"
        "adcs x18, x18, x14         \n\t"
        "adcs x26, x26, x20         \n\t"

        "mul x13, x21, x24           \n\t"
        "umulh x14, x21, x24         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x20, x20         \n\t"

        "mul x13, x22, x23           \n\t"
        "umulh x14, x22, x23         \n\t"

        "adds x18, x18, x13         \n\t"
        "adcs x26, x26, x14         \n\t"
        "adcs x27, x27, x20         \n\t"

        "stp x25, x18, [%2, #64]    \n\t"   // c8, c9

        "mul x13, x22, x24           \n\t"
        "umulh x14, x22, x24         \n\t"

        "adds x26, x26, x13         \n\t"
        "adc x27, x27, x14          \n\t"

        "stp x26, x27, [%2, #80]    \n\t"   // c8, c9

        :
        :"r"(&a[0]),"r"(&b[0]),"r"(&c[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x18", "x20", "x21", "x22", "x25", "x26", "x27"

    );
}




void fpmul768_karatsuba(uint64_t *a, uint64_t *b, uint64_t *c)
{
    uint64_t rplus[14];

    fpmul2x384_mixed_arm64(a, b, c, rplus);
    fpsub768_arm64(rplus, c, rplus);
    fpsub768_arm64(rplus, c+12, rplus);
    fpadd832_arm64(c+6, rplus, c+6);
}

void fpsqr768_asm(uint64_t *res, uint64_t *a){
    
   asm volatile(	
        //step 1
        "ldr x0, [%1]		\n\t"
        "ldp x1, x2, [%1, #8]	\n\t"
    
        "mul x3, x0, x0		\n\t"
        "umulh x4, x0, x0	\n\t"
        "str x3, [%0]		\n\t"	//c0
    
        "mul x5, x1, x0		\n\t"
        "umulh x6, x1, x0	\n\t"
        "mul x7, x2, x0		\n\t"
        "umulh x8, x2, x0	\n\t"	
        "adds x4, x4, x5	\n\t"
        "adcs x6, x6, x6	\n\t"
        "adcs x8, x8, x8	\n\t"
    
        "ldp x1, x2, [%1, #8 + 16]	\n\t"
        "mul x9, x1, x0		\n\t"
        "umulh x10, x1, x0	\n\t"
        "mul x11, x2, x0	\n\t"
        "umulh x12, x2, x0	\n\t"
        "adcs x10, x10, x10	\n\t"
        "adcs x12, x12, x12	\n\t"
    
        "ldp x1, x2, [%1, #8 + 32]	\n\t"
        "mul x13, x1, x0	\n\t"
        "umulh x14, x1, x0	\n\t"
        "mul x15, x2, x0	\n\t"
        "umulh x16, x2, x0	\n\t"
        "adcs x14, x14, x14	\n\t"
        "adcs x16, x16, x16	\n\t"
        
        "ldp x1, x2, [%1, #8 + 48]	\n\t"
        "mul x17, x1, x0	\n\t"
        "umulh x18, x1, x0	\n\t"
        "mul x19, x2, x0	\n\t"
        "umulh x20, x2, x0	\n\t"
        "adcs x18, x18, x18	\n\t"
        "adcs x20, x20, x20	\n\t"
      
        "ldp x1, x2, [%1, #8 + 64]	\n\t"
        "mul x21, x1, x0	\n\t"
        "umulh x22, x1, x0	\n\t"
        "mul x23, x2, x0	\n\t"
        "umulh x24, x2, x0	\n\t"
        "adcs x22, x22, x22	\n\t"
        "adcs x24, x24, x24	\n\t"
    
        "ldr x1, [%1, #8 + 80]		\n\t"
        "mul x25, x1, x0	\n\t"
        "umulh x26, x1, x0	\n\t"
        "adcs x26, x26, x26	\n\t"
        "eor x1, x1, x1		\n\t"
        "eor x3, x3, x3		\n\t"
        "adc x3, x3, x1		\n\t"	//carry 
        //addition and carry propagation
        "adds x4, x4, x5	\n\t"
        "str x4, [%0, #8 * 1]	\n\t"//c1
        "adcs x7, x7, x7	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x11, x11, x11	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x15, x15, x15	\n\t"
        "adcs x17, x17, x17	\n\t"
        "adcs x19, x19, x19	\n\t"
        "adcs x21, x21, x21	\n\t"
        "adcs x23, x23, x23	\n\t"
        "adcs x25, x25, x25	\n\t"
        "adcs x26, x26, x1	\n\t"
        "adc x3, x3, x1		\n\t"
        
        "adds x6, x6, x7	\n\t"
        "adcs x8, x8, x9	\n\t"
        "adcs x10, x10, x11	\n\t"
        "adcs x12, x12, x13	\n\t"
        "adcs x14, x14, x15	\n\t"
        "adcs x16, x16, x17	\n\t"
        "adcs x18, x18, x19	\n\t"
        "adcs x20, x20, x21	\n\t"
        "adcs x22, x22, x23	\n\t"
        "adcs x24, x24, x25	\n\t"
        "adcs x26, x26, x1	\n\t"
        "adc x3, x3, x1		\n\t"
        //step 2
        
        "ldr x0, [%1, #8]		\n\t"
        "ldp x1, x2, [%1, #16]	\n\t"
    
        "mul x4, x0, x0		\n\t"
        "umulh x5, x0, x0	\n\t"
        "adds x6, x6, x4	\n\t"
        "str x6, [%0, #8 * 2]	\n\t"	//c2
    
        "mul x7, x1, x0		\n\t"
        "umulh x9, x1, x0	\n\t"
        "mul x11, x2, x0	\n\t"
        "umulh x13, x2, x0	\n\t"	
        "adcs x5, x5, x7	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x13, x13, x13	\n\t"
    
        "ldp x1, x2, [%1, #16 + 16]	\n\t"
        "mul x15, x1, x0		\n\t"
        "umulh x17, x1, x0	\n\t"
        "mul x19, x2, x0	\n\t"
        "umulh x21, x2, x0	\n\t"
        "adcs x17, x17, x17	\n\t"
        "adcs x21, x21, x21	\n\t"
        
        "ldp x1, x2, [%1, #16 + 32]	\n\t"
        "mul x6, x1, x0		\n\t"
        "umulh x4, x1, x0	\n\t"
        "adcs x4, x4, x4	\n\t"
        "eor x23, x23, x23	\n\t"
        "adc x23, x23, x23	\n\t"//save middle carry
        
        "adds x5, x5, x7	\n\t"
        "adcs x11, x11, x11	\n\t"
        "adcs x15, x15, x15	\n\t"
        "adcs x19, x19, x19	\n\t"
        "adcs x6, x6, x6	\n\t"
        "eor x25, x25, x25	\n\t"
        "adc x25, x25, x25	\n\t"//save second carry
    
        "adds x5, x5, x8	\n\t"
        "str x5, [%0, #8 * 3]	\n\t"	//c3
        "adcs x9, x9, x11	\n\t"
        "adcs x13, x13, x15	\n\t"
        "adcs x17, x17, x19	\n\t"
        "adcs x21, x21, x6	\n\t"
        "eor x7, x7, x7		\n\t"
        "adc x25, x25, x7	\n\t"
    
        "adds x10, x10, x9	\n\t"
        "adcs x12, x12, x13	\n\t"
        "adcs x14, x14, x17	\n\t"
        "adcs x16, x16, x21	\n\t"
        "adc x25, x25, x7	\n\t"
    
        "mul x7, x2, x0		\n\t"
        "umulh x9, x2, x0	\n\t"
        "adds x25, x25, x4	\n\t"
        "adcs x9, x9, x9	\n\t"
    
        "ldp x1, x2, [%1, #16 + 48]	\n\t"
        "mul x11, x1, x0	\n\t"
        "umulh x13, x1, x0	\n\t"
        "mul x15, x2, x0	\n\t"
        "umulh x17, x2, x0	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x17, x17, x17	\n\t"
        
        "ldp x1, x2, [%1, #16 + 64]	\n\t"
        "mul x19, x1, x0	\n\t"
        "umulh x21, x1, x0	\n\t"
        "mul x5, x2, x0		\n\t"
        "umulh x6, x2, x0	\n\t"
        "adcs x21, x21, x21	\n\t"
        "adcs x6, x6, x6	\n\t"
        "eor x4, x4, x4		\n\t"
        "adc x4, x4, x4		\n\t"//final carry 
    
        "adds x7, x7, x7	\n\t"
        "adcs x11, x11, x11	\n\t"
        "adcs x15, x15, x15	\n\t"
        "adcs x19, x19, x19	\n\t"
        "adcs x5, x5, x5	\n\t"
        "eor x0, x0, x0		\n\t"
        "adcs x6, x6, x0	\n\t"
        "adc x4, x4, x0		\n\t"
        
        "adds x25, x25, x7	\n\t"
        "adcs x9, x9, x11	\n\t"
        "adcs x13, x13, x15	\n\t"
        "adcs x17, x17, x19	\n\t"
        "adcs x21, x21, x5	\n\t"
        "eor x0, x0, x0		\n\t"
        "adcs x6, x6, x0	\n\t"
        "adc x4, x4, x0		\n\t"
        
        "adds x18, x18, x25	\n\t"
        "adcs x23, x23, x9	\n\t"
        "adcs x22, x22, x13	\n\t"
        "adcs x24, x24, x17	\n\t"
        "adcs x26, x26, x21	\n\t"
        "adcs x3, x3, x6	\n\t"
        "adc x4, x4, x0		\n\t"
        
        "eor x0, x0, x0		\n\t"
        "adds x20, x20, x23	\n\t"
        "adcs x22, x22, x0	\n\t"
        "adcs x24, x24, x0	\n\t"
        "adcs x26, x26, x0	\n\t"
        "adcs x3, x3, x0	\n\t"
        "adc x4, x4, x0		\n\t"
        //step 3 
        "ldr x0, [%1, #16]	\n\t"
        "ldp x1, x2, [%1, #24]	\n\t"
        "mul x5, x0, x0		\n\t"
        "umulh x6, x0, x0	\n\t"
        
        "adds x10, x10, x5	\n\t"
        "str x10, [%0, #8 * 4]	\n\t"	//c4 	
        "mul x7, x1, x0		\n\t"
        "umulh x8, x1, x0	\n\t"
        "mul x9, x2, x0		\n\t"
        "umulh x11, x2, x0	\n\t"
        "adcs x6, x6, x7	\n\t"
        "adcs x8, x8, x8	\n\t"
        "adcs x11, x11, x11	\n\t"
    
        "ldp x1, x2, [%1, #24 + 16]	\n\t"
        "mul x13, x1, x0	\n\t"
        "umulh x15, x1, x0	\n\t"
        "mul x17, x2, x0	\n\t"
        "umulh x19, x2, x0	\n\t"
        "adcs x15, x15, x15	\n\t"
        "adcs x19, x19, x19	\n\t"
        
        "ldp x1, x2, [%1, #24 + 32]	\n\t"
        "mul x21, x1, x0	\n\t"
        "umulh x23, x1, x0	\n\t"
        "adcs x23, x23, x23	\n\t"
        "eor x25, x25, x25	\n\t"
        "adc x25, x25, x25	\n\t"//first carry
    
        "adds x6, x6, x7	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x17, x17, x17	\n\t"
        "adcs x21, x21, x21	\n\t"
        "eor x10, x10, x10	\n\t"
        "adc x10, x10, x10	\n\t"//second carry
    
        "adds x12, x12, x6	\n\t"
        "str x12, [%0, #8 * 5]	\n\t"	//c5
        "adcs x8, x8, x9	\n\t"
        "adcs x11, x11, x13	\n\t"
        "adcs x15, x15, x17	\n\t"
        "adcs x19, x19, x21	\n\t"
        "adcs x10, x10, x23	\n\t"
        "eor x6, x6, x6		\n\t"
        "adc x25, x25, x6	\n\t"
        
        "adds x14, x14, x8	\n\t"
        "adcs x16, x16, x11	\n\t"
        "adcs x18, x18, x15	\n\t"
        "adcs x20, x20, x19	\n\t"
        "adcs x10, x10, x6	\n\t"
        "adc x25, x25, x6	\n\t"
    
        "mul x5, x2, x0		\n\t"
        "umulh x6, x2, x0	\n\t"
        "adds x10, x10, x5	\n\t"
        "adcs x6, x6, x6	\n\t"
        
        "ldp x1, x2, [%1, #24 + 48]	\n\t"
        "mul x7, x1, x0		\n\t"
        "umulh x8, x1, x0	\n\t"
        "mul x9, x2, x0		\n\t"
        "umulh x11, x2, x0	\n\t"
        "adcs x8, x8, x8	\n\t"
        "adcs x11, x11, x11	\n\t"
    
        "ldr x1, [%1, #24 + 64]	\n\t"
        "mul x13, x1, x0	\n\t"
        "umulh x15, x1, x0	\n\t"
        "adcs x15, x15, x15	\n\t"
        "eor x17, x17, x17	\n\t"
        "adc x17, x17, x17	\n\t"
    
        "adds x10, x10, x5	\n\t"
        "adcs x7, x7, x7	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x13, x13, x13	\n\t"
        "eor x5, x5, x5		\n\t"
        "adcs x15, x15, x5	\n\t"
        "adc x17, x17, x5	\n\t"
        
        "adds x22, x22, x10	\n\t"
        "adcs x6, x6, x7	\n\t"
        "adcs x8, x8, x9	\n\t"
        "adcs x11, x11, x13	\n\t"
        "adcs x15, x15, x5	\n\t"
        "adc x17, x17, x5	\n\t"
    
        "eor x5, x5, x5		\n\t"
        "adds x25, x25, x6	\n\t"
        "adcs x26, x26, x8	\n\t"
        "adcs x3, x3, x11	\n\t"
        "adcs x4, x4, x15	\n\t"
        "adc x17, x17, x5	\n\t"
    
        "eor x5, x5, x5		\n\t"	
        "adds x24, x24, x25	\n\t"
        "adcs x26, x26, x5	\n\t"
        "adcs x3, x3, x5	\n\t"
        "adcs x4, x4, x5	\n\t"
        "adc x17, x17, x5	\n\t"
        //step 4 
        "ldr x0, [%1, #24]	\n\t"
        "ldp x1, x2, [%1, #32]	\n\t"
        
        "mul x5, x0, x0		\n\t"
        "umulh x6, x0, x0	\n\t"
        "adds x14, x14, x5	\n\t"
        "str x14, [%0, #8 * 6]	\n\t"	//c6
        
        "mul x7, x1, x0		\n\t"
        "umulh x8, x1, x0	\n\t"
        "mul x9, x2, x0		\n\t"
        "umulh x10, x2, x0	\n\t"
        "adcs x6, x6, x7	\n\t"
        "adcs x8, x8, x8	\n\t"
        "adcs x10, x10, x10	\n\t"
    
        "ldp x1, x2, [%1, #32 + 16]	\n\t"
        "mul x11, x1, x0	\n\t"
        "umulh x12, x1, x0	\n\t"
        "mul x13, x2, x0	\n\t"
        "umulh x15, x2, x0	\n\t"
        "adcs x12, x12, x12	\n\t"
        "adcs x15, x15, x15	\n\t"
        
        "ldp x1, x2, [%1, #32 + 32]	\n\t"
        "mul x19, x1, x0	\n\t"
        "umulh x21, x1, x0	\n\t"
        "mul x23, x2, x0	\n\t"
        "umulh x25, x2, x0	\n\t"
        "adcs x21, x21, x21	\n\t"
        "adcs x25, x25, x25	\n\t"
        
        "ldp x1, x2, [%1, #32 + 48]	\n\t"
        "mul x14, x1, x0	\n\t"
        "umulh x5, x1, x0	\n\t"
        "mul x1, x2, x0		\n\t"
        "umulh x2, x2, x0	\n\t"
        "adcs x5, x5, x5	\n\t"
        "adcs x2, x2, x2	\n\t"
        "eor x0, x0, x0		\n\t"
        "adc x0, x0, x0		\n\t"
            
        "adds x6, x6, x7	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x11, x11, x11	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x19, x19, x19	\n\t"
        "adcs x23, x23, x23	\n\t"
        "adcs x14, x14, x14	\n\t"
        "adcs x1, x1, x1	\n\t"
        "eor x7, x7, x7		\n\t"
        "adcs x2, x2, x7	\n\t"
        "adc x0, x0, x7		\n\t"
    
        "adds x16, x16, x6	\n\t"
        "str x16, [%0, #8 * 7]	\n\t"	//c7
        "adcs x8, x8, x9	\n\t"
        "adcs x10, x10, x11	\n\t"
        "adcs x12, x12, x13	\n\t"
        "adcs x15, x15, x19	\n\t"
        "adcs x21, x21, x23	\n\t"
        "adcs x25, x25, x14	\n\t"
        "adcs x5, x5, x1	\n\t"
        "eor x7, x7, x7		\n\t"
        "adcs x2, x2, x7	\n\t"
        "adc x0, x0, x7		\n\t"
    
        "adds x18, x18, x8	\n\t"
        "adcs x20, x20, x10	\n\t"
        "adcs x22, x22, x12	\n\t"
        "adcs x24, x24, x15	\n\t"
        "adcs x26, x26, x21	\n\t"
        "adcs x3, x3, x25	\n\t"
        "adcs x4, x4, x5	\n\t"
        "adcs x17, x17, x2	\n\t"
        "eor x7, x7, x7		\n\t"
        "adc x7, x7, x0		\n\t"
        //step 5	
        "ldr x0, [%1, #32]	\n\t"
        "ldp x1, x2, [%1, #40]	\n\t"
            
        "mul x5, x0, x0		\n\t"
        "umulh x6, x0, x0	\n\t"
        "adds x18, x18, x5	\n\t"
        "str x18, [%0, #8 * 8]	\n\t"	//c8
        "mul x8, x1, x0		\n\t"
        "umulh x9, x1, x0	\n\t"
        "mul x10, x2, x0	\n\t"
        "umulh x11, x2, x0	\n\t"
        "adcs x6, x6, x8	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x11, x11, x11	\n\t"
        
        "ldp x1, x2, [%1, #40 + 16]	\n\t"
        "mul x12, x1, x0	\n\t"
        "umulh x13, x1, x0	\n\t"
        "mul x14, x2, x0	\n\t"
        "umulh x15, x2, x0	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x15, x15, x15	\n\t"
        
        "ldp x1, x2, [%1, #40 + 32]	\n\t"
        "mul x16, x1, x0	\n\t"
        "umulh x18, x1, x0	\n\t"
        "mul x19, x2, x0	\n\t"
        "umulh x21, x2, x0	\n\t"
        "adcs x18, x18, x18	\n\t"
        "adcs x21, x21, x21	\n\t"
        
        "ldr x1, [%1, #40 + 48]		\n\t"
        "mul x23, x1, x0	\n\t"
        "umulh x25, x1, x0	\n\t"
        "adcs x25, x25, x25	\n\t"
        "eor x5, x5, x5		\n\t"
        "adc x5, x5, x5		\n\t"
        
        "adds x6, x6, x8	\n\t"
        "adcs x10, x10, x10	\n\t"
        "adcs x12, x12, x12	\n\t"
        "adcs x14, x14, x14	\n\t"
        "adcs x16, x16, x16	\n\t"
        "adcs x19, x19, x19	\n\t"
        "adcs x23, x23, x23	\n\t"
        "eor x8, x8, x8		\n\t"
        "adcs x25, x25, x8	\n\t"
        "adc x5, x5, x8		\n\t"
    
        
        "adds x20, x20, x6	\n\t"
        "str x20, [%0, #8 * 9]	\n\t"	//c9
        "adcs x9, x9, x10	\n\t"
        "adcs x11, x11, x12	\n\t"
        "adcs x13, x13, x14	\n\t"
        "adcs x15, x15, x16	\n\t"
        "adcs x18, x18, x19	\n\t"
        "adcs x21, x21, x23	\n\t"
        "eor x6, x6, x6		\n\t"
        "adcs x7, x7, x25	\n\t"
        "adc x5, x5, x6		\n\t"
    
        "adcs x22, x22, x9	\n\t"
        "adcs x24, x24, x11	\n\t"
        "adcs x26, x26, x13	\n\t"
        "adcs x3, x3, x15	\n\t"
        "adcs x4, x4, x18	\n\t"
        "adcs x17, x17, x21	\n\t"
        "adcs x7, x7, x6	\n\t"
        "adc x5, x5, x6		\n\t"	
        //step 6
        "ldr x0, [%1, #40]	\n\t"
        "ldp x1, x2, [%1, #48]	\n\t"
        
        "mul x6, x0, x0		\n\t"
        "umulh x8, x0, x0	\n\t"
        "adds x22, x22, x6	\n\t"
        "str x22, [%0, #8 * 10]	\n\t"	//c10
        "mul x9, x1, x0		\n\t"
        "umulh x10, x1, x0	\n\t"
        "mul x11, x2, x0	\n\t"
        "umulh x12, x2, x0	\n\t"
        "adcs x8, x8, x9	\n\t"
        "adcs x10, x10, x10	\n\t"
        "adcs x12, x12, x12	\n\t"
        
        "ldp x1, x2, [%1, #48 + 16]	\n\t"
        "mul x13, x1, x0	\n\t"
        "umulh x14, x1, x0	\n\t"
        "mul x15, x2, x0	\n\t"
        "umulh x16, x2, x0	\n\t"
        "adcs x14, x14, x14	\n\t"
        "adcs x16, x16, x16	\n\t"
        
        "ldp x1, x2, [%1, #48 + 32]	\n\t"
        "mul x18, x1, x0	\n\t"
        "umulh x19, x1, x0	\n\t"
        "mul x20, x2, x0	\n\t"
        "umulh x21, x2, x0	\n\t"
        "adcs x19, x19, x19	\n\t"
        "adcs x21, x21, x21	\n\t"
        "eor x23, x23, x23	\n\t"
        "adc x23, x23, x23	\n\t"
    
        "adds x8, x8, x9	\n\t"
        "adcs x11, x11, x11	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x15, x15, x15	\n\t"
        "adcs x18, x18, x18	\n\t"
        "adcs x20, x20, x20	\n\t"
        "eor x6, x6, x6		\n\t"
        "adcs x21, x21, x6	\n\t"
        "adc x23, x23, x6	\n\t"
        
        "adds x24, x24, x8	\n\t"
        "str x24, [%0, #8 * 11]	\n\t"	//c11	
        "adcs x10, x10, x11	\n\t"
        "adcs x12, x12, x13	\n\t"
        "adcs x14, x14, x15	\n\t"
        "adcs x16, x16, x18	\n\t"
        "adcs x19, x19, x20	\n\t"
        "adcs x21, x21, x6	\n\t"
        "adc x23, x23, x6	\n\t"
        
        "adds x26, x26, x10	\n\t"
        "adcs x3, x3, x12	\n\t"
        "adcs x4, x4, x14	\n\t"
        "adcs x17, x17, x16	\n\t"
        "adcs x7, x7, x19	\n\t"
        "adcs x5, x5, x21	\n\t"
        "adc x23, x23, x6	\n\t"
        //step 7
        "ldr x0, [%1, #48]	\n\t"
        "ldp x1, x2, [%1, #56]	\n\t"
        "mul x6, x0, x0		\n\t"
        "umulh x8, x0, x0	\n\t"
        "adds x26, x26, x6	\n\t"
        "str x26, [%0, #8 * 12]	\n\t"	//c12
        "mul x9, x1, x0		\n\t"
        "umulh x10, x1, x0	\n\t"
        "mul x11, x2, x0	\n\t"
        "umulh x12, x2, x0	\n\t"
        "adcs x8, x8, x9	\n\t"
        "adcs x10, x10, x10	\n\t"
        "adcs x12, x12, x12	\n\t"
        
        "ldp x1, x2, [%1, #56 + 16]	\n\t"
        "mul x13, x1, x0	\n\t"
        "umulh x14, x1, x0	\n\t"
        "mul x15, x2, x0	\n\t"
        "umulh x16, x2, x0	\n\t"
        "adcs x14, x14, x14	\n\t"
        "adcs x16, x16, x16	\n\t"
        
        "ldr x1, [%1, #56 + 32]	\n\t"
        "mul x18, x1, x0	\n\t"
        "umulh x19, x1, x0	\n\t"
        "adcs x19, x19, x19	\n\t"
        "eor x20, x20, x20	\n\t"
        "adc x20, x20, x20	\n\t"
    
        "adds x8, x8, x9	\n\t"
        "adcs x11, x11, x11	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x15, x15, x15	\n\t"
        "adcs x18, x18, x18	\n\t"
        "eor x6, x6, x6		\n\t"
        "adcs x19, x19, x6	\n\t"
        "adc x20, x20, x6	\n\t"
    
        "adds x3, x3, x8	\n\t"
        "str x3, [%0, #8 * 13]	\n\t"	//c13
        "adcs x10, x10, x11	\n\t"
        "adcs x12, x12, x13	\n\t"
        "adcs x14, x14, x15	\n\t"
        "adcs x16, x16, x18	\n\t"
        "adcs x19, x19, x6	\n\t"
        "adc x20, x20, x6	\n\t"
        
        "adds x4, x4, x10	\n\t"
        "adcs x17, x17, x12	\n\t"
        "adcs x7, x7, x14	\n\t"
        "adcs x5, x5, x16	\n\t"
        "adcs x23, x23, x19	\n\t"
        "adc x20, x20, x6	\n\t"
        //step 8
        "ldr x0, [%1, #56]	\n\t"
        "ldp x1, x2, [%1, #64]	\n\t"
        "mul x3, x0, x0		\n\t"
        "umulh x6, x0, x0	\n\t"
        "adds x4, x4, x3	\n\t"
        "str x4, [%0, #8 * 14]	\n\t"	//c14
        "mul x8, x1, x0		\n\t"
        "umulh x9, x1, x0	\n\t"
        "mul x10, x2, x0	\n\t"
        "umulh x11, x2, x0	\n\t"
        "adcs x6, x6, x8	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x11, x11, x11	\n\t"
        
        "ldp x1, x2, [%1, #64 + 16]	\n\t"
        "mul x12, x1, x0	\n\t"
        "umulh x13, x1, x0	\n\t"
        "mul x14, x2, x0	\n\t"
        "umulh x15, x2, x0	\n\t"
        "adcs x13, x13, x13	\n\t"
        "adcs x15, x15, x15	\n\t"
        "eor x16, x16, x16	\n\t"
        "adc x16, x16, x16	\n\t"
    
        "adds x6, x6, x8	\n\t"
        "adcs x10, x10, x10	\n\t"
        "adcs x12, x12, x12	\n\t"
        "adcs x14, x14, x14	\n\t"
        "eor x3, x3, x3		\n\t"
        "adcs x15, x15, x3	\n\t"
        "adc x16, x16, x3	\n\t"
        
        "adds x17, x17, x6	\n\t"
        "str x17, [%0, #8 * 15]	\n\t"	//c15
        "adcs x9, x9, x10	\n\t"
        "adcs x11, x11, x12	\n\t"
        "adcs x13, x13, x14	\n\t"
        "adcs x15, x15, x3	\n\t"
        "adc x16, x16, x3	\n\t"
        
        "adds x7, x7, x9	\n\t"
        "adcs x5, x5, x11	\n\t"
        "adcs x23, x23, x13	\n\t"
        "adcs x20, x20, x15	\n\t"
        "adc x16, x16, x3	\n\t"
        //step 9
        "ldr x0, [%1, #64]	\n\t"
        "ldp x1, x2, [%1, #72]	\n\t"
        "mul x3, x0, x0		\n\t"
        "umulh x4, x0, x0	\n\t"
        "adds x7, x7, x3	\n\t"
        "str x7, [%0, #8 * 16]	\n\t"	//c16
        "mul x6, x1, x0		\n\t"
        "umulh x8, x1, x0	\n\t"
        "mul x9, x2, x0		\n\t"
        "umulh x10, x2, x0	\n\t"
        "adcs x4, x4, x6	\n\t"
        "adcs x8, x8, x8	\n\t"
        "adcs x10, x10, x10	\n\t"
        
        "ldr x1, [%1, #72 + 16]	\n\t"
        "mul x11, x1, x0	\n\t"
        "umulh x12, x1, x0	\n\t"
        "adcs x12, x12, x12	\n\t"
        "eor x13, x13, x13	\n\t"
        "adc x13, x13, x13	\n\t"
        
        "adds x4, x4, x6	\n\t"
        "adcs x9, x9, x9	\n\t"
        "adcs x11, x11, x11	\n\t"
        "eor x14, x14, x14	\n\t"
        "adcs x12, x12, x14	\n\t"
        "adc x13, x13, x14	\n\t"
        
        "adds x5, x5, x4	\n\t"
        "str x5, [%0, #8 * 17]	\n\t"	//c17
        "adcs x8, x8, x9	\n\t"
        "adcs x10, x10, x11	\n\t"
        "adcs x12, x12, x14	\n\t"
        "adc x13, x13, x14	\n\t"
        
        "adds x23, x23, x8	\n\t"
        "adcs x20, x20, x10	\n\t"
        "adcs x16, x16, x12	\n\t"
        "adc x13, x13, x14	\n\t"
        //step 10
        "ldr x0, [%1, #72]	\n\t"
        "ldp x1, x2, [%1, #80]	\n\t"
        "mul x3, x0, x0		\n\t"
        "umulh x4, x0, x0	\n\t"
        "adds x23, x23, x3	\n\t"
        "str x23, [%0, #8 * 18]	\n\t"	//c18 
        "mul x5, x1, x0		\n\t"
        "umulh x6, x1, x0	\n\t"
        "mul x7, x2, x0		\n\t"
        "umulh x8, x2, x0	\n\t"
        "adcs x4, x4, x5	\n\t"
        "adcs x6, x6, x6	\n\t"
        "adcs x8, x8, x8	\n\t"
        "eor x9, x9, x9		\n\t"
        "adc x9, x9, x9		\n\t"
        
        "adds x4, x4, x5	\n\t"
        "adcs x7, x7, x7	\n\t"
        "eor x10, x10, x10	\n\t"
        "adcs x8, x8, x10	\n\t"
        "adc x9, x9, x10	\n\t"
    
        "adds x20, x20, x4	\n\t"
        "str x20, [%0, #8 * 19]	\n\t"	//c19
        "adcs x6, x6, x7	\n\t"
        "adcs x8, x8, x10	\n\t"
        "adc x9, x9, x10	\n\t"
    
        "adds x16, x16, x6	\n\t"
        "adcs x13, x13, x8	\n\t"
        "adc x9, x9, x10	\n\t"
        //step 11
        "ldp x0, x1, [%1, #80]	\n\t"
        "mul x3, x0, x0		\n\t"
        "umulh x4, x0, x0	\n\t"
        "adds x16, x16, x3	\n\t"
        "str x16, [%0, #8 * 20]	\n\t"	//c20
        "mul x5, x1, x0		\n\t"
        "umulh x6, x1, x0	\n\t"
        "adcs x4, x4, x5	\n\t"
        "adcs x6, x6, x6	\n\t"
        "eor x7, x7, x7		\n\t"
        
        "adds x4, x4, x5	\n\t"
        "adcs x9, x9, x6	\n\t"
        "eor x10, x10, x10	\n\t"
        "adc x7, x7, x10	\n\t"
    
        "adds x13, x13, x4	\n\t"
        "str x13, [%0, #8 * 21]	\n\t"	//c21
        "adcs x9, x9, x10	\n\t"
        "adc x7, x7, x10	\n\t"
        //step 12
        "ldr x0, [%1, #88]	\n\t"
        "mul x3, x0, x0		\n\t"
        "umulh x4, x0, x0	\n\t"
        "adds x9, x9, x3	\n\t"
        "adcs x7, x7, x4	\n\t"
        "stp x9, x7, [%0, #8 * 22]	\n\t"	//c22, c23	
        :
        :"r"(&res[0]),"r"(&a[0])
        :"memory", "x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14",
        "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26"
       );
        
    }
    
__inline void fpadd751(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular addition, c = a+b mod p751.
  // Inputs: a, b in [0, 2*p751-1] 
  // Output: c in [0, 2*p751-1]

    fpadd751_asm(a, b, c);
} 


__inline void fpsub751(const digit_t* a, const digit_t* b, digit_t* c)
{ // Modular subtraction, c = a-b mod p751.
  // Inputs: a, b in [0, 2*p751-1] 
  // Output: c in [0, 2*p751-1] 

    fpsub751_asm(a, b, c);
}


__inline void fpneg751(digit_t* a)
{ // Modular negation, a = -a mod p751.
  // Input/output: a in [0, 2*p751-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, ((digit_t*)p751x2)[i], a[i], borrow, a[i]); 
    }
}


void fpdiv2_751(const digit_t* a, digit_t* c)
{ // Modular division by two, c = a/2 mod p751.
  // Input : a in [0, 2*p751-1] 
  // Output: c in [0, 2*p751-1] 
    unsigned int i, carry = 0;
    digit_t mask;
        
    mask = 0 - (digit_t)(a[0] & 1);    // If a is odd compute a+p751
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(carry, a[i], ((digit_t*)p751)[i] & mask, carry, c[i]); 
    }

    mp_shiftr1(c, NWORDS_FIELD);
} 


void fpcorrection751(digit_t* a)
{ // Modular correction to reduce field element a in [0, 2*p751-1] to [0, p751-1].
    unsigned int i, borrow = 0;
    digit_t mask;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(borrow, a[i], ((digit_t*)p751)[i], borrow, a[i]);
    }
    mask = 0 - (digit_t)borrow;

    borrow = 0;
    for (i = 0; i < NWORDS_FIELD; i++) {
        ADDC(borrow, a[i], ((digit_t*)p751)[i] & mask, borrow, a[i]);
    }
}


void mp_mul(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{ // Multiprecision multiply, c = a*b, where lng(a) = lng(b) = nwords.

	UNREFERENCED_PARAMETER(nwords);
	//mul751_asm(a, b, c);
	fpmul768_karatsuba(a, b, c);
}

void mp_mul_mixed(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{
 // Multiprecision multiply, c = a*b, where lng(a) = lng(b) = nwords.

	UNREFERENCED_PARAMETER(nwords);
	fpmul768_karatsuba(a, b, c);
}


void rdc_mont(const digit_t* ma, digit_t* mc)
{ // Efficient Montgomery reduction using comba and exploiting the special form of the prime p751.
  // mc = ma*R^-1 mod p751x2, where R = 2^768.
  // If ma < 2^768*p751, the output mc is in the range [0, 2*p751-1].
  // ma is assumed to be in Montgomery representation.
  
    rdc751_asm(ma, mc);
}
