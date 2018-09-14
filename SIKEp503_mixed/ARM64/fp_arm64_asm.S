//*******************************************************************************************
// Supersingular Isogeny Key Encapsulation Library
//
// Abstract: field arithmetic in ARM64 assembly for P503 on Linux
//*******************************************************************************************
.text

.macro preserve_caller_registers
    sub   sp,  sp, #80
    stp   x19, x20, [sp]
    stp   x21, x22, [sp, #16]
    stp   x23, x24, [sp, #32]
    stp   x25, x26, [sp, #48]
    stp   x27, x28, [sp, #64]
.endm

// restore_caller_registers(): Restore x18-x28 on stack (sp)
.macro restore_caller_registers
    ldp   x19, x20, [sp]
    ldp   x21, x22, [sp, #16]
    ldp   x23, x24, [sp, #32]
    ldp   x25, x26, [sp, #48]
    ldp   x27, x28, [sp, #64]
    add   sp,  sp,  #80
.endm

.macro add_13_14_to_10_12
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, x12, xzr
.endm

.macro p503_mul_x3_x_x4_to_8
    mul     x18, x3, x4     //                              x18
    umulh   x19, x3, x4     //                          x19
    mul     x20, x3, x5     //                          x20
    umulh   x21, x3, x5     //                      x21
    mul     x22, x3, x6     //                      x22
    umulh   x23, x3, x6     //                  x23
    mul     x24, x3, x7     // +                x24
    umulh   x25, x3, x7     //              x25
    mul     x26, x3, x8     //              x26
.endm


.macro  p503_add_move_x9_to_x26
    adds    x9,  x10, x18
    adcs    x10, x11, x19
    adcs    x11, x12, x21
    adcs    x12, x13, x23
    adcs    x13, x14, x25
    adcs    x14, xzr, xzr

    adds    x10, x10, x20
    adcs    x11, x11, x22
    adcs    x12, x12, x24
    adcs    x13, x13, x26
    adcs    x14, x14, xzr
.endm

.macro p503_mul_x3_x_x4_to_7
    mul     x18, x3, x4     //                              x18
    umulh   x19, x3, x4     //                          x19
    mul     x20, x3, x5     //                          x20
    umulh   x21, x3, x5     //                      x21
    mul     x22, x3, x6     //                      x22
    umulh   x23, x3, x6     //                  x23
    mul     x24, x3, x7     // +                x24
    umulh   x25, x3, x7     //              x25
.endm


.macro  p503_add_move_x9_to_x25
    adds    x9,  x10, x18
    adcs    x10, x11, x20
    adcs    x11, x12, x22
    adcs    x12, x13, x24
    adcs    x13, x25, xzr

    adds    x10, x10, x19
    adcs    x11, x11, x21
    adcs    x12, x12, x23
    adcs    x13, x13, xzr
.endm

.macro  p503_add_x3_x10_and_x11_x18
    adds    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adcs    x6,  x6,  x14
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adcs    x10, x10, x18
.endm

.macro  load_a
    ldp     x3,  x4,  [x0]
    ldp     x5,  x6,  [x0, #16]
    ldp     x7,  x8,  [x0, #32]
    ldp     x9,  x10, [x0, #48]
.endm

.macro  load_b
    // load register of b 
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1, #16]
    ldp     x15, x16, [x1, #32]
    ldp     x17, x18, [x1, #48]
.endm

//***********************************************************************
// Basic addition
// Operation: c[x2]  = a[x0] + b[x1]
//**********************************************************************/
.global mp_add503_asm
mp_add503_asm:
    // load register of a into x3-x10
    load_a
    // load register of b into x11-x18
    load_b
    p503_add_x3_x10_and_x11_x18
    stp     x3,  x4,  [x2]
    stp     x5,  x6,  [x2, #16]
    stp     x7,  x8,  [x2, #32]
    stp     x9,  x10, [x2, #48]
    ret

//***********************************************************************
// Basic addition
// Operation: c[x2]  = a[x0] + b[x1]
//**********************************************************************/
.global mp_add503x2_asm
mp_add503x2_asm:
    // load register of a into x3-x10
    load_a
    // load register of b into x11-x18
    load_b
    p503_add_x3_x10_and_x11_x18
    stp     x3,  x4,  [x2]
    stp     x5,  x6,  [x2, #16]
    stp     x7,  x8,  [x2, #32]
    stp     x9,  x10, [x2, #48]
    // load register of a 
    ldp     x3,  x4,  [x0, #64]
    ldp     x5,  x6,  [x0, #80]
    ldp     x7,  x8,  [x0, #96]
    ldp     x9,  x10, [x0, #112]
    // load register of b 
    ldp     x11, x12, [x1, #64]
    ldp     x13, x14, [x1, #80]
    ldp     x15, x16, [x1, #96]
    ldp     x17, x18, [x1, #112]
    adcs    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adcs    x6,  x6,  x14
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adcs    x10, x10, x18
    stp     x3,  x4,  [x2, #64]
    stp     x5,  x6,  [x2, #80]
    stp     x7,  x8,  [x2, #96]
    stp     x9,  x10, [x2, #112]
    ret

//***********************************************************************
//  2x503-bit multiprecision subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]. Returns borrow mask
//***********************************************************************
.global mp_sub503x2_asm
mp_sub503x2_asm:
    // load register of a into x3-x10
    load_a
    // load register of b into x11-x18
    load_b
    subs    x3,  x3,  x11
    sbcs    x4,  x4,  x12
    sbcs    x5,  x5,  x13
    sbcs    x6,  x6,  x14
    sbcs    x7,  x7,  x15
    sbcs    x8,  x8,  x16
    sbcs    x9,  x9,  x17
    sbcs    x10, x10, x18
    stp     x3,  x4,  [x2]
    stp     x5,  x6,  [x2, #16]
    stp     x7,  x8,  [x2, #32]
    stp     x9,  x10, [x2, #48]
    // load register of a
    ldp     x3,  x4,  [x0, #64]
    ldp     x5,  x6,  [x0, #80]
    ldp     x7,  x8,  [x0, #96]
    ldp     x9,  x10, [x0, #112]
    // load register of b
    ldp     x11, x12, [x1, #64]
    ldp     x13, x14, [x1, #80]
    ldp     x15, x16, [x1, #96]
    ldp     x17, x18, [x1, #112]
    sbcs    x3,  x3,  x11
    sbcs    x4,  x4,  x12
    sbcs    x5,  x5,  x13
    sbcs    x6,  x6,  x14
    sbcs    x7,  x7,  x15
    sbcs    x8,  x8,  x16
    sbcs    x9,  x9,  x17
    sbcs    x10, x10, x18
    mov     x0,  #0
    sbc     x0,  x0,  x0
    stp     x3,  x4,  [x2, #64]
    stp     x5,  x6,  [x2, #80]
    stp     x7,  x8,  [x2, #96]
    stp     x9,  x10, [x2, #112]
    ret

// 2*p503
.align 16
p503x2: .quad 0xfffffffffffffffe, 0xffffffffffffffff, 0x57ffffffffffffff, 0x2610b7b44423cf41, 0x3737ed90f6fcfb5e,0xc08b8d7bb4ef49a0, 0x0080cdea83023c3c
//***********************************************************************
//  Field addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2]
//***********************************************************************
.global fpadd503_asm
fpadd503_asm:
    // load register of a into x3-x10
    load_a
    // load register of b into x11-x18
    load_b
    p503_add_x3_x10_and_x11_x18
    // load up the p503x2
    adr     x18, p503x2
    ldp     x11, x12, [x18]
    ldp     x13, x14, [x18, #16]
    ldp     x15, x16, [x18, #32]
    ldr     x17, [x18, #48]
    // perform the subtraction
    subs    x3,  x3,  x11
    sbcs    x4,  x4,  x12
    sbcs    x5,  x5,  x12
    sbcs    x6,  x6,  x13
    sbcs    x7,  x7,  x14
    sbcs    x8,  x8,  x15
    sbcs    x9,  x9,  x16
    sbcs    x10, x10, x17
    sbcs    x18, xzr, xzr

    and     x11, x11, x18
    and     x12, x12, x18
    and     x13, x13, x18
    and     x14, x14, x18
    and     x15, x15, x18
    and     x16, x16, x18
    and     x17, x17, x18

    adds    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x12
    adcs    x6,  x6,  x13
    adcs    x7,  x7,  x14
    adcs    x8,  x8,  x15
    adcs    x9,  x9,  x16
    adcs    x10, x10, x17

    // store the result 
    stp     x3,  x4,  [x2]
    stp     x5,  x6,  [x2, #16]
    stp     x7,  x8,  [x2, #32]
    stp     x9,  x10, [x2, #48]
    // ok done 
    ret

//***********************************************************************
//  Field subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]
//*********************************************************************** 
.global fpsub503_asm
fpsub503_asm:
    // load register of a into x3-x10
    load_a
    // load register of b into x11-x18
    load_b
    // perform the subtraction 
    subs    x3,  x3,  x11
    sbcs    x4,  x4,  x12
    sbcs    x5,  x5,  x13
    sbcs    x6,  x6,  x14
    sbcs    x7,  x7,  x15
    sbcs    x8,  x8,  x16
    sbcs    x9,  x9,  x17
    sbcs    x10, x10, x18
    sbcs    x1,  xzr, xzr
    // load up the p503x2 
    adr     x18, p503x2
    ldp     x11, x12, [x18]
    ldp     x13, x14, [x18, #16]
    ldp     x15, x16, [x18, #32]
    ldr     x17, [x18, #48]

    and     x11, x11, x1
    and     x12, x12, x1
    and     x13, x13, x1
    and     x14, x14, x1
    and     x15, x15, x1
    and     x16, x16, x1
    and     x17, x17, x1
    // add these values back in
    adds    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x12
    adcs    x6,  x6,  x13
    adcs    x7,  x7,  x14
    adcs    x8,  x8,  x15
    adcs    x9,  x9,  x16
    adcs    x10, x10, x17
    // do the storing 
    stp     x3,  x4,  [x2], #16
    stp     x5,  x6,  [x2], #16
    stp     x7,  x8,  [x2], #16
    stp     x9,  x10, [x2], #16
    // and get out 
    ret

// p503+1
.align 16
p503p1: .quad  0xac00000000000000, 0x13085bda2211e7a0, 0x1b9bf6c87b7e7daf, 0x6045c6bdda77a4d0, 0x004066f541811e1e
//***********************************************************************
//  Montgomery reduction
//  Based on comba method
//  Operation: c [reg_p2] = a [reg_p1]
//  NOTE: a=c is not allowed
//*********************************************************************** 
.global rdc503_asm
rdc503_asm:
    preserve_caller_registers

    // load p503-1 into x[15-19]
    adr     x20, p503p1
    ldp     x15, x16, [x20]
    ldp     x17, x18, [x20, #16]
    ldr     x19, [x20, #32]

    ldp     x2,  x3,  [x0]          // load a[0-7] 
    ldp     x4,  x5,  [x0, #16]
    ldp     x6,  x7,  [x0, #32]
    ldp     x8,  x9,  [x0, #48]

    // i = 3 to 7, (z = 3), (s = 8)
    // j = 0 < i - z + 1 = 1       
    mul     x10, x2,  x15
    umulh   x11, x2,  x15
    adds    x23, x10, x5          // x23 = c3
    adcs    x10, x11, xzr
    // i = 4
    // if (j < i - z + 1 = 2)
    //      (t, u, v) = c_j * p'i-j + (t, u, v)
    //  (t, u, v) = (t, u, v) + a_i
    //  c_i  = v; v = u, u = t, t = 0 */
    mul     x13, x2,  x16
    umulh   x14, x2,  x16
    adds    x10, x10, x13
    adcs    x11, x14, xzr
    mul     x13, x3,  x15
    umulh   x14, x3,  x15
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    adds    x24, x10, x6        // x24 = c4
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    // i = 5
    //  if (j < i - z + 1 = 3)
    //          (t, u, v) = c_j * p'i-j + (t, u, v)
    //  (t, u, v) = (t, u, v) + a_i
    //  c_i = v; v = u, u = t, t = 0 */
    mul     x13, x2,  x17
    umulh   x14, x2,  x17
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    mul     x13, x3,  x16
    umulh   x14, x3,  x16
    add_13_14_to_10_12
    mul     x13, x4,  x15
    umulh   x14, x4,  x15
    add_13_14_to_10_12
    adds    x25, x10, x7        //  x25 = c5
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    // i = 6
    // if (j < i - z + 1 = 4)
    //      (t, u, v) = c_j * p'i-j + (t, u, v)
    //  (t, u, v) = (t, u, v) + a_i
    //  c_i = v; v = u, u = t, t = 0 */
    mul     x13, x2,  x18
    umulh   x14, x2,  x18
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    mul     x13, x3,  x17
    umulh   x14, x3,  x17
    add_13_14_to_10_12
    mul     x13, x4,  x16
    umulh   x14, x4,  x16
    add_13_14_to_10_12
    mul     x13, x23, x15
    umulh   x14, x23, x15
    add_13_14_to_10_12
    adds    x26, x10, x8        //  x26 = c6
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    // i = 7
    // if (j < i - z + 1 = 4)
    //      (t, u, v) = c_j * p'i-j + (t, u, v)
    //  (t, u, v) = (t, u, v) + a_i
    //  c_i  = v; v = u, u = t, t = 0 */
    mul     x13, x2,  x19
    umulh   x14, x2,  x19
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    mul     x13, x3,  x18
    umulh   x14, x3,  x18
    add_13_14_to_10_12
    mul     x13, x4,  x17
    umulh   x14, x4,  x17
    add_13_14_to_10_12
    mul     x13, x23, x16
    umulh   x14, x23, x16
    add_13_14_to_10_12
    mul     x13, x24, x15
    umulh   x14, x24, x15
    add_13_14_to_10_12
    adds    x27, x10, x9        // x27 = c7
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    // for i = 8  to 14 ,
    //  i = 8
    //  z = z - 1 = 2
    //  for j == i - 8 + 1  to s - 1 = 7
    //      if (j < 8 - z)
    //          (t, u, v) = c_j * p'i-j + (t, u, v)
    //  (t, u, v) = (t, u, v) + a_i
    //  c_i-s = v; v = u, u = t, t = 0
    
    // z = 2, j = 8 - 8 + 1 = 1 to 7 */
    // j = 1 < 6 
    mul     x13, x3,  x19   // i - j = 7 
    umulh   x14, x3,  x19
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    // j = 2 
    mul     x13, x4,  x18
    umulh   x14, x4,  x18
    add_13_14_to_10_12
    // j = 3 
    mul     x13, x23, x17
    umulh   x14, x23, x17
    add_13_14_to_10_12
    // j = 4 
    mul     x13, x24,  x16
    umulh   x14, x24,  x16
    add_13_14_to_10_12
    // j = 5 
    ldp     x2, x3, [x0, #64]   // load a[8][9]
    mul     x13, x25, x15
    umulh   x14, x25, x15
    add_13_14_to_10_12
    adds    x20, x10, x2        // x20 = c0
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    // i = 9 
    // z = 1 
    // j = 9 - 8 + 1 = 2 < 8 - 1 = 7 
    // j = 2 
    mul     x13, x4,  x19   // i - j =
    umulh   x14, x4,  x19
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    // j = 3 
    mul     x13, x23, x18
    umulh   x14, x23, x18
    add_13_14_to_10_12
    // j = 4 
    mul     x13, x24, x17
    umulh   x14, x24, x17
    add_13_14_to_10_12
    // j = 5 
    mul     x13, x25, x16
    umulh   x14, x25, x16
    add_13_14_to_10_12
    // j = 6 
    mul     x13, x26, x15
    umulh   x14, x26, x15
    add_13_14_to_10_12
    adds    x21, x10, x3        // x21 = c1
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    ldp     x4, x5, [x0, #80]   // a[10-11]
    // i = 10 
    // z = 0 
    // j = 10 - 8 + 1 = 3 < 8 
    // j = 3 
    mul     x13, x23, x19   // i - j = 7
    umulh   x14, x23, x19
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    // j = 4 
    mul     x13, x24, x18
    umulh   x14, x24, x18
    add_13_14_to_10_12
    // j = 5 
    mul     x13, x25, x17
    umulh   x14, x25, x17
    add_13_14_to_10_12
    // j = 6 
    mul     x13, x26, x16
    umulh   x14, x26, x16
    add_13_14_to_10_12
    // j = 7 
    mul     x13, x27, x15
    umulh   x14, x27, x15
    add_13_14_to_10_12
    adds    x22, x10, x4        // x22 = c2
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    // i = 11 
    // z = 0 
    // j = 11 - 8 + 1 = 4 < 8 
    // j = 4 
    mul     x13, x24, x19 // i - j = 7
    umulh   x14, x24, x19
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    // j = 5 
    mul     x13, x25, x18
    umulh   x14, x25, x18
    add_13_14_to_10_12
    // j = 6 
    mul     x13, x26, x17
    umulh   x14, x26, x17
    add_13_14_to_10_12
    // j = 7 
    mul     x13, x27, x16
    umulh   x14, x27, x16
    add_13_14_to_10_12
    adds    x23, x10, x5       // x23 = c3
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    ldp     x6, x7, [x0, #96]   // a[12-13]
    // i = 12 
    // z = 0 
    // j = 12 - 8 + 1 = 5 < 8 
    // j = 5 
    mul     x13, x25, x19   // i - j = 7
    umulh   x14, x25, x19
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    /// j = 6 
    mul     x13, x26, x18
    umulh   x14, x26, x18
    add_13_14_to_10_12
    // j = 7 
    mul     x13, x27, x17
    umulh   x14, x27, x17
    add_13_14_to_10_12
    adds    x24, x10, x6       // x24 = c4
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    // i = 13 
    // z = 0 
    // j = 13 - 8 + 1 = 6 < 8 
    // j = 6 
    mul     x13, x26, x19      // i - j = 7
    umulh   x14, x26, x19
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    // j = 7 
    mul     x13, x27, x18
    umulh   x14, x27, x18
    add_13_14_to_10_12
    adds    x25, x10, x7       // x25 = c5
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr
    ldp     x8, x9, [x0, #112]   // x0[14-15]
    // i = 14 
    // z = 0 
    // j = 14 - 8 + 1 = 7 < 8 
    // j = 7 
    mul     x13, x27, x19      // i - j = 7
    umulh   x14, x27, x19
    adds    x10, x10, x13
    adcs    x14, x14, xzr
    adds    x26, x10, x8       // x26 = c6
    adcs    x10, x11, x14
    adcs    x11, x12, xzr

    add     x27, x10, x9       //  x27 = c7

    stp     x20, x21, [x1]
    stp     x22, x23, [x1, #16]
    stp     x24, x25, [x1, #32]
    stp     x26, x27, [x1, #48]

    restore_caller_registers
    ret

