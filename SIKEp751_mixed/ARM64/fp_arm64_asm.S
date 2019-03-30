//*******************************************************************************************
// Supersingular Isogeny Key Encapsulation Library
//
// Abstract: field arithmetic in ARM64 assembly for P751 on Linux
//
//
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

.macro mul_x3_x_x4_to_7
    mul     x19, x3, x4     //                              x19
    umulh   x20, x3, x4     //                          x20
    mul     x21, x3, x5     //                          x21
    umulh   x22, x3, x5     //                      x22
    mul     x23, x3, x6     //                      x23
    umulh   x24, x3, x6     //                  x24
    mul     x25, x3, x7     // +                x25
    umulh   x26, x3, x7     //              x26
.endm

.macro mul_x3_x_x8_to_10
    mul     x19, x3, x8     //              x19
    umulh   x20, x3, x8     //          x20
    mul     x21, x3, x9     //          x21
    umulh   x22, x3, x9     //      x22
    mul     x23, x3, x10    //      x23
.endm

.macro add_move_x11_to_x26
    adds    x11, x12, x19
    adcs    x12, x13, x20
    adcs    x13, x14, x22
    adcs    x14, x15, x24
    adc     x26, x26, xzr

    adds    x12, x12, x21
    adcs    x13, x13, x23
    adcs    x14, x14, x25
    adcs    x15, x16, x26
    adcs    x16, x17, xzr
    adcs    x17, x18, xzr
    adc     x18, xzr, xzr
.endm

.macro  add_x15_to_x23
    adds    x15, x15, x19
    adcs    x16, x16, x20
    adcs    x17, x17, x22
    adcs    x18, xzr, xzr

    adds    x16, x16, x21
    adcs    x17, x17, x23
    adcs    x18, x18, xzr
.endm

.macro add_x12_2_x17_into_x11
    adds    x11, x12, x19
    adcs    x12, x13, x20
    adcs    x13, x14, x22
    adcs    x24, x24, xzr

    adds    x12, x12, x21
    adcs    x13, x13, x23
    adcs    x14, x15, x24
    adcs    x15, x16, xzr
    adcs    x16, x17, xzr
    adcs    x17, xzr, xzr
.endm

.macro add_x14_2_x23
    adds    x14, x14, x19
    adcs    x15, x15, x20
    adcs    x16, x16, x22
    adcs    x24, x24, xzr

    adds    x15, x15, x21
    adcs    x16, x16, x23
    adcs    x17, x17, x24
.endm

.macro mul_x3_x_x7_to_9
    mul     x19, x3, x7     // +                x19
    umulh   x20, x3, x7     //              x20
    mul     x21, x3, x8     //              x21
    umulh   x22, x3, x8     //          x22
    mul     x23, x3, x9     //          x23
    umulh   x24, x3, x9     //      x24
.endm

.macro mul_x3_x_x4_to_6
    mul     x19, x3, x4     //                              x19
    umulh   x20, x3, x4     //                          x20
    mul     x21, x3, x5     //                          x21
    umulh   x22, x3, x5     //                      x22
    mul     x23, x3, x6     //                      x23
    umulh   x24, x3, x6     //                  x24
.endm

//***********************************************************************
// Basic addition
// Operation: c[x2]  = a[x0] + b[x1]
//*********************************************************************
.global mp_add751_asm
mp_add751_asm:
    // load register of a 
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0, #16]
    ldp     x7, x8,   [x0, #32]
    ldp     x9, x10,  [x0, #48]
    // load register of b 
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1, #16]
    ldp     x15, x16, [x1, #32]
    ldp     x17, x18, [x1, #48]
    adds    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adcs    x6,  x6,  x14
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adcs    x10, x10, x18
    stp     x3,  x4,  [x2]
    stp     x5,  x6,  [x2, #16]
    stp     x7,  x8,  [x2, #32]
    stp     x9,  x10, [x2, #48]
    // load register of a 
    ldp     x3, x4,   [x0, #64]
    ldp     x5, x6,   [x0, #80]
    // load register of b 
    ldp     x11, x12, [x1, #64]
    ldp     x13, x14, [x1, #80]
    adcs    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adc     x6,  x6,  x14
    stp     x3,  x4,  [x2, #64]
    stp     x5,  x6,  [x2, #80]
    ret
    
//***********************************************************************
// Basic addition
// Operation: c[x2]  = a[x0] + b[x1]
//*********************************************************************
.global mp_add751x2_asm
mp_add751x2_asm:
    // load register of a 
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0, #16]
    ldp     x7, x8,   [x0, #32]
    ldp     x9, x10,  [x0, #48]
    // load register of b 
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1, #16]
    ldp     x15, x16, [x1, #32]
    ldp     x17, x18, [x1, #48]
    adds    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adcs    x6,  x6,  x14
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adcs    x10, x10, x18
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
    // load register of a 
    ldp     x3,  x4,  [x0, #128]
    ldp     x5,  x6,  [x0, #144]
    ldp     x7,  x8,  [x0, #160]
    ldp     x9,  x10, [x0, #176]
    // load register of b 
    ldp     x11, x12, [x1, #128]
    ldp     x13, x14, [x1, #144]
    ldp     x15, x16, [x1, #160]
    ldp     x17, x18, [x1, #176]
    adcs    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adcs    x6,  x6,  x14
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adc     x10, x10, x18
    stp     x3,  x4,  [x2, #128]
    stp     x5,  x6,  [x2, #144]
    stp     x7,  x8,  [x2, #160]
    stp     x9,  x10, [x2, #176]
    ret

//***********************************************************************
//  2x751-bit multiprecision subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]. Returns borrow mask
//***********************************************************************
.global mp_sub751x2_asm
mp_sub751x2_asm:
    // load register of a
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0, #16]
    ldp     x7, x8,   [x0, #32]
    ldp     x9, x10,  [x0, #48]
    // load register of b
    ldp     x11, x12, [x1]
    ldp     x13, x14, [x1, #16]
    ldp     x15, x16, [x1, #32]
    ldp     x17, x18, [x1, #48]
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
    stp     x3,  x4,  [x2, #64]
    stp     x5,  x6,  [x2, #80]
    stp     x7,  x8,  [x2, #96]
    stp     x9,  x10, [x2, #112]
    // load register of a
    ldp     x3,  x4,  [x0, #128]
    ldp     x5,  x6,  [x0, #144]
    ldp     x7,  x8,  [x0, #160]
    ldp     x9,  x10, [x0, #176]
    // load register of b
    ldp     x11, x12, [x1, #128]
    ldp     x13, x14, [x1, #144]
    ldp     x15, x16, [x1, #160]
    ldp     x17, x18, [x1, #176]
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
    stp     x3,  x4,  [x2, #128]
    stp     x5,  x6,  [x2, #144]
    stp     x7,  x8,  [x2, #160]
    stp     x9,  x10, [x2, #176]
    ret
// 2*p751
.align 16
p751x2: .quad 0xfffffffffffffffe, 0xffffffffffffffff, 0xdd5fffffffffffff, 0xc7d92d0a93f0f151, 0xb52b363427ef98ed,0x109d30cfadd7d0ed, 0x0ac56a08b964ae90, 0x1c25213f2f75b8cd, 0x0000dfcbaa83ee38

//***********************************************************************
//  Field addition
//  Operation: c [reg_p3] = a [reg_p1] + b [reg_p2]
//********************************************************************* 
.global fpadd751_asm
fpadd751_asm:
    // load register of a 
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0, #16]
    ldp     x7, x8,   [x0, #32]
    ldp     x9, x10,  [x0, #48]
    ldp     x11, x12, [x0, #64]
    ldp     x13, x14, [x0, #80]
    // load register of b 
    ldp     x15, x16, [x1]
    ldp     x17, x18, [x1, #16]
    // perform the additions 
    adds    x3,  x3,  x15
    adcs    x4,  x4,  x16
    adcs    x5,  x5,  x17
    adcs    x6,  x6,  x18
    ldp     x15, x16, [x1, #32]
    ldp     x17, x18, [x1, #48]
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adcs    x10, x10, x18
    ldp     x15, x16, [x1, #64]
    ldp     x17, x18, [x1, #80]
    adcs    x11, x11, x15
    adcs    x12, x12, x16
    adcs    x13, x13, x17
    adcs    x14, x14, x18
    // load up the pp1x2 
    adr     x0, p751x2
    ldp     x15, x16, [x0]
    ldp     x17, x18, [x0, #16]
    
    // perform the subtraction 
    subs    x3,  x3,  x15
    sbcs    x4,  x4,  x16
    sbcs    x5,  x5,  x16
    sbcs    x6,  x6,  x16
    sbcs    x7,  x7,  x16
    sbcs    x8,  x8,  x17
    sbcs    x9,  x9,  x18
    ldp     x15, x16, [x0, #32]
    ldp     x17, x18, [x0, #48]
    ldr     x1, [x0, #64]
    
    sbcs    x10, x10, x15
    sbcs    x11, x11, x16
    sbcs    x12, x12, x17
    sbcs    x13, x13, x18
    sbcs    x14, x14, x1
    sbcs    x1,  xzr, xzr

    ldp     x15, x16, [x0]
    ldp     x17, x18, [x0, #16]
    and     x15, x15, x1
    and     x16, x16, x1
    and     x17, x17, x1
    and     x18, x18, x1
    
    adds    x3,  x3,  x15
    adcs    x4,  x4,  x16
    adcs    x5,  x5,  x16
    adcs    x6,  x6,  x16
    adcs    x7,  x7,  x16
    adcs    x8,  x8,  x17
    adcs    x9,  x9,  x18
    
    stp     x3,  x4,  [x2]
    stp     x5,  x6,  [x2, #16]
    stp     x7,  x8,  [x2, #32]
    
    ldp     x15, x16, [x0, #32]
    ldp     x17, x18, [x0, #48]
    ldr     x3,  [x0, #64]
    and     x15, x15, x1
    and     x16, x16, x1
    and     x17, x17, x1
    and     x18, x18, x1
    and     x3,  x3,  x1
    
    adcs    x10, x10, x15
    adcs    x11, x11, x16
    adcs    x12, x12, x17
    adcs    x13, x13, x18
    adcs    x14, x14, x3

    // store the result 
    stp     x9,  x10, [x2, #48]
    stp     x11, x12, [x2, #64]
    stp     x13, x14, [x2, #80]
    ret
    
//***********************************************************************
//  Field subtraction
//  Operation: c [reg_p3] = a [reg_p1] - b [reg_p2]
//********************************************************************* 
.global fpsub751_asm
fpsub751_asm:
    // load register of a 
    ldp     x3, x4,   [x0]
    ldp     x5, x6,   [x0, #16]
    ldp     x7, x8,   [x0, #32]
    ldp     x9, x10,  [x0, #48]
    ldp     x11, x12, [x0, #64]
    ldp     x13, x14, [x0, #80]
    // load register of b 
    ldp     x15, x16, [x1]
    ldp     x17, x18, [x1, #16]
    // perform the subtraction 
    subs    x3,  x3,  x15
    sbcs    x4,  x4,  x16
    sbcs    x5,  x5,  x17
    sbcs    x6,  x6,  x18
    ldp     x15, x16, [x1, #32]
    ldp     x17, x18, [x1, #48]
    sbcs    x7,  x7,  x15
    sbcs    x8,  x8,  x16
    sbcs    x9,  x9,  x17
    sbcs    x10, x10, x18
    ldp     x15, x16, [x1, #64]
    ldp     x17, x18, [x1, #80]
    sbcs    x11, x11, x15
    sbcs    x12, x12, x16
    sbcs    x13, x13, x17
    sbcs    x14, x14, x18
    sbcs    x1,  xzr, xzr
    // load up the p751x2 and and adds it up 
    adr     x0, p751x2
    ldp     x15, x16, [x0]
    ldp     x17, x18, [x0, #16]
    and     x15, x15, x1
    and     x16, x16, x1
    and     x17, x17, x1
    and     x18, x18, x1
    adds    x3,  x3,  x15
    adcs    x4,  x4,  x16
    adcs    x5,  x5,  x16
    adcs    x6,  x6,  x16
    adcs    x7,  x7,  x16
    adcs    x8,  x8,  x17
    adcs    x9,  x9,  x18
    stp     x3,  x4,  [x2]
    stp     x5,  x6,  [x2, #16]
    stp     x7,  x8,  [x2, #32]
    // do next half of load p751x2 and and and add 
    ldp     x15, x16, [x0, #32]
    ldp     x17, x18, [x0, #48]
    ldr     x3,  [x0, #64]
    and     x15, x15, x1
    and     x16, x16, x1
    and     x17, x17, x1
    and     x18, x18, x1
    and     x3,  x3,  x1
    adcs    x10, x10, x15
    adcs    x11, x11, x16
    adcs    x12, x12, x17
    adcs    x13, x13, x18
    adcs    x14, x14, x3
    stp     x9,  x10, [x2, #48]
    stp     x11, x12, [x2, #64]
    stp     x13, x14, [x2, #80]
    ret
 
// p751+1
.align 16
p751p1: .quad  0xeeb0000000000000, 0xe3ec968549f878a8, 0xda959b1a13f7cc76, 0x084e9867d6ebe876, 0x8562b5045cb25748, 0x0e12909f97badc66, 0x00006fe5d541f71c

//***********************************************************************
//  Montgomery reduction
//  Based on comba method
//  Operation: c [reg_p2] = a [reg_p1]
//  NOTE: a=c is not allowed
//********************************************************************* 
.global rdc751_asm
rdc751_asm:
    preserve_caller_registers

    adr     x28, p751p1
    ldp     x22, x23, [x28]
    ldp     x24, x25, [x28, #16]
    ldp     x26, x27, [x28, #32]
    ldr     x28, [x28, #48]

    ldp     x2,  x3,  [x0]
    ldp     x4,  x5,  [x0, #16]
    ldp     x6,  x7,  [x0, #32]
    ldp     x8,  x9,  [x0, #48]
    ldp     x20, x21, [x0, #64]

    mul     x10, x2,  x22
    umulh   x11, x2,  x22
    adds    x15, x10, x7          // x15 = c5
    adcs    x10, x11, xzr

    mul     x13, x2,  x23
    umulh   x14, x2,  x23
    adds    x10, x10, x13
    adcs    x11, x14, xzr

    mul     x13, x3,  x22
    umulh   x14, x3,  x22
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    adds    x16, x10, x8        // x16 = c6
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x2,  x24
    umulh   x14, x2,  x24
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x3,  x23
    umulh   x14, x3,  x23
    add_13_14_to_10_12

    mul     x13, x4,  x22
    umulh   x14, x4,  x22
    add_13_14_to_10_12

    adds    x17,  x10, x9        //  x17 = c7
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x2,  x25
    umulh   x14, x2,  x25
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x3,  x24
    umulh   x14, x3,  x24
    add_13_14_to_10_12

    mul     x13, x4,  x23
    umulh   x14, x4,  x23
    add_13_14_to_10_12

    mul     x13, x5,  x22
    umulh   x14, x5,  x22
    add_13_14_to_10_12

    adds    x18, x10, x20        //  x18 = c8
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    ldp     x8,  x9,  [x0, #80]

    mul     x13, x2,  x26
    umulh   x14, x2,  x26
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x3,  x25
    umulh   x14, x3,  x25
    add_13_14_to_10_12

    mul     x13, x4,  x24
    umulh   x14, x4,  x24
    add_13_14_to_10_12

    mul     x13, x5,  x23
    umulh   x14, x5,  x23
    add_13_14_to_10_12

    mul     x13, x6,  x22
    umulh   x14, x6,  x22
    add_13_14_to_10_12

    adds    x19, x10, x21        // x19 = c9
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x2,  x27
    umulh   x14, x2,  x27
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x3,  x26
    umulh   x14, x3,  x26
    add_13_14_to_10_12

    mul     x13, x4,  x25
    umulh   x14, x4,  x25
    add_13_14_to_10_12

    mul     x13, x5,  x24
    umulh   x14, x5,  x24
    add_13_14_to_10_12

    mul     x13, x6,  x23
    umulh   x14, x6,  x23
    add_13_14_to_10_12

    mul     x13, x15,  x22
    umulh   x14, x15,  x22
    add_13_14_to_10_12
    adds    x20, x10, x8        // x20 = c10
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x2,  x28
    umulh   x14, x2,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x3,  x27
    umulh   x14, x3,  x27
    add_13_14_to_10_12

    mul     x13, x4,  x26
    umulh   x14, x4,  x26
    add_13_14_to_10_12

    mul     x13, x5,  x25
    umulh   x14, x5,  x25
    add_13_14_to_10_12

    mul     x13, x6,  x24
    umulh   x14, x6,  x24
    add_13_14_to_10_12

    mul     x13, x15,  x23
    umulh   x14, x15,  x23
    add_13_14_to_10_12

    mul     x13, x16, x22
    umulh   x14, x16, x22
    add_13_14_to_10_12
    adds    x21, x10, x9        // x21 = c11 x0[11]
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x3,  x28
    umulh   x14, x3,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x4,  x27
    umulh   x14, x4,  x27
    add_13_14_to_10_12

    mul     x13, x5,  x26
    umulh   x14, x5,  x26
    add_13_14_to_10_12

    mul     x13, x6,  x25
    umulh   x14, x6,  x25
    add_13_14_to_10_12

    mul     x13, x15, x24
    umulh   x14, x15, x24
    add_13_14_to_10_12

    mul     x13, x16, x23
    umulh   x14, x16, x23
    add_13_14_to_10_12

    ldp     x8,  x9,  [x0, #96]

    mul     x13, x17,  x22
    umulh   x14, x17,  x22
    add_13_14_to_10_12
    adds    x8, x10, x8        // x8 = c0
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x4,  x28
    umulh   x14, x4,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x5,  x27
    umulh   x14, x5,  x27
    add_13_14_to_10_12

    mul     x13, x6,  x26
    umulh   x14, x6,  x26
    add_13_14_to_10_12

    mul     x13, x15,  x25
    umulh   x14, x15,  x25
    add_13_14_to_10_12
    
    mul     x13, x16,  x24
    umulh   x14, x16,  x24
    add_13_14_to_10_12

    mul     x13, x17,   x23
    umulh   x14, x17,   x23
    add_13_14_to_10_12

    mul     x13, x18,   x22
    umulh   x14, x18,   x22
    add_13_14_to_10_12

    adds    x9, x10, x9       // x9 = c1
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    stp     x8,  x9, [x1]
    ldp     x8,  x9, [x0, #112]

    mul     x13, x5,  x28
    umulh   x14, x5,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x6,  x27
    umulh   x14, x6,  x27
    add_13_14_to_10_12

    mul     x13, x15, x26
    umulh   x14, x15, x26
    add_13_14_to_10_12

    mul     x13, x16, x25
    umulh   x14, x16, x25
    add_13_14_to_10_12
    
    mul     x13, x17,  x24
    umulh   x14, x17,  x24
    add_13_14_to_10_12

    mul     x13, x18,  x23
    umulh   x14, x18,  x23
    add_13_14_to_10_12

    mul     x13, x19,  x22
    umulh   x14, x19,  x22
    add_13_14_to_10_12
    adds    x8, x10, x8       // x8 = c2
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x6,  x28
    umulh   x14, x6,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x15, x27
    umulh   x14, x15, x27
    add_13_14_to_10_12

    mul     x13, x16, x26
    umulh   x14, x16, x26
    add_13_14_to_10_12

    mul     x13, x17,  x25
    umulh   x14, x17,  x25
    add_13_14_to_10_12

    mul     x13, x18,  x24
    umulh   x14, x18,  x24
    add_13_14_to_10_12

    mul     x13, x19,  x23
    umulh   x14, x19,  x23
    add_13_14_to_10_12

    mul     x13, x20,  x22
    umulh   x14, x20,  x22
    add_13_14_to_10_12
    adds    x9, x10, x9       // x9 = c3
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    stp     x8, x9, [x1, #16]
    ldp     x8,  x9,  [x0, #128]

    mul     x13, x15, x28
    umulh   x14, x15, x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x16, x27
    umulh   x14, x16, x27
    add_13_14_to_10_12

    mul     x13, x17,  x26
    umulh   x14, x17,  x26
    add_13_14_to_10_12

    mul     x13, x18,  x25
    umulh   x14, x18,  x25
    add_13_14_to_10_12

    mul     x13, x19,  x24
    umulh   x14, x19,  x24
    add_13_14_to_10_12

    mul     x13, x20,  x23
    umulh   x14, x20,  x23
    add_13_14_to_10_12

    mul     x13, x21,  x22
    umulh   x14, x21,  x22
    add_13_14_to_10_12
    adds    x8, x10, x8       // x8 = c4
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x16, x28
    umulh   x14, x16, x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x17,  x27
    umulh   x14, x17,  x27
    add_13_14_to_10_12

    mul     x13, x18,  x26
    umulh   x14, x18,  x26
    add_13_14_to_10_12

    mul     x13, x19,  x25
    umulh   x14, x19,  x25
    add_13_14_to_10_12

    mul     x13, x20,  x24
    umulh   x14, x20,  x24
    add_13_14_to_10_12

    mul     x13, x21,  x23
    umulh   x14, x21,  x23
    add_13_14_to_10_12
    adds    x9, x10, x9       //  x9 = c5
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    stp     x8, x9, [x1, #32]
    ldp     x8,  x9, [x0, #144]    // START HERE

    mul     x13, x17,  x28
    umulh   x14, x17,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x18,  x27
    umulh   x14, x18,  x27
    add_13_14_to_10_12

    mul     x13, x19,  x26
    umulh   x14, x19,  x26
    add_13_14_to_10_12

    mul     x13, x20,  x25
    umulh   x14, x20,  x25
    add_13_14_to_10_12

    mul     x13, x21,  x24
    umulh   x14, x21,  x24
    add_13_14_to_10_12
    adds    x8, x10, x8       // x8 = c6
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x18,  x28
    umulh   x14, x18,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x19,  x27
    umulh   x14, x19,  x27
    add_13_14_to_10_12

    mul     x13, x20,  x26
    umulh   x14, x20,  x26
    add_13_14_to_10_12

    mul     x13, x21,  x25
    umulh   x14, x21,  x25
    add_13_14_to_10_12
    adds    x9, x10, x9       // x9 = c7
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    stp     x8, x9, [x1, #48]
    ldp     x8,  x9,  [x0, #160]

    mul     x13, x19,  x28
    umulh   x14, x19,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x20,  x27
    umulh   x14, x20,  x27
    add_13_14_to_10_12

    mul     x13, x21,  x26
    umulh   x14, x21,  x26
    add_13_14_to_10_12
    adds    x8, x10, x8       // x8 = c8
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    mul     x13, x20,  x28
    umulh   x14, x20,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr

    mul     x13, x21,  x27
    umulh   x14, x21,  x27
    add_13_14_to_10_12
    adds    x9, x10, x9       // x9 = c9
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    stp     x8, x9, [x1, #64]
    ldp     x8, x9, [x0, #176]

    mul     x13, x21,  x28
    umulh   x14, x21,  x28
    adds    x10, x10, x13
    adcs    x11, x11, x14
    adcs    x12, xzr, xzr
    adds    x8, x10, x8       // x8 = c10
    adcs    x10, x11, xzr
    adcs    x11, x12, xzr

    add     x9, x10, x9       //  x9 = c11

    stp     x8,  x9,  [x1, #80]

    restore_caller_registers
    ret
