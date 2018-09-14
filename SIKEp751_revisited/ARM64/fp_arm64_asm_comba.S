//**********************************************************************************
//
// Abstract: Assembly optimizations for finite field arithmetic on 64-bit ARM. 
//
// Microsoft Research Code	
// 
//**********************************************************************************


.data
// p751 + 1
p751p1:
/*
#define p751p1_5   0xEEB0000000000000
#define p751p1_6   0xE3EC968549F878A8
#define p751p1_7   0xDA959B1A13F7CC76
#define p751p1_8   0x084E9867D6EBE876
#define p751p1_9   0x8562B5045CB25748
#define p751p1_10  0x0E12909F97BADC66
#define p751p1_11  0x00006FE5D541F71C
*/
.quad   0xEEB0000000000000
.quad   0xE3EC968549F878A8
.quad   0xDA959B1A13F7CC76
.quad   0x084E9867D6EBE876
.quad   0x8562B5045CB25748
.quad   0x0E12909F97BADC66
.quad   0x00006FE5D541F71C

//  p751
#define p751_0     0xFFFFFFFFFFFFFFFF
#define p751_5     0xEEAFFFFFFFFFFFFF
#define p751_6     0xE3EC968549F878A8
#define p751_7     0xDA959B1A13F7CC76
#define p751_8     0x084E9867D6EBE876
#define p751_9     0x8562B5045CB25748
#define p751_10    0x0E12909F97BADC66
#define p751_11    0x00006FE5D541F71C

p751:
.quad     0xFFFFFFFFFFFFFFFF
.quad     0xEEAFFFFFFFFFFFFF
.quad     0xE3EC968549F878A8
.quad     0xDA959B1A13F7CC76
.quad     0x084E9867D6EBE876
.quad     0x8562B5045CB25748
.quad     0x0E12909F97BADC66
.quad     0x00006FE5D541F71C
neg_p751:
.quad  0x0000000000000001
.quad  0x0000000000000000
.quad  0x1150000000000000

.quad  0x1c13697ab6078757 
.quad  0x256a64e5ec083389
.quad  0xf7b1679829141789
.quad  0x7a9d4afba34da8b7
.quad  0xf1ed6f6068452399
.quad  0xffff901a2abe08e3



.text

.global fpadd751_arm64_asm
.global fpsub751_arm64_asm
.global fpmul751_arm64_asm
.global fpmul3by3_arm64_asm
.global rdc751_arm64_asm

fpadd751_arm64_asm:
    //  arguments are 3 pointers of type digit_t*, where the first two
    //  arguments are summands and the third is the result register. These arguments
    //  are stored in x0, x1, and x2 respectively.

    //  load first summand into x3 - x14
    ldp x3, x4,   [x0,#0]
    ldp x5, x6,   [x0,#16]
    ldp x7, x8,   [x0,#32]
    ldp x9, x10,  [x0,#48]
    ldp x11, x12, [x0,#64]
    ldp x13, x14, [x0,#80]


    //  add first summand and second summand and store result in x3 - x14
    ldp x15, x16,   [x1,#0]
    ldp x17, x18,   [x1,#16]
    adds x3, x3, x15
    adcs x4, x4, x16
    adcs x5, x5, x17
    adcs x6, x6, x18

    ldp x15, x16,   [x1,#32]
    ldp x17, x18,   [x1,#48]
    adcs x7, x7, x15
    adcs x8, x8, x16
    adcs x9, x9, x17
    adcs x10, x10, x18

    ldp x15, x16,   [x1,#64]
    ldp x17, x18,   [x1,#80]
    adcs x11, x11, x15
    adcs x12, x12, x16
    adcs x13, x13, x17
    adcs x14, x14, x18

    
    //  add the negative of p751 to the resut in x3 - x14 (i.e, subtract p751)

    ldr x15, neg_p751
    ldr x16, neg_p751 + 8
    adds x3, x3, x15
    adcs x4, x4, x16
    adcs x5, x5, x16
    adcs x6, x6, x16
    adcs x7, x7, x16
    ldr x17, neg_p751 + 16
    ldr x18, neg_p751 + 24
    adcs x8, x8, x17
    ldr x15, neg_p751 + 32
    adcs x9, x9, x18
    ldr x16, neg_p751 + 40
    adcs x10, x10, x15
    ldr x17, neg_p751 + 48
    adcs x11, x11, x16
    ldr x18, neg_p751 + 56
    adcs x12, x12, x17
    ldr x15, neg_p751 + 64
    adcs x13, x13, x18
    adcs x14, x14, x15


    sbc x15, xzr, xzr


    //  add p751 back but anded with the mask in x15


    ldr x16, p751
    and x16, x16, x15
    ldr x17, p751 + 8
    and x17, x17, x15
    ldr x18, p751 + 16
    and x18, x18, x15    

    adds x3, x3, x16
    adcs x4, x4, x16
    adcs x5, x5, x16
    adcs x6, x6, x16
    adcs x7, x7, x16
    adcs x8, x8, x17
    adcs x9, x9, x18

    ldr x16, p751 + 24
    and x16, x16, x15
    ldr x17, p751 + 32
    and x17, x17, x15
    ldr x18, p751 + 40
    and x18, x18, x15 

    adcs x10, x10, x16   
    adcs x11, x11, x17   
    adcs x12, x12, x18   

    ldr x16, p751 + 48
    and x16, x16, x15
    ldr x17, p751 + 56
    and x17, x17, x15

    adcs x13, x13, x16
    adcs x14, x14, x17



    stp x3, x4,   [x2,#0]
    stp x5, x6,   [x2,#16]
    stp x7, x8,   [x2,#32]
    stp x9, x10,  [x2,#48]
    stp x11, x12, [x2,#64]
    stp x13, x14, [x2,#80]

    ret


fpsub751_arm64_asm:

    ldp x3, x4,   [x0,#0]
    ldp x5, x6,   [x0,#16]
    ldp x7, x8,   [x0,#32]
    ldp x9, x10,  [x0,#48]
    ldp x11, x12, [x0,#64]
    ldp x13, x14, [x0,#80]

    ldp x15, x16, [x1, #0]
    subs x3, x3, x15
    sbcs x4, x4, x16
    ldp x15, x16, [x1, #16]
    sbcs x5, x5, x15
    sbcs x6, x6, x16
    ldp x15, x16, [x1, #32]
    sbcs x7, x7, x15
    sbcs x8, x8, x16
    ldp x15, x16, [x1, #48]
    sbcs x9, x9, x15
    sbcs x10, x10, x16
    ldp x15, x16, [x1, #64]
    sbcs x11, x11, x15
    sbcs x12, x12, x16
    ldp x15, x16, [x1, #80]
    sbcs x13, x13, x15
    sbcs x14, x14, x16

    sbc x17, xzr, xzr

    
    ldr x15, p751
    and x15, x15, x17
    ldr x16, p751 + 8
    and x16, x16, x17

    adds x3, x3, x15
    adcs x4, x4, x15
    adcs x5, x5, x15
    adcs x6, x6, x15
    adcs x7, x7, x15
    adcs x8, x8, x16

    ldr x15, p751 + 16
    and x15, x15, x17    
    ldr x16, p751 + 24
    and x16, x16, x17

    adcs x9, x9, x15
    adcs x10, x10, x16   

    ldr x15, p751 + 32
    and x15, x15, x17
    ldr x16, p751 + 40
    and x16, x16, x17 

    adcs x11, x11, x15   
    adcs x12, x12, x16   

    ldr x15, p751 + 48
    and x15, x15, x17
    ldr x16, p751 + 56
    and x16, x16, x17

    adcs x13, x13, x15
    adcs x14, x14, x16


    stp x3, x4,   [x2,#0]
    stp x5, x6,   [x2,#16]
    stp x7, x8,   [x2,#32]
    stp x9, x10,  [x2,#48]
    stp x11, x12, [x2,#64]
    stp x13, x14, [x2,#80]


    ret



fpmul751_arm64_asm:

    sub sp, sp, #80
    stp x19, x20, [sp]
    stp x21, x22, [sp, #16]
    stp x23, x24, [sp, #32]
    stp x25, x26, [sp, #48]
    stp x27, x28, [sp, #64]


    ldp x3, x4, [x0, #0]
    ldp x5, x6, [x1, #0]
    mul x18, x3, x5
    umulh x17, x3, x5
    //  c0 is now in x18

    //  a0 * b1
    mul x13, x3, x6
    umulh x14, x3, x6

    adds x17, x17, x13
    adcs x16, x14, xzr
    adcs x15, xzr, xzr

    //  b0 * a1
    mul x13, x4, x5
    umulh x14, x4, x5

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  store c0 and c1
    stp x18, x17, [x2, #0]



    //  load a2, a3, b2, b3
    ldp x7, x8, [x0, #16]
    ldp x9, x10, [x1, #16]

    //  a0 * b2
    mul x13, x3, x9
    umulh x14, x3, x9

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, xzr, xzr

    //  a1 * b1
    mul x13, x4, x6
    umulh x14, x4, x6

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a2 * b0
    mul x13, x7, x5
    umulh x14, x7, x5

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  c2 is now in x16

    //  a0 * b3
    mul x13, x3, x10
    umulh x14, x3, x10

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, xzr, xzr

    //  a1 * b2
    mul x13, x4, x9
    umulh x14, x4, x9

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a2 * b1
    mul x13, x7, x6
    umulh x14, x7, x6

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a3 * b0
    mul x13, x8, x5
    umulh x14, x8, x5

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  store c2 and c3
    stp x16, x15, [x2, #16]


    //  a1 * b3
    mul x13, x4, x10
    umulh x14, x4, x10

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, xzr, xzr

    //  a2 * b2
    mul x13, x7, x9
    umulh x14, x7, x9

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a3 * b1
    mul x13, x8, x6
    umulh x14, x8, x6

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  load a4, a5
    ldp x11, x12, [x0, #32]
    
    //  a4 * b0
    mul x13, x11, x5
    umulh x14, x11, x5

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  load b4, b5
    ldp x19, x20, [x1, #32]

    //  a0 * b4
    mul x13, x3, x19
    umulh x14, x3, x19

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  c4 is now in x18

    //  a0 * b5
    mul x13, x3, x20
    umulh x14, x3, x20

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, xzr, xzr

    //  a1 * b4
    mul x13, x4, x19
    umulh x14, x4, x19

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a2 * b3
    mul x13, x7, x10
    umulh x14, x7, x10

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a3 * b2
    mul x13, x8, x9
    umulh x14, x8, x9

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a4 * b1
    mul x13, x11, x6
    umulh x14, x11, x6

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a5 * b0
    mul x13, x12, x5
    umulh x14, x12, x5

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  store c4 and c5
    stp x18, x17, [x2, #32]

    //  load a6, a7
    ldp x21, x22, [x0, #48]

    //  a6 * b0
    mul x13, x21, x5
    umulh x14, x21, x5

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, xzr, xzr

    //  a5 * b1
    mul x13, x12, x6
    umulh x14, x12, x6

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr
    
    //  a4 * b2
    mul x13, x11, x9
    umulh x14, x11, x9

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a3 * b3
    mul x13, x8, x10
    umulh x14, x8, x10

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a2 * b4
    mul x13, x7, x19
    umulh x14, x7, x19

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a1 * b5
    mul x13, x4, x20
    umulh x14, x4, x20

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  load b6, b7
    ldp x23, x24, [x1, #48]

    //  a0 * b6
    mul x13, x3, x23
    umulh x14, x3, x23

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  c6 is now in x16

    //  a0 * b7
    mul x13, x3, x24
    umulh x14, x3, x24

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, xzr, xzr

    //  a1 * b6
    mul x13, x4, x23
    umulh x14, x4, x23

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a2 * b5
    mul x13, x7, x20
    umulh x14, x7, x20

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a3 * b4
    mul x13, x8, x19
    umulh x14, x8, x19

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a4 * b3
    mul x13, x11, x10
    umulh x14, x11, x10

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a5 * b2
    mul x13, x12, x9
    umulh x14, x12, x9

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a6 * b1
    mul x13, x21, x6
    umulh x14, x21, x6

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a7 * b0
    mul x13, x22, x5
    umulh x14, x22, x5

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  store c6 and c7
    stp x16, x15, [x2, #48]

    //  load a8, a9
    ldp x25, x26, [x0, #64]

    //  a8 * b0
    mul x13, x25, x5
    umulh x14, x25, x5

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, xzr, xzr

    //  a7 * b1
    mul x13, x22, x6
    umulh x14, x22, x6

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a6 * b2
    mul x13, x21, x9
    umulh x14, x21, x9

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a5 * b3
    mul x13, x12, x10
    umulh x14, x12, x10

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a4 * b4
    mul x13, x11, x19
    umulh x14, x11, x19

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a3 * b5
    mul x13, x8, x20
    umulh x14, x8, x20

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a2 * b6
    mul x13, x7, x23
    umulh x14, x7, x23

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a1 * b7
    mul x13, x4, x24
    umulh x14, x4, x24

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  load b8, b9
    ldp x27, x28, [x1, #64]

    //  a0 * b8
    mul x13, x3, x27
    umulh x14, x3, x27

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  c8 is now in x18

    //  a0 * b9
    mul x13, x3, x28
    umulh x14, x3, x28

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, xzr, xzr

    //  a1 * b8
    mul x13, x4, x27
    umulh x14, x4, x27

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a2 * b7
    mul x13, x7, x24
    umulh x14, x7, x24

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a3 * b6
    mul x13, x8, x23
    umulh x14, x8, x23

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a4 * b5
    mul x13, x11, x20
    umulh x14, x11, x20

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a5 * b4
    mul x13, x12, x19
    umulh x14, x12, x19

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a6 * b3
    mul x13, x21, x10
    umulh x14, x21, x10

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a7 * b2
    mul x13, x22, x9
    umulh x14, x22, x9

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a8 * b1
    mul x13, x25, x6
    umulh x14, x25, x6

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a9 * b0
    mul x13, x26, x5
    umulh x14, x26, x5

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  store c8 and c9
    stp x18, x17, [x2, #64]

    //  load a10, a11; a0 and a1 unloaded
    ldp x3, x4, [x0, #80]

    //  a10 * b0
    mul x13, x3, x5
    umulh x14, x3, x5

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, xzr, xzr

    //  a9 * b1
    mul x13, x26, x6
    umulh x14, x26, x6

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a8 * b2
    mul x13, x25, x9
    umulh x14, x25, x9

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a7 * b3
    mul x13, x22, x10
    umulh x14, x22, x10

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a6 * b4
    mul x13, x21, x19
    umulh x14, x21, x19

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a5 * b5
    mul x13, x12, x20
    umulh x14, x12, x20

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a4 * b6
    mul x13, x11, x23
    umulh x14, x11, x23

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a3 * b7
    mul x13, x8, x24
    umulh x14, x8, x24

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a2 * b8
    mul x13, x7, x27
    umulh x14, x7, x27

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  load a0, a1; b0 and b1 unloaded
    ldp x5, x6, [x0, #0]

    //  a1 * b9
    mul x13, x6, x28
    umulh x14, x6, x28

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  load b10, b11; a10 and a11 unloaded
    ldp x3, x4, [x1, #80]

    //  a0 * b10
    mul x13, x3, x5
    umulh x14, x3, x5

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  c10 now in x16

    //  a0 * b11
    mul x13, x4, x5
    umulh x14, x4, x5

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, xzr, xzr

    //  a1 * b10
    mul x13, x3, x6
    umulh x14, x3, x6

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a2 * b9
    mul x13, x7, x28
    umulh x14, x7, x28

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a3 * b8
    mul x13, x8, x27
    umulh x14, x8, x27

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a4 * b7
    mul x13, x11, x24
    umulh x14, x11, x24

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a5 * b6
    mul x13, x12, x23
    umulh x14, x12, x23

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a6 * b5
    mul x13, x21, x20
    umulh x14, x21, x20

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a7 * b4
    mul x13, x22, x19
    umulh x14, x22, x19

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a8 * b3
    mul x13, x25, x10
    umulh x14, x25, x10

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a9 * b2
    mul x13, x26, x9
    umulh x14, x26, x9

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  load a10, a11; b10 and b11 unloaded
    ldp x3, x4, [x0, #80]
    //  load b0, b1; a0 and a1 unloaded
    ldp x5, x6, [x1, #0]

    //  a10 * b1
    mul x13, x3, x6
    umulh x14, x3, x6

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a11 * b0
    mul x13, x4, x5
    umulh x14, x4, x5

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  store c10 and c11
    stp x16, x15, [x2, #80]

    //  a11 * b1
    mul x13, x4, x6
    umulh x14, x4, x6

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, xzr, xzr
    
    //  a10 * b2
    mul x13, x9, x3
    umulh x14, x9, x3

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a9 * b3
    mul x13, x26, x10
    umulh x14, x26, x10

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a8 * b4
    mul x13, x25, x19
    umulh x14, x25, x19

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a7 * b5
    mul x13, x22, x20
    umulh x14, x22, x20

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a6 * b6
    mul x13, x21, x23
    umulh x14, x21, x23

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a5 * b7
    mul x13, x12, x24
    umulh x14, x12, x24

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a4 * b8
    mul x13, x11, x27
    umulh x14, x11, x27

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a3 * b9
    mul x13, x8, x28
    umulh x14, x8, x28

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  load b10, b11; a10 and a11 unloaded
    ldp x3, x4, [x1, #80]
    //  load a0, a1; b0 and b1 unloaded
    ldp x5, x6, [x0, #0]

    //  a2 * b10
    mul x13, x7, x3
    umulh x14, x7, x3

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a1 * b11
    mul x13, x6, x4
    umulh x14, x6, x4

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  c12 now in x18

    //  a2 * b11
    mul x13, x7, x4
    umulh x14, x7, x4

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, xzr, xzr

    //  a3 * b10
    mul x13, x8, x3
    umulh x14, x8, x3

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a4 * b9
    mul x13, x11, x28
    umulh x14, x11, x28

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a5 * b8
    mul x13, x12, x27
    umulh x14, x12, x27

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a6 * b7
    mul x13, x21, x24
    umulh x14, x21, x24

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a7 * b6
    mul x13, x22, x23
    umulh x14, x22, x23

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a8 * b5
    mul x13, x25, x20
    umulh x14, x25, x20

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a9 * b4
    mul x13, x26, x19
    umulh x14, x26, x19

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  load a10, a11; a0 and a1 unloaded
    ldp x5, x6, [x0, #80]

    //  a10 * b3
    mul x13, x5, x10
    umulh x14, x5, x10

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a11 * b2
    mul x13, x6, x9
    umulh x14, x6, x9

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  store c12 and c13
    stp x18, x17, [x2, #96]

    //  a11 * b3
    mul x13, x6, x10
    umulh x14, x6, x10

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, xzr, xzr

    //  a10 * b4
    mul x13, x5, x19
    umulh x14, x5, x19

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a9 * b5
    mul x13, x26, x20
    umulh x14, x26, x20

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a8 * b6
    mul x13, x25, x23
    umulh x14, x25, x23

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a7 * b7
    mul x13, x22, x24
    umulh x14, x22, x24

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a6 * b8
    mul x13, x21, x27
    umulh x14, x21, x27

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a5 * b9
    mul x13, x12, x28
    umulh x14, x12, x28

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a4 * b10
    mul x13, x11, x3
    umulh x14, x11, x3

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a3 * b11
    mul x13, x8, x4
    umulh x14, x8, x4

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  c14 is now in x16

    //  a4 * b11
    mul x13, x11, x4
    umulh x14, x11, x4

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, xzr, xzr

    //  a5 * b10
    mul x13, x12, x3
    umulh x14, x12, x3

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a6 * b9
    mul x13, x21, x28
    umulh x14, x21, x28

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a7 * b8
    mul x13, x22, x27
    umulh x14, x22, x27

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a8 * b7
    mul x13, x25, x24
    umulh x14, x25, x24

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a9 * b6
    mul x13, x26, x23
    umulh x14, x26, x23

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a10 * b5
    mul x13, x5, x20
    umulh x14, x5, x20

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a11 * b4
    mul x13, x6, x19
    umulh x14, x6, x19

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  c15 is now in x15

    //  store c14 and c15
    stp x16, x15, [x2, #112]

    //  a11 * b5
    mul x13, x6, x20
    umulh x14, x6, x20

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, xzr, xzr

    //  a10 * b6
    mul x13, x5, x23
    umulh x14, x5, x23

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a9 * b7
    mul x13, x26, x24
    umulh x14, x26, x24

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a8 * b8
    mul x13, x25, x27
    umulh x14, x25, x27

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a7 * b9
    mul x13, x22, x28
    umulh x14, x22, x28

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a6 * b10
    mul x13, x21, x3
    umulh x14, x21, x3

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a5 * b11
    mul x13, x12, x4
    umulh x14, x12, x4

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  c16 is now in x18

    //  a6 * b11
    mul x13, x21, x4
    umulh x14, x21, x4

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, xzr, xzr

    //  a7 * b10
    mul x13, x22, x3
    umulh x14, x22, x3

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a8 * b9
    mul x13, x25, x28
    umulh x14, x25, x28

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a9 * b8
    mul x13, x26, x27
    umulh x14, x26, x27

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a10 * b7
    mul x13, x5, x24
    umulh x14, x5, x24

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  a11 * b6
    mul x13, x6, x23
    umulh x14, x6, x23

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  store c16 and c17
    stp x18, x17, [x2, #128]

    //  a11 * b7
    mul x13, x6, x24
    umulh x14, x6, x24

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, xzr, xzr

    //  a10 * b8
    mul x13, x5, x27
    umulh x14, x5, x27

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a9 * b9
    mul x13, x26, x28
    umulh x14, x26, x28

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a8 * b10
    mul x13, x25, x3
    umulh x14, x25, x3

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  a7 * b11
    mul x13, x22, x4
    umulh x14, x22, x4

    adds x16, x16, x13
    adcs x15, x15, x14
    adcs x18, x18, xzr

    //  c18 is now in x16

    //  a8 * b11
    mul x13, x25, x4
    umulh x14, x25, x4

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, xzr, xzr

    //  a9 * b10
    mul x13, x26, x3
    umulh x14, x26, x3

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a10 * b9
    mul x13, x5, x28
    umulh x14, x5, x28

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  a11 * b8
    mul x13, x6, x27
    umulh x14, x6, x27

    adds x15, x15, x13
    adcs x18, x18, x14
    adcs x17, x17, xzr

    //  store c18 and c19
    stp x16, x15, [x2, #144]

    //  a11 * b9
    mul x13, x6, x28
    umulh x14, x6, x28

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, xzr, xzr

    //  a10 * b10
    mul x13, x5, x3
    umulh x14, x5, x3

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  a9 * b11
    mul x13, x26, x4
    umulh x14, x26, x4

    adds x18, x18, x13
    adcs x17, x17, x14
    adcs x16, x16, xzr

    //  c20 is now in x18

    //  a10 * b11
    mul x13, x5, x4
    umulh x14, x5, x4

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, xzr, xzr

    //  a11 * b10
    mul x13, x6, x3
    umulh x14, x6, x3

    adds x17, x17, x13
    adcs x16, x16, x14
    adcs x15, x15, xzr

    //  store c20 and c21
    stp x18, x17, [x2, #160]

    //  a11 * b11
    mul x13, x4, x6
    umulh x14, x4, x6

    adds x16, x16, x13
    adcs x15, x15, x14

    //  store c22 and c23
    stp x16, x15, [x2, #176]

    ldp x19, x20, [sp]
    ldp x21, x22, [sp, #16]
    ldp x23, x24, [sp, #32]
    ldp x25, x26, [sp, #48]
    ldp x27, x28, [sp, #64]
    add sp, sp, #80

    ret


//  three arguments:
//  array a in x0 which contains a 192-bit integer
//  array b in x1 which contains a 192-bit integer
//  array c in x2 which is the 384-bit result
fpmul3by3_arm64_asm:

    ldp x9, x10, [x0, #0]
    ldr x11, [x0, #16]

    ldp x12, x13, [x1, #0]
    ldr x14, [x1, #16]

    //  c0 computation

    mul x3, x9, x12
    umulh x4, x9, x12

    //  c1 computation

    mul x15, x9, x13
    umulh x5, x9, x13
    adds x4, x4, x15
    adcs x5, x5, xzr

    mul x15, x10, x12
    umulh x16, x10, x12
    adds x4, x4, x15
    adcs x5, x5, x16
    adcs x6, xzr, xzr

    //  c2 computation

    mul x15, x9, x14
    umulh x16, x9, x14
    adds x5, x5, x15
    adcs x6, x6, x16
    adcs x7, xzr, xzr

    mul x15, x10, x13
    umulh x16, x10, x13
    adds x5, x5, x15
    adcs x6, x6, x16
    adcs x7, x7, xzr

    mul x15, x11, x12
    umulh x16, x11, x12
    adds x5, x5, x15
    adcs x6, x6, x16
    adcs x7, x7, xzr

    //  c3 computation

    mul x15, x10, x14
    umulh x16, x10, x14
    adds x6, x6, x15
    adcs x7, x7, x16
    adcs x8, xzr, xzr

    mul x15, x11, x13
    umulh x16, x11, x13
    adds x6, x6, x15
    adcs x7, x7, x16
    adcs x8, x8, xzr

    //  c4 computation

    mul x15, x11, x14
    umulh x16, x11, x14

    adds x7, x7, x15
    adcs x8, x8, x16

    stp x3, x4, [x2, #0]
    stp x5, x6, [x2, #16]
    stp x7, x8, [x2, #32]

    ret


rdc751_arm64_asm:

    //  ma is in x0
    //  mc is in x1

    sub sp, sp, #80
    stp x19, x20, [sp]
    stp x21, x22, [sp, #16]
    stp x23, x24, [sp, #32]
    stp x25, x26, [sp, #48]
    stp x27, x28, [sp, #64]

    //  load the prime values into x22 through x28
    ldr x14, p751p1 + 0
    ldr x15, p751p1 + 8
    ldr x16, p751p1 + 16
    ldr x17, p751p1 + 24
    ldr x18, p751p1 + 32
    ldr x19, p751p1 + 40
    ldr x20, p751p1 + 48

    //  the values mc[0] through mc[11] will be held in x2 through x13
    //  until the very end when they will be stored

    //  load mc[0] through mc[4] and ma[5]
    ldp x2, x3, [x0, #0]
    ldp x4, x5, [x0, #16]
    ldp x6, x21, [x0, #32]

    //  ma[5] iteration
    mul x22, x2, x14   
    umulh x23, x2, x14 
    adds x24, x22, x21
    adcs x25, x23, xzr

    add x7, x24, xzr    //  set mc[5]

    //  ma[6] iteration

    ldr x21, [x0, #48]  //  ma[6]

    mul x22, x2, x15
    umulh x23, x2, x15
    adds x25, x25, x22
    adcs x26, x23, xzr

    mul x22, x3, x14
    umulh x23, x3, x14
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, xzr, xzr

    adds x25, x25, x21
    adcs x26, x26, xzr
    adcs x24, x24, xzr

    add x8, x25, xzr    //  set mc[6]

    //  ma[7] iteration

    ldr x21, [x0, #56]

    mul x22, x2, x16
    umulh x23, x2, x16
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, xzr, xzr

    mul x22, x3, x15
    umulh x23, x3, x15
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x4, x14
    umulh x23, x4, x14
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    adds x26, x26, x21
    adcs x24, x24, xzr
    adcs x25, x25, xzr

    add x9, x26, xzr    //  set mc[7] 

    //  ma[8] iteration

    ldr x21, [x0, #64]

    mul x22, x2, x17
    umulh x23, x2, x17
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, xzr, xzr

    mul x22, x3, x16
    umulh x23, x3, x16
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x4, x15
    umulh x23, x4, x15
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x5, x14
    umulh x23, x5, x14
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    adds x24, x24, x21
    adcs x25, x25, xzr
    adcs x26, x26, xzr

    add x10, x24, xzr   //  set mc[8]

    //  ma[9] iteration

    ldr x21, [x0, #72]

    mul x22, x2, x18
    umulh x23, x2, x18
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, xzr, xzr

    mul x22, x3, x17
    umulh x23, x3, x17
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x4, x16
    umulh x23, x4, x16
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x5, x15
    umulh x23, x5, x15
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x6, x14
    umulh x23, x6, x14
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    adds x25, x25, x21
    adcs x26, x26, xzr
    adcs x24, x24, xzr

    add x11, x25, xzr   //  set mc[9]

    //  ma[10] iteration

    ldr x21, [x0, #80]

    mul x22, x2, x19
    umulh x23, x2, x19
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, xzr, xzr

    mul x22, x3, x18
    umulh x23, x3, x18
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x4, x17
    umulh x23, x4, x17
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x5, x16
    umulh x23, x5, x16
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x6, x15
    umulh x23, x6, x15
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x7, x14
    umulh x23, x7, x14
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    adds x26, x26, x21
    adcs x24, x24, xzr
    adcs x25, x25, xzr

    add x12, x26, xzr   //  set mc[10]

    //  ma[11] iteration

    ldr x21, [x0, #88]

    mul x22, x2, x20
    umulh x23, x2, x20
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, xzr, xzr

    mul x22, x3, x19
    umulh x23, x3, x19
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x4, x18
    umulh x23, x4, x18
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x5, x17
    umulh x23, x5, x17
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x6, x16
    umulh x23, x6, x16
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x7, x15
    umulh x23, x7, x15
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x8, x14
    umulh x23, x8, x14
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    adds x24, x24, x21
    adcs x25, x25, xzr
    adcs x26, x26, xzr

    add x13, x24, xzr   //  set mc[11]

    //  ma[12] iteration

    ldr x21, [x0, #96]

    mul x22, x3, x20
    umulh x23, x3, x20
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, xzr, xzr

    mul x22, x4, x19
    umulh x23, x4, x19
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x5, x18
    umulh x23, x5, x18
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x6, x17
    umulh x23, x6, x17
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x7, x16
    umulh x23, x7, x16
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x8, x15
    umulh x23, x8, x15
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x9, x14
    umulh x23, x9, x14
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    adds x25, x25, x21
    adcs x26, x26, xzr
    adcs x24, x24, xzr

    add x2, x25, xzr   //  set mc[0]

    //  ma[13] iteration

    ldr x21, [x0, #104]

    mul x22, x4, x20
    umulh x23, x4, x20
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, xzr, xzr

    mul x22, x5, x19
    umulh x23, x5, x19
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x6, x18
    umulh x23, x6, x18
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x7, x17
    umulh x23, x7, x17
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x8, x16
    umulh x23, x8, x16
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x9, x15
    umulh x23, x9, x15
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x10, x14
    umulh x23, x10, x14
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    adds x26, x26, x21
    adcs x24, x24, xzr
    adcs x25, x25, xzr

    add x3, x26, xzr   //  set mc[1]

    //  ma[14] iteration

    ldr x21, [x0, #112]

    mul x22, x5, x20
    umulh x23, x5, x20
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, xzr, xzr

    mul x22, x6, x19
    umulh x23, x6, x19
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x7, x18
    umulh x23, x7, x18
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x8, x17
    umulh x23, x8, x17
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x9, x16
    umulh x23, x9, x16
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x10, x15
    umulh x23, x10, x15
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x11, x14
    umulh x23, x11, x14
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    adds x24, x24, x21
    adcs x25, x25, xzr
    adcs x26, x26, xzr

    add x4, x24, xzr   //  set mc[2]

    //  ma[15] iteration

    ldr x21, [x0, #120]

    mul x22, x6, x20
    umulh x23, x6, x20
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, xzr, xzr

    mul x22, x7, x19
    umulh x23, x7, x19
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x8, x18
    umulh x23, x8, x18
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x9, x17
    umulh x23, x9, x17
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x10, x16
    umulh x23, x10, x16
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x11, x15
    umulh x23, x11, x15
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x12, x14
    umulh x23, x12, x14
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    adds x25, x25, x21
    adcs x26, x26, xzr
    adcs x24, x24, xzr

    add x5, x25, xzr   //  set mc[3]

    //  ma[16] iteration

    ldr x21, [x0, #128]

    mul x22, x7, x20
    umulh x23, x7, x20
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, xzr, xzr

    mul x22, x8, x19
    umulh x23, x8, x19
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x9, x18
    umulh x23, x9, x18
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x10, x17
    umulh x23, x10, x17
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x11, x16
    umulh x23, x11, x16
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x12, x15
    umulh x23, x12, x15
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x13, x14
    umulh x23, x13, x14
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    adds x26, x26, x21
    adcs x24, x24, xzr
    adcs x25, x25, xzr

    add x6, x26, xzr   //  set mc[4]

    //  ma[17] iteration

    ldr x21, [x0, #136]

    mul x22, x8, x20
    umulh x23, x8, x20
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, xzr, xzr

    mul x22, x9, x19
    umulh x23, x9, x19
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x10, x18
    umulh x23, x10, x18
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x11, x17
    umulh x23, x11, x17
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x12, x16
    umulh x23, x12, x16
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x13, x15
    umulh x23, x13, x15
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    adds x24, x24, x21
    adcs x25, x25, xzr
    adcs x26, x26, xzr

    add x7, x24, xzr   //  set mc[5]

    //  ma[18] iteration

    ldr x21, [x0, #144]

    mul x22, x9, x20
    umulh x23, x9, x20
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, xzr, xzr

    mul x22, x10, x19
    umulh x23, x10, x19
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x11, x18
    umulh x23, x11, x18
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x12, x17
    umulh x23, x12, x17
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    mul x22, x13, x16
    umulh x23, x13, x16
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    adds x25, x25, x21
    adcs x26, x26, xzr
    adcs x24, x24, xzr

    add x8, x25, xzr   //  set mc[6]

    //  ma[19] iteration

    ldr x21, [x0, #152]

    mul x22, x10, x20
    umulh x23, x10, x20
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, xzr, xzr

    mul x22, x11, x19
    umulh x23, x11, x19
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x12, x18
    umulh x23, x12, x18
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    mul x22, x13, x17
    umulh x23, x13, x17
    adds x26, x26, x22
    adcs x24, x24, x23
    adcs x25, x25, xzr

    adds x26, x26, x21
    adcs x24, x24, xzr
    adcs x25, x25, xzr

    add x9, x26, xzr   //  set mc[7]

    //  ma[20] iteration

    ldr x21, [x0, #160]

    mul x22, x11, x20
    umulh x23, x11, x20
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, xzr, xzr

    mul x22, x12, x19
    umulh x23, x12, x19
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    mul x22, x13, x18
    umulh x23, x13, x18
    adds x24, x24, x22
    adcs x25, x25, x23
    adcs x26, x26, xzr

    adds x24, x24, x21
    adcs x25, x25, xzr
    adcs x26, x26, xzr

    add x10, x24, xzr   //  set mc[8]

    //  ma[21] iteration

    ldr x21, [x0, #168]

    mul x22, x12, x20
    umulh x23, x12, x20
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, xzr, xzr

    mul x22, x13, x19
    umulh x23, x13, x19
    adds x25, x25, x22
    adcs x26, x26, x23
    adcs x24, x24, xzr

    adds x25, x25, x21
    adcs x26, x26, xzr
    adcs x24, x24, xzr

    add x11, x25, xzr   //  set mc[9]

    //  ma[22] iteration

    ldr x21, [x0, #176]

    mul x22, x13, x20
    umulh x23, x13, x20
    adds x26, x26, x22
    adcs x24, x24, x23

    adds x26, x26, x21

    ldr x21, [x0, #184]

    adcs x24, x24, x21

    add x12, x26, xzr   //  set mc[10]
    add x13, x24, xzr   //  set mc[11]


    //  add the negative of p751 to the resut in x2 - x13 (i.e, subtract p751)

    ldr x14, neg_p751
    ldr x15, neg_p751 + 8
    adds x2, x2, x14
    adcs x3, x3, x15
    adcs x4, x4, x15
    adcs x5, x5, x15
    adcs x6, x6, x15
    ldr x16, neg_p751 + 16
    ldr x17, neg_p751 + 24
    adcs x7, x7, x16
    ldr x14, neg_p751 + 32
    adcs x8, x8, x17
    ldr x15, neg_p751 + 40
    adcs x9, x9, x14
    ldr x16, neg_p751 + 48
    adcs x10, x10, x15
    ldr x17, neg_p751 + 56
    adcs x11, x11, x16
    ldr x14, neg_p751 + 64
    adcs x12, x12, x17
    adcs x13, x13, x14

    sbc x14, xzr, xzr


    //  add p751 back but anded with the mask in x14

    ldr x15, p751
    and x15, x15, x14
    ldr x16, p751 + 8
    and x16, x16, x14
    ldr x17, p751 + 16
    and x17, x17, x14    

    adds x2, x2, x15
    adcs x3, x3, x15
    adcs x4, x4, x15
    adcs x5, x5, x15
    adcs x6, x6, x15
    adcs x7, x7, x16
    adcs x8, x8, x17

    ldr x15, p751 + 24
    and x15, x15, x14
    ldr x16, p751 + 32
    and x16, x16, x14
    ldr x17, p751 + 40
    and x17, x17, x14 

    adcs x9, x9, x15   
    adcs x10, x10, x16   
    adcs x11, x11, x17   

    ldr x15, p751 + 48
    and x15, x15, x14
    ldr x16, p751 + 56
    and x16, x16, x14

    adcs x12, x12, x15
    adcs x13, x13, x16


    stp x2, x3, [x1, #0]
    stp x4, x5, [x1, #16]
    stp x6, x7, [x1, #32]
    stp x8, x9, [x1, #48]
    stp x10, x11, [x1, #64]
    stp x12, x13, [x1, #80]

    
    ldp x19, x20, [sp]
    ldp x21, x22, [sp, #16]
    ldp x23, x24, [sp, #32]
    ldp x25, x26, [sp, #48]
    ldp x27, x28, [sp, #64]
    add sp, sp, #80

    ret    
