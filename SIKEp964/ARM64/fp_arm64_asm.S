/**********************************************************************************
*
*	Abstract: 964-bit addition/subtraction and Comba-based Montgomery reduction 
*	The motivation behind developing this library is to use inside SIDH key-exchange
*	protocol at 160-bit quantum security level
*
*	Author: Amir Jalali					Date: November 2016
*
*	                Florida Atlantic University
*	                    All rights reserved
***********************************************************************************/

.text
.global fpadd964_arm64_asm
.global fpsub964_arm64_asm
.global rdc964_arm64_asm
.global sub1024
.global fpadd1024_arm64
.global fpsub1024_arm64

.global mp_add964_asm
.global mp_add964x2_asm
.global mp_sub964x2_asm

mp_add964_asm:
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
    ldp     x7, x8,   [x0, #96]
    ldp     x9, x10,  [x0, #112]
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
    adc     x10, x10, x18
    stp     x3,  x4,  [x2, #64]
    stp     x5,  x6,  [x2, #80]
    stp     x7,  x8,  [x2, #96]
    stp     x9,  x10,  [x2, #112]

    ret

mp_add964x2_asm:
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
    ldp     x7, x8,   [x0, #96]
    ldp     x9, x10,  [x0, #112]
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
    stp     x9,  x10,  [x2, #112]

    ldp     x3, x4,   [x0, #128]
    ldp     x5, x6,   [x0, #144]
    ldp     x7, x8,   [x0, #160]
    ldp     x9, x10,  [x0, #176]
    // load register of b 
    ldp     x11, x12, [x1, #128]
    ldp     x13, x14, [x1, #144]
    ldp     x15, x16, [x1, #160]
    ldp     x17, x18, [x1, #176]
    adds    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adcs    x6,  x6,  x14
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adcs    x10, x10, x18
    stp     x3,  x4,  [x2, #128]
    stp     x5,  x6,  [x2, #144]
    stp     x7,  x8,  [x2, #160]
    stp     x9,  x10, [x2, #176]
    // load register of a 
    ldp     x3, x4,   [x0, #192]
    ldp     x5, x6,   [x0, #208]
    ldp     x7, x8,   [x0, #224]
    ldp     x9, x10,  [x0, #240]
    // load register of b 
    ldp     x11, x12, [x1, #192]
    ldp     x13, x14, [x1, #208]
    ldp     x15, x16, [x1, #224]
    ldp     x17, x18, [x1, #240]
    adcs    x3,  x3,  x11
    adcs    x4,  x4,  x12
    adcs    x5,  x5,  x13
    adcs    x6,  x6,  x14
    adcs    x7,  x7,  x15
    adcs    x8,  x8,  x16
    adcs    x9,  x9,  x17
    adc     x10, x10, x18
    stp     x3,  x4,  [x2, #192]
    stp     x5,  x6,  [x2, #208]
    stp     x7,  x8,  [x2, #224]
    stp     x9,  x10,  [x2, #240]

    ret

mp_sub964x2_asm:
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
    stp     x3,  x4,  [x2, #128]
    stp     x5,  x6,  [x2, #144]
    stp     x7,  x8,  [x2, #160]
    stp     x9,  x10, [x2, #176]

    ldp     x3, x4,   [x0, #192]
    ldp     x5, x6,   [x0, #208]
    ldp     x7, x8,   [x0, #224]
    ldp     x9, x10,  [x0, #240]
    // load register of b 
    ldp     x11, x12, [x1, #192]
    ldp     x13, x14, [x1, #208]
    ldp     x15, x16, [x1, #224]
    ldp     x17, x18, [x1, #240]
    sbcs    x3,  x3,  x11
    sbcs    x4,  x4,  x12
    sbcs    x5,  x5,  x13
    sbcs    x6,  x6,  x14
    sbcs    x7,  x7,  x15
    sbcs    x8,  x8,  x16
    sbcs    x9,  x9,  x17
    sbcs     x10, x10, x18
    mov     x0,  #0
    sbc     x0,  x0,  x0
    stp     x3,  x4,  [x2, #192]
    stp     x5,  x6,  [x2, #208]
    stp     x7,  x8,  [x2, #224]
    stp     x9,  x10,  [x2, #240]

    ret    



fpadd1024_arm64:

    sub sp, sp, #112
    stp x15, x16, [sp]
    stp x17, x18, [sp, #16]
    stp x19, x20, [sp, #32]
    stp x21, x22, [sp, #48]
    stp x23, x24, [sp, #64]
    stp x25, x26, [sp, #80]
    stp x27, x28, [sp, #96]
	ldp x3, x4, [x0]
	ldp x5, x6, [x0, #16]
	ldp x7, x8, [x0, #32]
	ldp x9, x10, [x0, #48]
	ldp x11, x12, [x1]
	ldp x13, x14, [x1, #16]
	ldp x15, x16, [x1, #32]
	ldp x17, x18, [x1, #48]
	adds x3, x3, x11
	adcs x4, x4, x12
	adcs x5, x5, x13
	adcs x6, x6, x14
	adcs x7, x7, x15
	adcs x8, x8, x16
	adcs x9, x9, x17
	adcs x10, x10, x18
	stp x3, x4, [x2]
	stp x5, x6, [x2, #16]
	stp x7, x8, [x2, #32]
	stp x9, x10, [x2, #48]
	ldp x3, x4, [x0, #64]
	ldp x5, x6, [x0, #80]
	ldp x7, x8, [x0, #96]
	ldp x9, x10, [x0, #112]
	ldp x11, x12, [x1, #64]
	ldp x13, x14, [x1, #80]
	ldp x15, x16, [x1, #96]
	ldp x17, x18, [x1, #112]
	adcs x3, x3, x11
	adcs x4, x4, x12
	adcs x5, x5, x13
	adcs x6, x6, x14
	adcs x7, x7, x15
	adcs x8, x8, x16
	adcs x9, x9, x17
	adcs x10, x10, x18
	stp x3, x4, [x2, #64]
	stp x5, x6, [x2, #80]
	stp x7, x8, [x2, #96]
	stp x9, x10, [x2, #112]
    ldp x15, x16, [sp]
    ldp x17, x18, [sp, #16]
    ldp x19, x20, [sp, #32]
    ldp x21, x22, [sp, #48]
    ldp x23, x24, [sp, #64]
    ldp x25, x26, [sp, #80]
    ldp x27, x28, [sp, #96]
    add sp, sp, #112
ret

fpsub1024_arm64:

    sub sp, sp, #112
    stp x15, x16, [sp]
    stp x17, x18, [sp, #16]
    stp x19, x20, [sp, #32]
    stp x21, x22, [sp, #48]
    stp x23, x24, [sp, #64]
    stp x25, x26, [sp, #80]
    stp x27, x28, [sp, #96]
	ldp x3, x4, [x0]
	ldp x5, x6, [x0, #16]
	ldp x7, x8, [x0, #32]
	ldp x9, x10, [x0, #48]
	ldp x11, x12, [x1]
	ldp x13, x14, [x1, #16]
	ldp x15, x16, [x1, #32]
	ldp x17, x18, [x1, #48]
	subs x3, x3, x11
	sbcs x4, x4, x12
	sbcs x5, x5, x13
	sbcs x6, x6, x14
	sbcs x7, x7, x15
	sbcs x8, x8, x16
	sbcs x9, x9, x17
	sbcs x10, x10, x18
	stp x3, x4, [x2]
	stp x5, x6, [x2, #16]
	stp x7, x8, [x2, #32]
	stp x9, x10, [x2, #48]
	ldp x3, x4, [x0, #64]
	ldp x5, x6, [x0, #80]
	ldp x7, x8, [x0, #96]
	ldp x9, x10, [x0, #112]
	ldp x11, x12, [x1, #64]
	ldp x13, x14, [x1, #80]
	ldp x15, x16, [x1, #96]
	ldp x17, x18, [x1, #112]
	sbcs x3, x3, x11
	sbcs x4, x4, x12
	sbcs x5, x5, x13
	sbcs x6, x6, x14
	sbcs x7, x7, x15
	sbcs x8, x8, x16
	sbcs x9, x9, x17
	sbcs x10, x10, x18
	stp x3, x4, [x2, #64]
	stp x5, x6, [x2, #80]
	stp x7, x8, [x2, #96]
	stp x9, x10, [x2, #112]
    ldp x15, x16, [sp]
    ldp x17, x18, [sp, #16]
    ldp x19, x20, [sp, #32]
    ldp x21, x22, [sp, #48]
    ldp x23, x24, [sp, #64]
    ldp x25, x26, [sp, #80]
    ldp x27, x28, [sp, #96]
    add sp, sp, #112
ret

fpadd964_arm64_asm:

    sub sp, sp, #112
    stp x15, x16, [sp]
    stp x17, x18, [sp, #16]
    stp x19, x20, [sp, #32]
    stp x21, x22, [sp, #48]
    stp x23, x24, [sp, #64]
    stp x25, x26, [sp, #80]
    stp x27, x28, [sp, #96]
    ldp x3, x4,   [x0,#0]
    ldp x5, x6,   [x0,#16]
    ldp x7, x8,   [x0,#32]
    ldp x9, x10,  [x0,#48]
    ldp x11, x12, [x0,#64]
    ldp x13, x14, [x0,#80]
    ldp x15, x16, [x0,#96]
    ldp x17, x18, [x0,#112]
    ldp x19, x20,   [x1,#0]
    ldp x21, x22,   [x1,#16]
    adds x3, x3, x19
    adcs x4, x4, x20
    adcs x5, x5, x21
    adcs x6, x6, x22
    ldp x19, x20,   [x1,#32]
    ldp x21, x22,   [x1,#48]
    adcs x7, x7, x19
    adcs x8, x8, x20
    adcs x9, x9, x21
    adcs x10, x10, x22
    ldp x19, x20,   [x1,#64]
    ldp x21, x22,   [x1,#80]
    adcs x11, x11, x19
    adcs x12, x12, x20
    adcs x13, x13, x21
    adcs x14, x14, x22
    ldp x19, x20,   [x1,#96]
    ldp x21, x22,   [x1,#112]
    adcs x15, x15, x19
    adcs x16, x16, x20
    adcs x17, x17, x21
    adcs x18, x18, x22
    ldr x19, =0x0000000000000001
    ldr x20, =0x0000000000000000
    adds x3, x3, x19
    adcs x4, x4, x20
    adcs x5, x5, x20
    adcs x6, x6, x20
    adcs x7, x7, x20
    adcs x8, x8, x20
    adcs x9, x9, x20
    ldr x21, =0xbae32b4000000000
    ldr x22, =0x544c7154b985321a
    adcs x10, x10, x21
    ldr x19, =0x1a910955a856b8b6
    adcs x11, x11, x22
    ldr x20, =0xf6c0d47252ad7e18
    adcs x12, x12, x19
    ldr x21, =0x3454dbaecab96454
    adcs x13, x13, x20
    ldr x22, =0xaf34558a5d5e05bb
    adcs x14, x14, x21
    ldr x19, =0xeffd7db752b03b4e
    adcs x15, x15, x22
    ldr x20, =0x94a40089bc39b085
    adcs x16, x16, x19
    ldr x21, =0xfffffffffffffff7
    adcs x17, x17, x20
    adcs x18, x18, x21 
    sbc x19, xzr, xzr
    ldr x20, p964
    and x20, x20, x19
    ldr x21, =0x451cd4bfffffffff
    and x21, x21, x19
    ldr x22, =0xabb38eab467acde5
    and x22, x22, x19    
    adds x3, x3, x20
    adcs x4, x4, x20
    adcs x5, x5, x20
    adcs x6, x6, x20
    adcs x7, x7, x20
    adcs x8, x8, x20
    adcs x9, x9, x20
    adcs x10, x10, x21
    adcs x11, x11, x22
    ldr x23,=0xe56ef6aa57a94749
    and x23, x23, x19
    ldr x24, =0x093f2b8dad5281e7
    and x24, x24, x19
    ldr x25, =0xcbab245135469bab
    and x25, x25, x19 
    adcs x12, x12, x23   
    adcs x13, x13, x24   
    adcs x14, x14, x25   
    ldr x26, =0x50cbaa75a2a1fa44
    and x26, x26, x19
    ldr x27, =0x10028248ad4fc4b1
    and x27, x27, x19
    adcs x15, x15, x26
    adcs x16, x16, x27
    ldr x28, =0x6b5bff7643c64f7a
    and x28, x28, x19
    ldr x20, =0x0000000000000008
    and x20, x20, x19
    adcs x17, x17, x28
    adcs x18, x18, x20
    stp x3, x4,   [x2,#0]
    stp x5, x6,   [x2,#16]
    stp x7, x8,   [x2,#32]
    stp x9, x10,  [x2,#48]
    stp x11, x12, [x2,#64]
    stp x13, x14, [x2,#80]
    stp x15, x16, [x2,#96]
    stp x17, x18, [x2,#112]
    ldp x15, x16, [sp]
    ldp x17, x18, [sp, #16]
    ldp x19, x20, [sp, #32]
    ldp x21, x22, [sp, #48]
    ldp x23, x24, [sp, #64]
    ldp x25, x26, [sp, #80]
    ldp x27, x28, [sp, #96]
    add sp, sp, #112
ret

fpsub964_arm64_asm:
    sub sp, sp, #112
    stp x15, x16, [sp]
    stp x17, x18, [sp, #16]
    stp x19, x20, [sp, #32]
    stp x21, x22, [sp, #48]
    stp x23, x24, [sp, #64]
    stp x25, x26, [sp, #80]
    stp x27, x28, [sp, #96]
    ldp x3, x4,   [x0,#0]
    ldp x5, x6,   [x0,#16]
    ldp x7, x8,   [x0,#32]
    ldp x9, x10,  [x0,#48]
    ldp x11, x12, [x0,#64]
    ldp x13, x14, [x0,#80]
    ldp x15, x16, [x0,#96]
    ldp x17, x18, [x0,#112]
    ldp x19, x20, [x1, #0]
    ldp x21, x22, [x1, #16]
    subs x3, x3, x19
    sbcs x4, x4, x20
    sbcs x5, x5, x21
    sbcs x6, x6, x22
    ldp x19, x20, [x1, #32]
    ldp x21, x22, [x1, #48]
    sbcs x7, x7, x19
    sbcs x8, x8, x20
    sbcs x9, x9, x21
    sbcs x10, x10, x22
    ldp x19, x20, [x1, #64]
    ldp x21, x22, [x1, #80]
    sbcs x11, x11, x19
    sbcs x12, x12, x20
    sbcs x13, x13, x21
    sbcs x14, x14, x22
    ldp x19, x20, [x1, #96]
    ldp x21, x22, [x1, #112]
    sbcs x15, x15, x19
    sbcs x16, x16, x20
    sbcs x17, x17, x21
    sbcs x18, x18, x22
    sbc x19, xzr, xzr
    ldr x20, =0xFFFFFFFFFFFFFFFF
    and x20, x20, x19
    ldr x21, =0x451cd4bfffffffff
    and x21, x21, x19
    adds x3, x3, x20
    adcs x4, x4, x20
    adcs x5, x5, x20
    adcs x6, x6, x20
    adcs x7, x7, x20
    adcs x8, x8, x20
    adcs x9, x9, x20
    adcs x10, x10, x21
    ldr x22, =0xabb38eab467acde5
    and x22, x22, x19    
    ldr x23,=0xe56ef6aa57a94749
    and x23, x23, x19
    adcs x11, x11, x22
    adcs x12, x12, x23   
    ldr x20, =0x093f2b8dad5281e7
    and x20, x20, x19
    ldr x21, =0xcbab245135469bab
    and x21, x21, x19 
    adcs x13, x13, x20   
    adcs x14, x14, x21   
    ldr x22, =0x50cbaa75a2a1fa44
    and x22, x22, x19
    ldr x23, =0x10028248ad4fc4b1
    and x23, x23, x19
    adcs x15, x15, x22
    adcs x16, x16, x23
    ldr x20, =0x6b5bff7643c64f7a
    and x20, x20, x19
    ldr x21, =0x0000000000000008
    and x21, x21, x19
    adcs x17, x17, x20
    adc x18, x18, x21
    stp x3, x4,   [x2,#0]
    stp x5, x6,   [x2,#16]
    stp x7, x8,   [x2,#32]
    stp x9, x10,  [x2,#48]
    stp x11, x12, [x2,#64]
    stp x13, x14, [x2,#80]
    stp x15, x16, [x2,#96]
    stp x17, x18, [x2,#112]
    ldp x15, x16, [sp]
    ldp x17, x18, [sp, #16]
    ldp x19, x20, [sp, #32]
    ldp x21, x22, [sp, #48]
    ldp x23, x24, [sp, #64]
    ldp x25, x26, [sp, #80]
    ldp x27, x28, [sp, #96]
    add sp, sp, #112

ret

rdc964_arm64_asm:
    sub sp, sp, #112
    stp x15, x16, [sp]
    stp x17, x18, [sp, #16]
    stp x19, x20, [sp, #32]
    stp x21, x22, [sp, #48]
    stp x23, x24, [sp, #64]
    stp x25, x26, [sp, #80]
    stp x27, x28, [sp, #96]
    ldr x18, =0x451cd4c000000000
    ldr x19, =0xabb38eab467acde5
    ldr x20, =0xe56ef6aa57a94749
    ldr x21, =0x093f2b8dad5281e7
    ldr x22, =0xcbab245135469bab
    ldp x2, x3, [x0, #0]
    ldp x4, x5, [x0, #16]
    ldp x6, x7, [x0, #32]
    ldp x8, x26, [x0, #48]
    mul x27, x2, x18   	
    umulh x28, x2, x18 
    adds x27, x27, x26
    adcs x23, x28, xzr
    add x9, x27, xzr// x9 -> mc[7]
    ldr x26, [x0, #64]  //  ma[8]
    mul x27, x2, x19
    umulh x28, x2, x19
    adds x23, x23, x27
    adcs x25, x28, xzr
    mul x27, x3, x18
    umulh x28, x3, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x10, x23, xzr  //  x10 -> mc[8]
    ldr x26, [x0, #72]
    mul x27, x2, x20
    umulh x28, x2, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, xzr, xzr
    mul x27, x3, x19
    umulh x28, x3, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x4, x18
    umulh x28, x4, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    adds x25, x25, x26
    adcs x24, x24, xzr
    adcs x23, x23, xzr
    add x11, x25, xzr  //  x11 -> mc[9] 
    ldr x26, [x0, #80]
    mul x27, x2, x21
    umulh x28, x2, x21
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, xzr, xzr
    mul x27, x3, x20
    umulh x28, x3, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x4, x19
    umulh x28, x4, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x5, x18
    umulh x28, x5, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    adds x24, x24, x26
    adcs x23, x23, xzr
    adcs x25, x25, xzr
    add x12, x24, xzr  //  x12 ->  mc[10]
    ldr x26, [x0, #88]
    mul x27, x2, x22
    umulh x28, x2, x22
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    mul x27, x3, x21
    umulh x28, x3, x21
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x4, x20
    umulh x28, x4, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x5, x19
    umulh x28, x5, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x6, x18
    umulh x28, x6, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x13, x23, xzr  //  x13 -> set mc[11]
    ldr x26, [x0, #96]
    ldr x18, =0x50cbaa75a2a1fa44
    mul x27, x2, x18
    umulh x28, x2, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, xzr, xzr
    mul x27, x3, x22
    umulh x28, x3, x22
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x4, x21
    umulh x28, x4, x21
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x5, x20
    umulh x28, x5, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x6, x19
    umulh x28, x6, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x18, =0x451cd4c000000000
    mul x27, x7, x18
    umulh x28, x7, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    adds x25, x25, x26
    adcs x24, x24, xzr
    adcs x23, x23, xzr
    add x14, x25, xzr  //  x14 -> mc[12]
    ldr x26, [x0, #104]
    ldr x18, =0x10028248ad4fc4b1;
    mul x27, x2, x18
    umulh x28, x2, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, xzr, xzr
    ldr x18, =0x50cbaa75a2a1fa44
    mul x27, x3, x18
    umulh x28, x3, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x4, x22
    umulh x28, x4, x22
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x5, x21
    umulh x28, x5, x21
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x6, x20
    umulh x28, x6, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x7, x19
    umulh x28, x7, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x18, =0x451cd4c000000000
    mul x27, x8, x18
    umulh x28, x8, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    adds x24, x24, x26
    adcs x23, x23, xzr
    adcs x25, x25, xzr
    add x15, x24, xzr  //  x15 -> set mc[13]
    ldr x26, [x0, #112]
    mul x27, x9, x18
    umulh x28, x9, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    mul x27, x8, x19
    umulh x28, x8, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x7, x20
    umulh x28, x7, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x6, x21
    umulh x28, x6, x21
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x5, x22
    umulh x28, x5, x22
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x18, =0x50cbaa75a2a1fa44
    mul x27, x4, x18
    umulh x28, x4, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x18, =0x10028248ad4fc4b1
    mul x27, x3, x18
    umulh x28, x3, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x19, =0x6b5bff7643c64f7a
    mul x27, x2, x19
    umulh x28, x2, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x16, x23, xzr //  x16 ->  mc[14]
    ldr x26, [x0, #120]
    mul x27, x3, x19
    umulh x28, x3, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, xzr, xzr
    mul x27, x4, x18
    umulh x28, x4, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x18, =0x0000000000000008
    mul x27, x2, x18
    umulh x28, x2, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x19, =0x50cbaa75a2a1fa44
    mul x27, x5, x19
    umulh x28, x5, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x6, x22
    umulh x28, x6, x22
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x7, x21
    umulh x28, x7, x21
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x8, x20
    umulh x28, x8, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x19, =0xabb38eab467acde5
    mul x27, x9, x19
    umulh x28, x9, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x18, =0x451cd4c000000000
    mul x27, x10, x18
    umulh x28, x10, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    adds x25, x25, x26
    adcs x24, x24, xzr
    adcs x23, x23, xzr
    add x17, x25, xzr  // x17 -> set mc[15]
    ldr x26, [x0, #128]
    mul x27, x11, x18
    umulh x28, x11, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, xzr, xzr
    mul x27, x10, x19
    umulh x28, x10, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x9, x20
    umulh x28, x9, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x8, x21
    umulh x28, x8, x21
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x7, x22
    umulh x28, x7, x22
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x18, =0x50cbaa75a2a1fa44
    mul x27, x6, x18
    umulh x28, x6, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x18, =0x10028248ad4fc4b1
    mul x27, x5, x18
    umulh x28, x5, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x19, =0x6b5bff7643c64f7a
    mul x27, x4, x19
    umulh x28, x4, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x20, =0x0000000000000008
    mul x27, x3, x20
    umulh x28, x3, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    adds x24, x24, x26
    adcs x23, x23, xzr
    adcs x25, x25, xzr
    add x2, x24, xzr  //  x2 ->  mc[0]
    ldr x26, [x0, #136]
    mul x27, x4, x20
    umulh x28, x4, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    mul x27, x5, x19
    umulh x28, x5, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x6, x18
    umulh x28, x6, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x20, =0x50cbaa75a2a1fa44
    mul x27, x7, x20
    umulh x28, x7, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x8, x22
    umulh x28, x8, x22
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x9, x21
    umulh x28, x9, x21
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x20, =0xe56ef6aa57a94749
    mul x27, x10, x20
    umulh x28, x10, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x19, =0xabb38eab467acde5
    mul x27, x11, x19
    umulh x28, x11, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x18, =0x451cd4c000000000
    mul x27, x12, x18
    umulh x28, x12, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x3, x23, xzr  //  x3 ->  mc[1]
    ldr x26, [x0, #144]
    mul x27, x13, x18
    umulh x28, x13, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, xzr, xzr
    mul x27, x12, x19
    umulh x28, x12, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x11, x20
    umulh x28, x11, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x10, x21
    umulh x28, x10, x21
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x9, x22
    umulh x28, x9, x22
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x18, =0x50cbaa75a2a1fa44
    mul x27, x8, x18
    umulh x28, x8, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x18, =0x10028248ad4fc4b1
    mul x27, x7, x18
    umulh x28, x7, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x19, =0x6b5bff7643c64f7a
    mul x27, x6, x19
    umulh x28, x6, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x20, =0x0000000000000008
    mul x27, x5, x20
    umulh x28, x5, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    adds x25, x25, x26
    adcs x24, x24, xzr
    adcs x23, x23, xzr
    add x4, x25, xzr  // x4 ->  mc[2]
    ldr x26, [x0, #152]
    mul x27, x6, x20
    umulh x28, x6, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, xzr, xzr
    mul x27, x7, x19
    umulh x28, x7, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x8, x18
    umulh x28, x8, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x20, =0x50cbaa75a2a1fa44
    mul x27, x9, x20
    umulh x28, x9, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x10, x22
    umulh x28, x10, x22
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x11, x21
    umulh x28, x11, x21
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x20, =0xe56ef6aa57a94749
    mul x27, x12, x20
    umulh x28, x12, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x19, =0xabb38eab467acde5
    mul x27, x13, x19
    umulh x28, x13, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x18, =0x451cd4c000000000
    mul x27, x14, x18
    umulh x28, x14, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    adds x24, x24, x26
    adcs x23, x23, xzr
    adcs x25, x25, xzr
    add x5, x24, xzr	//  x5 ->  mc[3]
    ldr x26, [x0, #160]
    mul x27, x15, x18
    umulh x28, x15, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    mul x27, x14, x19
    umulh x28, x14, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x13, x20
    umulh x28, x13, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x12, x21
    umulh x28, x12, x21
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x11, x22
    umulh x28, x11, x22
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x18, =0x50cbaa75a2a1fa44
    mul x27, x10, x18
    umulh x28, x10, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x18, =0x10028248ad4fc4b1
    mul x27, x9, x18
    umulh x28, x9, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x19, =0x6b5bff7643c64f7a
    mul x27, x8, x19
    umulh x28, x8, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x20, =0x0000000000000008
    mul x27, x7, x20
    umulh x28, x7, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x6, x23, xzr   //  set mc[4]
    ldr x26, [x0, #168]
    mul x27, x8, x20
    umulh x28, x8, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, xzr, xzr
    mul x27, x9, x19
    umulh x28, x9, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x10, x18
    umulh x28, x10, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x20, =0x50cbaa75a2a1fa44
    mul x27, x11, x20
    umulh x28, x11, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x12, x22
    umulh x28, x12, x22
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x13, x21
    umulh x28, x13, x21
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x20, =0xe56ef6aa57a94749
    mul x27, x14, x20
    umulh x28, x14, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x19, =0xabb38eab467acde5
    mul x27, x15, x19
    umulh x28, x15, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x18, =0x451cd4c000000000
    mul x27, x16, x18
    umulh x28, x16, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    adds x25, x25, x26
    adcs x24, x24, xzr
    adcs x23, x23, xzr
    add x7, x25, xzr   //  set mc[5]
    ldr x26, [x0, #176]
    mul x27, x17, x18
    umulh x28, x17, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, xzr, xzr
    mul x27, x16, x19
    umulh x28, x16, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x15, x20
    umulh x28, x15, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x14, x21
    umulh x28, x14, x21
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x13, x22
    umulh x28, x13, x22
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x18, =0x50cbaa75a2a1fa44
    mul x27, x12, x18
    umulh x28, x12, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x18, =0x10028248ad4fc4b1
    mul x27, x11, x18
    umulh x28, x11, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x19, =0x6b5bff7643c64f7a
    mul x27, x10, x19
    umulh x28, x10, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x20, =0x0000000000000008
    mul x27, x9, x20
    umulh x28, x9, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    adds x24, x24, x26
    adcs x23, x23, xzr
    adcs x25, x25, xzr
    add x8, x24, xzr   //  set mc[6]
    ldr x26, [x0, #184]
    mul x27, x10, x20
    umulh x28, x10, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    mul x27, x11, x19
    umulh x28, x11, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x12, x18
    umulh x28, x12, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x20, =0x50cbaa75a2a1fa44
    mul x27, x13, x20
    umulh x28, x13, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x14, x22
    umulh x28, x14, x22
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x15, x21
    umulh x28, x15, x21
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x20, =0xe56ef6aa57a94749
    mul x27, x16, x20
    umulh x28, x16, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    ldr x19, =0xabb38eab467acde5
    mul x27, x17, x19
    umulh x28, x17, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x9, x23, xzr   //  set mc[7]
    ldr x26, [x0, #192]
    mul x27, x17, x20
    umulh x28, x17, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, xzr, xzr
    mul x27, x16, x21
    umulh x28, x16, x21
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x15, x22
    umulh x28, x15, x22
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x19, =0x50cbaa75a2a1fa44
    mul x27, x14, x19
    umulh x28, x14, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x13, x18
    umulh x28, x13, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x19, =0x6b5bff7643c64f7a
    mul x27, x12, x19
    umulh x28, x12, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    ldr x20, =0x0000000000000008
    mul x27, x11, x20
    umulh x28, x11, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    adds x25, x25, x26
    adcs x24, x24, xzr
    adcs x23, x23, xzr
    add x10, x25, xzr// set mc[8]
    ldr x26, [x0, #200]
    mul x27, x12, x20
    umulh x28, x12, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, xzr, xzr
    mul x27, x13, x19
    umulh x28, x13, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x14, x18
    umulh x28, x14, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x16, x22
    umulh x28, x16, x22
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x17, x21
    umulh x28, x17, x21
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    ldr x21, =0x50cbaa75a2a1fa44
    mul x27, x15, x21
    umulh x28, x15, x21
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    adds x24, x24, x26
    adcs x23, x23, xzr
    adcs x25, x25, xzr
    add x11, x24, xzr	//set mc[9]
    ldr x26, [x0, #208]
    mul x27, x13, x20
    umulh x28, x13, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    mul x27, x14, x19
    umulh x28, x14, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x15, x18
    umulh x28, x15, x18
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x16, x21
    umulh x28, x16, x21
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    mul x27, x17, x22
    umulh x28, x17, x22
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x12, x23, xzr	//set mc[10]
    ldr x26, [x0, #216]
    mul x27, x14, x20
    umulh x28, x14, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, xzr, xzr
    mul x27, x15, x19
    umulh x28, x15, x19
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x16, x18
    umulh x28, x16, x18
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    mul x27, x17, x21
    umulh x28, x17, x21
    adds x25, x25, x27
    adcs x24, x24, x28
    adcs x23, x23, xzr
    adds x25, x25, x26
    adcs x24, x24, xzr
    adcs x23, x23, xzr
    add x13, x25, xzr	//set mc[11]
    ldr x26, [x0, #224]
    mul x27, x15, x20
    umulh x28, x15, x20
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, xzr, xzr
    mul x27, x16, x19
    umulh x28, x16, x19
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    mul x27, x17, x18
    umulh x28, x17, x18
    adds x24, x24, x27
    adcs x23, x23, x28
    adcs x25, x25, xzr
    adds x24, x24, x26
    adcs x23, x23, xzr
    adcs x25, x25, xzr
    add x14, x24, xzr	//set mc[12]
    ldr x26, [x0, #232]
    mul x27, x16, x20
    umulh x28, x16, x20
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, xzr, xzr
    mul x27, x17, x19
    umulh x28, x17, x19
    adds x23, x23, x27
    adcs x25, x25, x28
    adcs x24, x24, xzr
    adds x23, x23, x26
    adcs x25, x25, xzr
    adcs x24, x24, xzr
    add x15, x23, xzr	//set mc[13]
    ldr x26, [x0, #240]
    mul x27, x17, x20
    umulh x28, x17, x20
    adds x25, x25, x27
    adcs x24, x24, x28
    adds x25, x25, x26
    ldr x26, [x0, #248]
    adcs x24, x24, x26
    add x16, x25, xzr // set mc[14]
    add x17, x24, xzr // set mc[15]
    ldr x19, =0x0000000000000001
    ldr x20, =0x0000000000000000
    adds x2, x2, x19
    adcs x3, x3, x20
    adcs x4, x4, x20
    adcs x5, x5, x20
    adcs x6, x6, x20
    adcs x7, x7, x20
    adcs x8, x8, x20
    ldr x21, =0xbae32b4000000000
    ldr x22, =0x544c7154b985321a
    adcs x9, x9, x21
    ldr x19, =0x1a910955a856b8b6
    adcs x10, x10, x22
    ldr x20, =0xf6c0d47252ad7e18
    adcs x11, x11, x19
    ldr x21, =0x3454dbaecab96454
    adcs x12, x12, x20
    ldr x22, =0xaf34558a5d5e05bb
    adcs x13, x13, x21
    ldr x19, =0xeffd7db752b03b4e
    adcs x14, x14, x22
    ldr x20, =0x94a40089bc39b085
    adcs x15, x15, x19
    ldr x21, =0xfffffffffffffff7
    adcs x16, x16, x20
    adcs x17, x17, x21 
    sbc x19, xzr, xzr
    ldr x20, p964
    and x20, x20, x19
    ldr x21, =0x451cd4bfffffffff
    and x21, x21, x19
    ldr x22, =0xabb38eab467acde5
    and x22, x22, x19    
    adds x2, x2, x20
    adcs x3, x3, x20
    adcs x4, x4, x20
    adcs x5, x5, x20
    adcs x6, x6, x20
    adcs x7, x7, x20
    adcs x8, x8, x20
    adcs x9, x9, x21
    adcs x10, x10, x22
    ldr x20,=0xe56ef6aa57a94749
    and x20, x20, x19
    ldr x21, =0x093f2b8dad5281e7
    and x21, x21, x19
    ldr x22, =0xcbab245135469bab
    and x22, x22, x19 
    adcs x11, x11, x20   
    adcs x12, x12, x21   
    adcs x13, x13, x22   
    ldr x20, =0x50cbaa75a2a1fa44
    and x20, x20, x19
    ldr x21, =0x10028248ad4fc4b1
    and x21, x21, x19
    adcs x14, x14, x20
    adcs x15, x15, x21
    ldr x22, =0x6b5bff7643c64f7a
    and x22, x22, x19
    ldr x20, =0x0000000000000008
    and x20, x20, x19
    adcs x16, x16, x22
    adcs x17, x17, x20
    stp x2, x3,   [x1,#0]
    stp x4, x5,   [x1,#16]
    stp x6, x7,   [x1,#32]
    stp x8, x9,   [x1,#48]
    stp x10, x11, [x1,#64]
    stp x12, x13, [x1,#80]
    stp x14, x15, [x1,#96]
    stp x16, x17, [x1,#112]
    ldp x15, x16, [sp]
    ldp x17, x18, [sp, #16]
    ldp x19, x20, [sp, #32]
    ldp x21, x22, [sp, #48]
    ldp x23, x24, [sp, #64]
    ldp x25, x26, [sp, #80]
    ldp x27, x28, [sp, #96]
    add sp, sp, #112
ret    

