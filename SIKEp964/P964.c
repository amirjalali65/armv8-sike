/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P964
*********************************************************************************************/

#include "P964_internal.h"

// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format).
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position.
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 964-bit field element is represented with Ceil(964 / 64) = 12 64-bit digits or Ceil(964 / 32) = 24 32-bit digits.

//
// Curve isogeny system "SIDHp964". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p964^2), where A=0, B=1, C=1 and p964 = 2^486*3^301-1
//

const uint64_t p964[NWORDS64_FIELD]   = {0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF,
                                         0xFFFFFFFFFFFFFFFF, 0x451CD4BFFFFFFFFF, 0xABB38EAB467ACDE5, 0xE56EF6AA57A94749, 0x093F2B8DAD5281E7, 0xCBAB245135469BAB,
                                         0x50CBAA75A2A1FA44, 0x10028248AD4FC4B1, 0x6B5BFF7643C64F7A, 0x0000000000000008};
const uint64_t p964p1[NWORDS64_FIELD] = {0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                         0x0000000000000000, 0x451CD4C000000000, 0xABB38EAB467ACDE5, 0xE56EF6AA57A94749, 0x093F2B8DAD5281E7, 0xCBAB245135469BAB,
                                         0x50CBAA75A2A1FA44, 0x10028248AD4FC4B1, 0x6B5BFF7643C64F7A, 0x0000000000000008};
const uint64_t p964x2[NWORDS64_FIELD] = {0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF,
                                         0xFFFFFFFFFFFFFFFF, 0x8A39A97FFFFFFFFF, 0x57671D568CF59BCA, 0xCADDED54AF528E93, 0x127E571B5AA503CF, 0x975648A26A8D3756,
                                         0xA19754EB4543F489, 0x200504915A9F8962, 0xD6B7FEEC878C9EF4, 0x0000000000000010};
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER] = {0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                              0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000004000000000};
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER] = {0xAD19EB3795147353, 0xA95EA51D26AECE3A, 0x36B54A079F95BBDA, 0x44D51A6EAC24FCAE, 
                                            0xD68A87E9132EAC91, 0x22B53F12C5432EA9, 0xD90F193DE8400A09, 0x0000000021AD6FFD};

// Alice's generator values {XPA0 + XPA1*i, XQA0, XRA0 + XRA1*i} in GF(p964^2), expressed in Montgomery representation
const uint64_t A_gen[5 * NWORDS64_FIELD] = {0x3205A18001E2A24B, 0x087DA50D0112D59D, 0x69DFC3FFB0C63715, 0xA00CAE7EF30E6DDC,
                                            0xEEDE43545DB0B42A, 0xC790BD6BF5619353, 0xA99D13238E101887, 0x960621FDA2636C00,
                                            0xDD218BC3785B5C5D, 0x97E517EFA176994D, 0xAEF106832644ACCF, 0x74FF91133E945E06,
                                            0xB5E6760EA6854AC0, 0x2D4044C2EFD7EA27, 0x21A6ECE47D629E3E, 0x0000000000000003, // XPA0
                                            0x95020EF3612F6103, 0x09D16814DC5A38FD, 0x69B5DCC5E35A438C, 0xF248F36527CDE8F2,
                                            0x82169CA2D8262FFD, 0x61F3DB2299A84BC8, 0xC9411CAB9B45A9D2, 0xF91224145D30FCBE,
                                            0xE0F91721D081755A, 0x844760C44893F4FD, 0xF1810657C2B7E041, 0x129EA85F8E342073,
                                            0xADDECBE88C1DF985, 0x1F12073934C874B3, 0xE4E9629F2FD1765F, 0x0000000000000006, // XPA1
                                            0x7F14FA2E03845DAE, 0x45EA605E1B86A894, 0x783698A71544E10A, 0x21ECA6EC4D3086A7,
                                            0x85EC8F42D4AD1505, 0x7DBCDC3DFDB60BA3, 0x5F2BDAF703E83BB9, 0x8615254DF178C1F7,
                                            0x7310B6D02E78DB88, 0xFE141B85D591243E, 0xF8F5F9B19C369F7A, 0xFDFF3A867BBB267B,
                                            0xEA94116F1E45C834, 0xB712A45EED7720D3, 0xCD5AEB4E02DBE2CE, 0x0000000000000001, // XQA0
                                            0xFBBC718642FB32A8, 0x79D7185E700AC201, 0x4F92A4C89A01CB59, 0x56AA1CCF17FBC00B,
                                            0x645B4995B9BFF792, 0x89636A809A2ADBF8, 0x196D235E0A3C2C47, 0x5D79C903F64060BC,
                                            0x9CC125E995B36C51, 0x45594E396A23678F, 0x7DF03D96C554F021, 0x6A19627476971A52,
                                            0x071F0B5C0ED97BD6, 0xE949811AB51AC64A, 0x320D0F8FD7DB6145, 0x0000000000000002, // XRA0
                                            0x7D035370A2C9FDA5, 0xB6AE3E236ABB05D6, 0xA6D2701A022FC5B0, 0x9D9D0A9386746015,
                                            0x6D69B5823F2021AE, 0xE4F3058E9F186B2F, 0x267ADC560581B29D, 0x50C2D1C409F0835C,
                                            0xACA19A3673421638, 0x07DE707A62554D4A, 0x71B87313ACF6AEF8, 0x0A20C18B53536924,
                                            0xD0798EA692B30404, 0x1F9234DB715779E5, 0xD15E38C73113FB67, 0x0000000000000004}; // XRA1
// Bob's generator values {XPB0 + XPB1*i, XQB0, XRB0 + XRB1*i} in GF(p964^2), expressed in Montgomery representation
const uint64_t B_gen[5 * NWORDS64_FIELD] = {0x2AB74EBDA9B62855, 0x9F86F8167A469032, 0xC795B5804514A045, 0x527C8F62AF7E0C25,
                                            0x0DB1E1420BCBC268, 0xFAA6B0E4AAEC2203, 0x3A670D66E561AA2E, 0xD1990F715F4AF064,
                                            0xC40C859D3A34EA79, 0xCD8726F217580B01, 0xE1859571885F2FEF, 0xE02B8A4175CEE50B,
                                            0x6ECBD05AF09A66A4, 0xBE2B6F3A9D8C8BF4, 0x138541CB804679E0, 0x0000000000000008, // XPB0
                                            0xA7E06E6859BB2967, 0xBFEDEA48A5861BCB, 0x7AEC252293B335C5, 0x6EA37528967BF5E0,
                                            0xFE02F925A4548B9D, 0x67118E065A992294, 0x76B7CD7B8C6E67B6, 0xFF630AF9BFC08A5B,
                                            0x695782EEA8A7AEC4, 0xDC4F028E37D61B50, 0x3F0A9046697B2368, 0x6003776BF0B1C9A8,
                                            0x946397CD3A3A3E58, 0x5ECD73D587810018, 0x82CD4A5231698510, 0x0000000000000005, // XPB1
                                            0xFDCD1FCF2FB36A8B, 0x8D6070AD207085D2, 0x5D7734640F8D64BD, 0xEA6DCC0D163A892F,
                                            0x9E48DAB491576BEA, 0xCF02834C9288ADFD, 0xFB3B031CBC2AA497, 0xF2AB6C7AE889C92C,
                                            0xC667688B2C89CE84, 0x10652320E1863721, 0x34AA3614751ED239, 0xF38D7983E47428E6,
                                            0x8E7011C52BE8F8D2, 0x5D95BD3DB486C845, 0x1C20C5875BEF325E, 0x0000000000000000, // XQB0
                                            0xAE75C752E142F6FC, 0x248D1A6B3AA61378, 0xC95584F9BD9533CC, 0xF601E155B389E3D7,
                                            0xF741E8FF22DD0EE4, 0x4A5754F586CBFB26, 0x6367E81451ECA66C, 0x69A4A30D1C7AD925,
                                            0x04897B7C624B3361, 0x28C17D284F033135, 0x00959A17D3BC0E49, 0xAEC6A480B1463225,
                                            0x003E25956B612ECB, 0xB2015FECCACC8711, 0x8B94229769EE788A, 0x0000000000000004, // XRB0
                                            0xAFB27F2AE1FAFD55, 0x77E129B8204ABE62, 0xCE7A0C030B309120, 0x6C3CEE3D639CF528,
                                            0x8A013B7829CD958F, 0xB935D14F3C27048B, 0xF3F606801FA76991, 0x026900C9B5F1A30C,
                                            0xC0B71E83A4518B71, 0xE98C946971C8B47B, 0x9DF5D432823DF88D, 0x69B107709CB6F44F,
                                            0x1946FFD85FD24AA3, 0xD036522C10D6EFE0, 0x0536217E4CA915B3, 0x0000000000000004}; // XRB1
// Montgomery constant Montgomery_R2 = (2^1024)^2 mod p964
const uint64_t Montgomery_R2[NWORDS64_FIELD] = {0x49BEAB58E287A9F1, 0x86A242D12EA5A11A, 0xB12512C92DD30800, 0x83D81C6BF9BC092E,
                                                0x5AFA951F0F780370, 0x36BAB97634D25944, 0x130638EF637CB27F, 0x0CDB7689C4452519,
                                                0xEF58B1B27E0542F8, 0x6109A63CB14FD223, 0x7E60914AB79A734B, 0x282640CE02415051,
                                                0x8E76DB8220153984, 0xDE7F01881434E82D, 0xB0AA8606EF7C00DE, 0x0000000000000006};
// Value one in Montgomery representation
const uint64_t Montgomery_one[NWORDS64_FIELD] = {0x1E67F3F78730D81B, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x3E738FC000000000,
                                                 0x0A21FE854AE027A8, 0x853FB7B621CC75C3, 0xDB8515C38E354F0C, 0x6232C1D569A850B5,
                                                 0x377C1534E2EF4915, 0x964EE9BA62138CAD, 0xE011B4007E0E88AA, 0x0000000000000004};

// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice] = {
0, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 5, 6, 7, 8, 8, 9, 9, 9, 9,
9, 9, 9, 12, 11, 12, 12, 13, 14, 15, 16, 16, 16, 16, 16, 16, 17, 17, 18, 18, 17,
21, 17, 18, 21, 20, 21, 21, 21, 21, 21, 22, 25, 25, 25, 26, 27, 28, 28, 29, 30,
31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 35, 36, 36, 33, 36, 35, 36, 36, 35,
36, 36, 37, 38, 38, 39, 40, 41, 42, 38, 39, 40, 41, 42, 40, 46, 42, 43, 46, 46,
46, 46, 48, 48, 48, 48, 49, 49, 48, 53, 54, 51, 52, 53, 54, 55, 56, 57, 58, 59,
59, 60, 62, 62, 63, 64, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 66, 67,
65, 66, 67, 66, 69, 70, 66, 67, 66, 69, 70, 69, 70, 70, 71, 72, 71, 72, 72, 74,
74, 75, 72, 72, 74, 74, 75, 72, 72, 74, 75, 75, 72, 72, 74, 75, 75, 77, 77, 79,
80, 80, 82, 82, 84, 83, 83, 84, 85, 86, 87, 86, 87, 90, 87, 86, 87, 87, 88, 86,
88, 88, 87, 88, 86, 88, 88, 89, 90, 91, 93, 98, 94, 103, 96, 98, 98, 99, 100,
101, 103, 103, 104, 105, 106, 106, 106, 109, 108, 106, 106, 109, 109, 113, 111,
113, 113, 113, 118, 115, 120};

const unsigned int strat_Bob[MAX_Bob] = {
0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 8, 8, 8, 9,
9, 9, 9, 9, 10, 12, 12, 12, 12, 12, 12, 13, 14, 14, 15, 16, 16, 16, 16, 16, 17,
16, 16, 17, 19, 19, 20, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 24, 24, 25,
27, 27, 28, 28, 29, 28, 29, 28, 28, 28, 30, 28, 28, 28, 29, 30, 33, 33, 33, 33,
34, 35, 37, 37, 37, 37, 38, 38, 37, 38, 38, 38, 38, 38, 39, 43, 38, 38, 38, 38,
43, 40, 41, 42, 43, 48, 45, 46, 47, 47, 48, 49, 49, 49, 50, 51, 50, 49, 49, 49,
49, 51, 49, 53, 50, 51, 50, 51, 51, 51, 52, 55, 55, 55, 56, 56, 56, 56, 56, 58,
58, 61, 61, 61, 63, 63, 63, 64, 65, 65, 65, 65, 66, 66, 65, 65, 66, 66, 66, 66,
66, 66, 66, 71, 66, 73, 66, 66, 71, 66, 73, 66, 66, 71, 66, 73, 68, 68, 71, 71,
73, 73, 73, 75, 75, 78, 78, 78, 80, 80, 80, 81, 81, 82, 83, 84, 85, 86, 86, 86,
86, 86, 87, 86, 88, 86, 86, 86, 86, 88, 86, 88, 86, 86, 86, 88, 88, 86, 86, 86,
93, 90, 90, 92, 92, 92, 93, 93, 93, 93, 93, 97, 97, 97, 97, 97, 97, 97, 97, 97,
97, 97, 97, 102, 102, 102, 102, 102, 107, 107, 107, 107, 107, 112, 112, 112,
112, 112, 110, 111, 112, 113, 114, 115, 114, 114, 114, 115, 114, 115, 114, 115,
115, 116, 115, 115, 117, 116, 115, 122, 117, 116, 116, 115, 122, 117, 116, 122,
115, 122, 117, 116, 122, 115, 122, 117, 116, 122, 115};

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions

#define fpcopy fpcopy964
#define fpzero fpzero964
#define fpadd fpadd964
#define fpsub fpsub964
#define fpneg fpneg964
#define fpdiv2 fpdiv2_964
#define fpcorrection fpcorrection964
#define fpmul_mont fpmul964_mont
#define fpsqr_mont fpsqr964_mont
#define fpinv_mont fpinv964_mont
#define fpinv_chain_mont fpinv964_chain_mont
#define fpinv_mont_bingcd fpinv964_mont_bingcd
#define fp2copy fp2copy964
#define fp2zero fp2zero964
#define fp2add fp2add964
#define fp2sub fp2sub964
#define fp2neg fp2neg964
#define fp2div2 fp2div2_964
#define fp2correction fp2correction964
#define fp2mul_mont fp2mul964_mont
#define fp2sqr_mont fp2sqr964_mont
#define fp2inv_mont fp2inv964_mont
#define fp2inv_mont_bingcd fp2inv964_mont_bingcd
#define fpequal_non_constant_time fpequal964_non_constant_time
#define mp_add_asm mp_add964_asm
#define mp_addx2_asm mp_add964x2_asm
#define mp_subx2_asm mp_sub964x2_asm

#include "fpx.c"
#include "ec_isogeny.c"
#include "sidh.c"
#include "sike.c"
