# armv8-sike
Highly-optimized ARMv8 implementation of Supersingular Isogeny Key Encapsulation (SIKE).

The efficient implementation of SIKE protocol on ARMv8 high-performance processors. 
The finite field arithmetic implementation is developed by the state-of-the-art implementation techniques, taking advantage of ARMv8 64-bit general purpose registers combined with ASIMD vectorization. The field multiplication is designed and developed using one- and two-level additive Karatsuba method. The independent multiplications are implemented using both AArch64 and ASIMD hand-crafted assembly using an interleaved technique  to maximize the pipeline throuhput and efficiency of the library. 

The submitted SIKE proposal contains the optimized implementation of SIKEp503 and SIKEp751 on different platforms. This repositoy contains the highly-optimized implementation of SIKEp503, SIKEp751, and SIKEp964 on ARMv8 platforms. 

## Content
[SIKEp503](https://github.com/amirjalali65/armv8-sike/tree/master/SIKEp503): Optimized implementaion of SIKEp503 using only 64-bit general registers (Designed and developed by Matthew Campagna)

[SIKEp503_mixed](https://github.com/amirjalali65/armv8-sike/tree/master/SIKEp503_mixed): Optimized implementation of SIKEp503 using the mixture of general registers and ASIMD vectorization hand-written assembly.

[SIKEp751](https://github.com/amirjalali65/armv8-sike/tree/master/SIKEp751): Optimized implementaion of SIKEp751 using only 64-bit general registers (Designed and developed by Matthew Campagna)

[SIKEp751_mixed](https://github.com/amirjalali65/armv8-sike/tree/master/SIKEp751_mixed): Optimized implementation of SIKEp751 using the mixture of general registers and ASIMD vectorization hand-written assembly.

[SIKEp964_mixed](https://github.com/amirjalali65/armv8-sike/tree/master/SIKEp964_mixed): Optimized implementation of SIKEp964 using the mixture of general registers and ASIMD vectorization hand-written assembly.

## Builing Binaries

### Cross Compilation for ARMv8 on Linux
ARMv8 executables can be generated using cross-compilation on Linux. There are different methods for cross-compilation. An easy approach is to install `gcc-aarch64-linux-gnu` package by executing:
```sh
$  sudo apt-get install gcc-aarch64-linux-gnu
```
After installation, simply use the following command to generate the ARMv8 executables:
```sh
$ make CC=aarch64-linux-gnu-gcc ARCH=ARM64
```
Now, the generated binaries can be run on ARMv8-A cores. 

## Contributors
* All the mixed arithmetic libraries are designed and developed by [Amir Jalali]()
* The SIKE optimized implementation (portable) is designed and developed by Microsoft Research and the SIKE team.
* The arithmetic libraries using general registers are designed and developed by Matthew Campagna

