# nthPrime_P2
----
## 现代软件工程－中科－16秋第二次作业

###问题：
>For the given number N, you should give the Nth prime.  
>####input:  
One number N in decimal, 1 <= N <= 100,000,000  
>####output:  
The Nth prime in decimal (<=2,038,074,743, so 32 bit integer is enough, do not need large number arithmetic)  

###算法描述：  
主要用了[Sieve of Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes#Euler.27s_Sieve)
>Input: an integer n > 1
> 
>Let A be an array of Boolean values, indexed by integers 2 to n,
>initially all set to true.
> 
>for i = 2, 3, 4, ..., not exceeding √n:
>  if A[i] is true:
>    for j = i2, i2+i, i2+2i, i2+3i, ..., not exceeding n :
>      A[j] := false
> 
>Output: all i such that A[i] is true.  

改进了一部分：  
- 去除偶数
- 使用1位代表一个数降低了内存需求
- 增大了内存循环步长  

__时间复杂度 O(nlog(log(n)))__  
__空间复杂度 O(nlog(n))__  

###测试环境：  
开发环境 xcode7.0 编译器 APPLE LLVM7.1  
OSX10.11.6 g++ 测试通过  
无VS测试  

----
###update1:增加一份多线程实现，测试效率约为单线程2倍  

文件放在了 /mth 下，仅对4线程做了优化 1亿时耗时约7000ms
