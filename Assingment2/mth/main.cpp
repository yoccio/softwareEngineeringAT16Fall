//
//  main.cpp
//  test2
//
//  Created by ZhaoYoccio on 16/10/1.
//  Copyright © 2016年 ZhaoYoccio. All rights reserved.
//

#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <future>
#include <thread>
#include <fstream>
#include <chrono>

using namespace std;

#define MACHINE_LENGTH 64 //假定使用64位字长
#define DOUBLE_MACHINE_LENGTH 128 //两倍字长，因很多地方需要去除偶数处理后字长
#define MASK_TRUE 0x1LL //64位字长素数掩码
#define MASK_MACHINE_LENGTH 0xFFFFLL
#define INT_WORD_LENGTH uint64_t //64位字长类型定义
#define INT_TO_BITS(indexOfInt) ((indexOfInt) / (DOUBLE_MACHINE_LENGTH))
#define MASK_INT_TO_BITS(indexOfInt) (MASK_TRUE << ( ((indexOfInt) >> 1) & MASK_MACHINE_LENGTH))

int g_threadLock[4] = {1,1,1,1};

inline uint64_t find_start_NOT_prime(uint64_t num, uint64_t bgn){
    //使用num做筛子筛选bgn开始数组时，起始筛选数
    int base = (int)(bgn / num);
    int reminder = (int)(bgn % num);
    if (base < num) {
        return  num * num;
    }
    else{
        if (reminder == 0) {
            if (base % 2 == 1) {
                return bgn;
            }
            else{
                return bgn + num;
            }
        }
        else{
            if (base % 2 == 1) {
                return num * (base + 2);
            }
            else{
                return num * (base + 1);
            }
        }
    }
}
void do_sieve(uint64_t bgn, uint64_t end, uint64_t * isPrime,
                     uint64_t end_sm, uint64_t * prime){
    //用prime (0,end_sm)做筛子，筛除isPrime 中 [bgn,end）段
    for (INT_WORD_LENGTH num = 3;
         num < end_sm;
         num += 2) {
        if(prime[INT_TO_BITS(num)] & MASK_INT_TO_BITS(num)){
            INT_WORD_LENGTH t_double_num = 2 * num;
            
            for (INT_WORD_LENGTH isNOT_Prime = find_start_NOT_prime(num, bgn);
                 isNOT_Prime < end;
                 isNOT_Prime += t_double_num) {
                assert(isNOT_Prime - bgn >= 0);
                isPrime[INT_TO_BITS(isNOT_Prime - bgn)] &= ~ MASK_INT_TO_BITS(isNOT_Prime - bgn);
            }
        }
    }
}
uint32_t parallel_sieve(int nthread,uint64_t bgn,uint64_t end, uint64_t * isPrime, uint64_t end_sm,uint64_t *prime){
    do_sieve(bgn, end, isPrime, end_sm, prime);
    uint32_t count = 0;
    if (nthread != 4) {
        for (INT_WORD_LENGTH num = ((bgn % 2 == 1) ? bgn : bgn + 1);
             num < end;
             num += 2) {
            if (isPrime[INT_TO_BITS(num - bgn)] & MASK_INT_TO_BITS(num - bgn)) {
                (count)++;
            }
        }
    }
    g_threadLock[nthread - 1] = 0;
    return count;
}
uint64_t thread_sieve(INT_WORD_LENGTH nth){
    //处理小输入
    uint64_t smalln[10] = {2,3,5,7,11,13,17,19,23,29};
    if (nth < 10){
        return smalln[nth - 1];
    }
    
    //初始化工作
    //将数组大小换算为64位按bit存储的大小，
    //note:仅需保存奇数
    uint32_t maxNum = nth < 6 ? 25 : nth * log(nth) + nth * log(log(nth));
    
    //做筛子的小数组，内嵌在大数组内
    const uint32_t small_sieve_size = sqrt(maxNum) + 1;
    
    //首先计算各线程要算的数组范围
    uint64_t length_each_thread = (maxNum - small_sieve_size) / 4 + 1;
    uint64_t thread1_bgn = small_sieve_size;
    uint64_t thread1_end = thread1_bgn + length_each_thread;
    uint64_t thread2_bgn = thread1_end;
    uint64_t thread2_end = thread2_bgn + length_each_thread;
    uint64_t thread3_bgn = thread2_end;
    uint64_t thread3_end = thread3_bgn + length_each_thread;
    uint64_t thread4_bgn = thread3_end;
    uint64_t thread4_end = maxNum;
    
    const uint64_t arraySize = length_each_thread / (MACHINE_LENGTH << 1) + 1;
    INT_WORD_LENGTH * isPrime1;
    INT_WORD_LENGTH * isPrime2;
    INT_WORD_LENGTH * isPrime3;
    INT_WORD_LENGTH * isPrime4;
    if (arraySize == 1) {
        isPrime1 = new INT_WORD_LENGTH;
        isPrime2 = new INT_WORD_LENGTH;
        isPrime3 = new INT_WORD_LENGTH;
        isPrime4 = new INT_WORD_LENGTH;
    }
    else{
        isPrime1 = new INT_WORD_LENGTH[arraySize];
        isPrime2 = new INT_WORD_LENGTH[arraySize];
        isPrime3 = new INT_WORD_LENGTH[arraySize];
        isPrime4 = new INT_WORD_LENGTH[arraySize];
    }
    INT_WORD_LENGTH * small_prime1 = new INT_WORD_LENGTH[small_sieve_size];
    INT_WORD_LENGTH * small_prime2 = new INT_WORD_LENGTH[small_sieve_size];
    INT_WORD_LENGTH * small_prime3 = new INT_WORD_LENGTH[small_sieve_size];
    INT_WORD_LENGTH * small_prime4 = new INT_WORD_LENGTH[small_sieve_size];

    memset(isPrime1, 0xFF, arraySize * sizeof(INT_WORD_LENGTH));
    memset(isPrime2, 0xFF, arraySize * sizeof(INT_WORD_LENGTH));
    memset(isPrime3, 0xFF, arraySize * sizeof(INT_WORD_LENGTH));
    memset(isPrime4, 0xFF, arraySize * sizeof(INT_WORD_LENGTH));
    memset(small_prime1, 0xFF, small_sieve_size * sizeof(INT_WORD_LENGTH));

    uint32_t count = 1;
    
    uint64_t n = 0;
    
    //Eratosthenes 筛法求素数,先算出一个小筛子
    for(INT_WORD_LENGTH num = 3;
        num < (INT_WORD_LENGTH)sqrt(small_sieve_size) + 1;
        num += 2){
        if (small_prime1[INT_TO_BITS(num)] & MASK_INT_TO_BITS(num)) {
            //如果 num 是素数，做筛选
            //begin at num^2 ,not exceeding maxNum
            
            INT_WORD_LENGTH t_max_num = num * num;
            INT_WORD_LENGTH t_double_num = 2 * num;
            
            for (INT_WORD_LENGTH isNOT_Prime = t_max_num;
                 isNOT_Prime < small_sieve_size;
                 isNOT_Prime += t_double_num) {
                // 因已去除偶数，故只对合数为奇数时做处理
                small_prime1[INT_TO_BITS(isNOT_Prime)] &= ~ MASK_INT_TO_BITS(isNOT_Prime);
            }
        }
    }
    //统计小筛子内的素数个数
    for (INT_WORD_LENGTH num = 3;
         num < small_sieve_size;
         num += 2) {
        if (small_prime1[INT_TO_BITS(num)] & MASK_INT_TO_BITS(num)) {
            count++;
        }
    }
    memcpy(small_prime2, small_prime1, small_sieve_size * sizeof(INT_WORD_LENGTH));
    memcpy(small_prime3, small_prime1, small_sieve_size * sizeof(INT_WORD_LENGTH));
    memcpy(small_prime4, small_prime1, small_sieve_size * sizeof(INT_WORD_LENGTH));
    //四线程做筛法
    
    auto retVal1(std::async(parallel_sieve, 1, thread1_bgn, thread1_end, isPrime1, small_sieve_size, small_prime1));
    auto retVal2(std::async(parallel_sieve, 2, thread2_bgn, thread2_end, isPrime2, small_sieve_size, small_prime2));
    auto retVal3(std::async(parallel_sieve, 3, thread3_bgn, thread3_end, isPrime3, small_sieve_size, small_prime3));
    auto retVal4(std::async(parallel_sieve, 4, thread4_bgn, thread4_end, isPrime4, small_sieve_size, small_prime4));
    this_thread::yield();
    retVal4.wait();
    count += (retVal1.get() + retVal2.get() + retVal3.get());
    
    for (INT_WORD_LENGTH num = ((thread4_bgn % 2 == 1) ? thread4_bgn : thread4_bgn + 1);
         num < maxNum;
         num += 2) {
        if (isPrime4[INT_TO_BITS(num - thread4_bgn)] & MASK_INT_TO_BITS(num - thread4_bgn)) {
            count++;
            if (count == nth) {
                n = num;
                break;
            }
        }
    }
    
    delete [] isPrime1;
    delete [] isPrime2;
    delete [] isPrime3;
    delete [] isPrime4;
    delete [] small_prime1;
    delete [] small_prime2;
    delete [] small_prime3;
    delete [] small_prime4;
    
    return n;
}

int main(int argc, const char * argv[]) {
    uint64_t nth;
    if (std::cin >> nth && nth > 0 && nth <= 100000000) {
        std::cout << thread_sieve(nth) <<std::endl;
    }
    else{
        std::cout << "error input" << std::endl;
    }
    return 0;
}
