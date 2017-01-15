//
//  main.cpp
//  nthPrime_P2
//
//  Created by ZhaoYoccio on 16/10/12.
//  Copyright © 2016年 ZhaoYoccio. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdint>

#define INT_WORD_LENGTH uint64_t //64位字长类型定义
#define MACHINE_LENGTH 64 //假定使用64位字长
#define DOUBLE_MACHINE_LENGTH 128 //两倍字长，因很多地方需要去除偶数处理后字长
#define MASK_TRUE 0x1LL //64位字长素数掩码
#define MASK_MACHINE_LENGTH 0xFFFFLL

#define INT_TO_BITS(indexOfInt) (indexOfInt / (DOUBLE_MACHINE_LENGTH))
#define MASK_INT_TO_BITS(indexOfInt) (MASK_TRUE << ( (indexOfInt >> 1) & MASK_MACHINE_LENGTH))


uint64_t nthPrime(uint32_t nth){
    //处理输入1
    if (nth == 1) {
        return 2;
    }
    
    //初始化工作
    //将数组大小换算为64位按bit存储的大小，
    //note:仅需保存奇数
    
    //pn<n*ln(n)+n*ln(ln(n)) (n≥6)
    //https://en.wikipedia.org/wiki/Prime_number_theorem
    uint64_t maxNum = nth < 6 ? 25 : nth * log(nth) + nth * log(log(nth));
    
    const uint64_t arraySize = maxNum / (MACHINE_LENGTH << 1) + 1;
    INT_WORD_LENGTH * isPrime = new INT_WORD_LENGTH[arraySize];
    memset(isPrime, 0xFF, arraySize * sizeof(INT_WORD_LENGTH));
    
    uint32_t count = 1;
    uint64_t n = 0;
    
    //Eratosthenes 筛法求素数
    for(INT_WORD_LENGTH num = 3;
        num < (INT_WORD_LENGTH)sqrt(maxNum) + 1;
        num += 2){
        if (isPrime[INT_TO_BITS(num)] & MASK_INT_TO_BITS(num)) {
            //如果 num 是素数，做筛选
            //begin at num^2 ,not exceeding maxNum
            INT_WORD_LENGTH t_max_num = num * num;
            INT_WORD_LENGTH t_double_num = 2 * num;
            for (INT_WORD_LENGTH isNOT_Prime = t_max_num;
                 isNOT_Prime < maxNum;
                 isNOT_Prime += t_double_num) {
                // 因已去除偶数，故只对合数为奇数时做处理
                isPrime[INT_TO_BITS(isNOT_Prime)] &= ~ MASK_INT_TO_BITS(isNOT_Prime);
            }
        }
    }
    for (INT_WORD_LENGTH num = 3;
         num < maxNum;
         num += 2) {
        if (isPrime[INT_TO_BITS(num)] & MASK_INT_TO_BITS(num)) {
            count++;
            if (count == nth) {
                n = num;
                break;
            }
        }
    }
    delete [] isPrime;
    return n;
}
int main(int argc, const char * argv[]) {
    uint32_t nth;
    if (std::cin >> nth && nth > 0 && nth <= 100000000) {
        std::cout << nthPrime(nth) <<std::endl;
    }
    else{
        std::cout << "error input" << std::endl;
    }
    return 0;
}
