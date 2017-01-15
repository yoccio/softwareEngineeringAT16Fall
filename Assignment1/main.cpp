//
//  main.cpp
//  largeScaleNumber
//
//  Created by ZhaoYoccio on 16/9/23.
//  Copyright © 2016年 ZhaoYoccio. All rights reserved.
//
//
//
//  note: use GNU GMP LIB (V6.1.1). you can get it from https://gmplib.org/
//  remcommend to use default configuation to get a good enough efficiency if you are not familiar with it.
//  If you can't successfully install gmp lib, an alternative lib is boost/multiprecision.You can add '#define _USING_BOOST' and delate '#define _USING_GMP' to change lib .In the wikipedia,I found BOOST/MUTIPRECISION using gmp as it's backends.But maybe is less efficiency then gmp lib for some language properties of C++.
//  usage: largeScaleNumber [input file name] [output file name]
//
//#define _USING_BOOST
#define _USING_BOOST
#include <stdio.h>
//#ifdef _USING_GMP
//#include <gmp.h>
//#endif
#include <stdlib.h>
#ifdef _USING_BOOST
#include <boost/multiprecision/gmp.hpp>
using namespace boost::multiprecision;
#endif
#define OP_ADD 1
#define OP_SUB 2
#define OP_MUL 4
#define OP_DIV 8
FILE *fp;
FILE *fpo;
char temp1[100000000] = {0};
char temp2[100000000] = {0};
char temp3[100000000] = {0};
void error(void){
    fprintf(fpo,"Error input\n");
}
int parseOp(void){
    int retval;
    char op;
    sscanf(temp3,"%c",&op);
    switch (op) {
        case '+':
            retval = OP_ADD;
            break;
        case '-':
            retval = OP_SUB;
            break;
        case '*':
            retval = OP_MUL;
            break;
        case '/':
            retval = OP_DIV;
            break;
        default:
            retval = -1;
            break;
    }
    return retval;
}
int parseCmdArg(int argc,const char * argv[]){
    if (argc != 3) {
        printf("please input a file name as the input file such as \"input.txt\"\n");
        return -1;
    }
    return 0;
}
int u_fscanfs(FILE *,char temp[]){
    long long i = 0;
    char u_temp;
    //step1 get first character
    u_temp = fgetc(fp);
    if (u_temp == EOF) {
        return EOF;
    }
    if (u_temp != '\n') {
        temp[i] = u_temp;
        i++;
    }else{
        temp[i] = '\0';
        return 0;
    }
    //step2 get others
    while ((u_temp = fgetc(fp)) != '\n' && u_temp != EOF) {
        temp[i] = u_temp;
        i++;
    }
    temp[i] = '\0';
    return 1;
}
int main(int argc, const char * argv[]) {
    mpz_t lhs,rhs,ans,reminder;
    int retval = 1;
    int retL = 1;
    int retR = 1;
    mpz_init(lhs);
    mpz_init(rhs);
    mpz_init(ans);
    mpz_init(reminder);
    retval = parseCmdArg(argc, argv);
    if (retval != 0) {
        return -1;
    }
    fp = fopen(argv[1],"r");
    fpo = fopen(argv[2], "w");
    if (!fp) {
        printf("con't open the input file\n");
        return -1;
    }
    while (1) {
        retval = u_fscanfs(fp,temp1);
        if (retval == EOF) {
            break;
        }
        retL = u_fscanfs(fp,temp2);
        retR = u_fscanfs(fp,temp3);
        if(retL != 1 || retR != 1){
            //just to handle empty line
            error();
            continue;
        }
        // to handel space in the middle of line
        char has_space_worng_L = 0;
        char has_space_worng_R = 0;
        retL = gmp_sscanf(temp1,"%Zd%c",lhs,&has_space_worng_L);
        retR = gmp_sscanf(temp2,"%Zd%c",rhs,&has_space_worng_R);

        if(retL != 1 || retR != 1){
            error();
            continue;
        }

        retval = parseOp();
        if (retval == -1) {
            error();
            continue;
        }

        if (retval == OP_SUB) {
            mpz_sub(ans,lhs,rhs);
            gmp_fprintf(fpo,"%Zd\n",ans);
        }
        else if (retval == OP_ADD) {
            mpz_add(ans,lhs,rhs);
            gmp_fprintf(fpo,"%Zd\n",ans);
        }
        else if (retval == OP_MUL) {
            mpz_mul(ans,lhs,rhs);
            gmp_fprintf(fpo,"%Zd\n",ans);
        }
        else if (retval == OP_DIV) {
            if (!mpz_sgn(rhs)) {
                //if divisor is zero
                error();
                continue;
            }
            if(mpz_sgn(rhs) == -1 || mpz_sgn(lhs) == -1){
                error();
                continue;
            }
            mpz_tdiv_qr(ans,reminder,lhs,rhs);
            gmp_fprintf(fpo,"%Zd %Zd\n",ans,reminder);
        }
    }
    mpz_clear(ans);
    mpz_clear(lhs);
    mpz_clear(rhs);
    mpz_clear(reminder);
    fclose(fp);
    fclose(fpo);
    return  0;
}
