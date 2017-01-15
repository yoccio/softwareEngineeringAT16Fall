//
//  main.cpp
//  test
//
//  Created by ZhaoYoccio on 16/9/28.
//  Copyright © 2016年 ZhaoYoccio. All rights reserved.
//

#include "mini-gmp.h"

#include <stdio.h>
#include <stdlib.h>

#define OP_ADD 1
#define OP_SUB 2
#define OP_MUL 4
#define OP_DIV 8
FILE *fp;
FILE *fpo;
char temp1[1000000] = {0};
char temp2[1000000] = {0};
char temp3[1000000] = {0};
char lhs_str[1000000] = {0};
char rhs_str[1000000] = {0};
char ans_str[1000000] = {0};
char ans1_str[1000000] = {0};
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
        retL = sscanf(temp1,"%s%c",lhs_str,&has_space_worng_L);
        retR = sscanf(temp2,"%s%c",rhs_str,&has_space_worng_R);
        
        if(retL != 1 || retR != 1){
            error();
            continue;
        }
        retL = mpz_init_set_str(lhs, lhs_str, 10);
        retR = mpz_init_set_str(rhs, rhs_str, 10);
        if (retL == -1 || retR == -1) {
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
            fprintf(fpo,"%s\n",mpz_get_str(ans_str,10, ans));
        }
        else if (retval == OP_ADD) {
            mpz_add(ans,lhs,rhs);
            fprintf(fpo,"%s\n",mpz_get_str(ans_str,10, ans));
        }
        else if (retval == OP_MUL) {
            mpz_mul(ans,lhs,rhs);
            fprintf(fpo,"%s\n",mpz_get_str(ans_str,10, ans));
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
            fprintf(fpo,"%s %s\n",mpz_get_str(ans_str,10, ans),mpz_get_str(ans1_str, 10,reminder));
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
