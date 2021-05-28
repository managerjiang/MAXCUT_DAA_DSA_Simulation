#include <stdio.h>
#include "time.h"
#include <math.h>
#define DIM 18
double norm1_f(double* set, int num);
void main(int argc, char* argv[]) {

    int i, j;
    clock_t start, finish;     //定义第一次调用CPU时钟单位的实际，可以理解为定义一个计数器
    double Total_time;
    double C1[DIM][DIM] = { 0.0 };
    C1[0][2] = 1.0 / 3;
    C1[2][0] = 1.0 / 3;
    C1[0][12] = 0.2;
    C1[2][12] = 1.2;
    C1[12][0] = 0.2;
    C1[12][2] = 1.2;

    double C2[DIM][DIM] = { 0.0};
    C2[0][1] = 0.5;
    C2[0][3] = 0.25;
    C2[1][3] = 1.0 / 3;
    C2[1][0] = 0.5;
    C2[3][0] = 0.25;
    C2[3][1] = 1.0 / 3;
    C2[9][13] = 0.5;
    C2[13][9] = 0.5;
    C2[1][9] = 1.5;
    C2[9][1] = 1.5;
    double C3[DIM][DIM] = { 0.0 };
    C3[3][4] = 1.0 / 2;
    C3[4][3] = 1.0 / 2;
    C3[3][13] = 0.8;
    C3[13][3] = 0.8;
    C3[4][17] = 1.1;
    C3[17][4] = 1.1;
    double C4[DIM][DIM] = { 0.0 };
    C4[2][3] = 1.0 / 4;
    C4[3][2] = 1.0 / 4;
    C4[2][15] = 1.5;
    C4[15][2] = 1.5;
    C4[3][15] = 3.0;
    C4[15][3] = 3.0;
    double C5[DIM][DIM] = { 0.0 };
    C5[2][5] = 1.0 / 7;
    C5[2][6] = 1.0 / 2;
    C5[5][6] = 1.0 / 4;
    C5[5][2] = 1.0 / 7;
    C5[6][2] = 1.0 / 2;
    C5[6][5] = 1.0 / 4;
    C5[10][12] = 1.0 / 3;
    C5[12][10] = 1.0 / 3;
    C5[6][14] = 1.2;
    C5[14][6] = 1.2;
    double C6[DIM][DIM] = { 0.0 };
    C6[2][7] = 1.0 / 6;
    C6[7][2] = 1.0 / 6;
    C6[11][2] = 0.1;
    C6[2][11] = 0.1;
    C6[11][16] = 2.5;
    C6[16][11] = 2.5;
    C6[7][12] = 2.32;
    C6[12][7] = 2.32;
    C6[2][8] = 1.2;
    C6[8][2] = 1.2;
    double C[DIM][DIM] = { 0.0};
    for (i = 0;i < DIM;i++)
    {
        for (j = 0;j < DIM;j++)
        {
            C[i][j] = C1[i][j] + C2[i][j] + C3[i][j] + C4[i][j] + C5[i][j] + C6[i][j];
        }
    }
    /*集中式更新*/
    double Vv_var[DIM][DIM];
    double Vv_var_next[DIM][DIM];
    for (i = 0;i < DIM;i++)
    {
        for (j = 0;j < DIM;j++)
        {
            if (i == j)
                Vv_var[i][j] = 1.0;
            else
                Vv_var[i][j] = 0.0;
        }
    }
    for (i = 0;i < DIM;i++)
    {
        for (j = 0;j < DIM;j++)
        {
            Vv_var_next[i][j] = Vv_var[i][j];
        }
    }
    start = clock();
    int iteration_c = 600;
    double theta_c =  0.1854; 
    /*迭代开始*/
    for (int ite = 0;ite < iteration_c;ite++)
    {
        for (i = 0;i < DIM;i++)
        {
            double sumj[DIM] = { 0 };
            for (j = 0;j < DIM;j++)
            {
                if (j < i)
                {
                    for (int jj = 0;jj < DIM;jj++)
                        sumj[jj] = sumj[jj] + C[i][j] * Vv_var_next[jj][j];
                }
                else
                {
                    for (int jj = 0;jj < DIM;jj++)
                        sumj[jj] = sumj[jj] + C[i][j] * Vv_var[jj][j];
                }
            }
            double gv_c[DIM] = { 0.0 };
            for (int jj = 0;jj < DIM;jj++)
            {
                gv_c[jj] = Vv_var[jj][i] - theta_c * sumj[jj];
            }
            for (int jj = 0;jj < DIM;jj++)
            {
                Vv_var[jj][i] = Vv_var_next[jj][i];
            }
            double* gvt = gv_c;
            for (int jj = 0;jj < DIM;jj++)
            {
                Vv_var_next[jj][i] = gv_c[jj] * 1.0 / norm1_f(gvt, DIM);
            }
        }

    }
    double Xm[DIM][DIM];
    for (i = 0;i < DIM;i++) {
        for (j = 0;j < DIM;j++) {
            double sumn = 0.0;
            for (int k = 0;k < DIM;k++) {

                sumn = sumn + Vv_var_next[k][i] * Vv_var_next[k][j];////* (*(result_matrix + i) + j)
            }
            Xm[i][j] = sumn;
        }
    }
    printf("X matrix is################\n");
    for (i = 0;i < DIM;i++)
    {
        for (j = 0;j < DIM;j++)
        {
            printf("%f, ", Xm[i][j]);
        }
        printf("\n");
    }
    double Rs[DIM][DIM];
    for (i = 0;i < DIM;i++) {
        for (j = 0;j < DIM;j++) {
            double sumn = 0.0;
            for (int k = 0;k < DIM;k++) {

                sumn = sumn + C[k][i] * Xm[k][j];////* (*(result_matrix + i) + j)
            }
            Rs[i][j] = sumn;
        }
    }
    double fvalue = 0.0;
    for (i = 0;i < DIM;i++)
    {
        for (j = 0;j < DIM;j++)
        {
            if (i == j)
                fvalue = fvalue + Rs[i][j];
        }
    }
    double flim = -36.98024;//-4.880952;
    printf("center f function value is %f!!!!!!!!!!!\n", fvalue);
    printf("center f-f* function value is %f!!!!!!!!!!!\n", fvalue - flim);
    finish = clock();
    Total_time = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\n集中式算法函数运行时间：%0.3f毫秒 \n", Total_time);
}
double norm1_f(double* set, int num)
{
    double norm_value = 0.0;
    for (int i = 0;i < num;i++)
    {
        double a = *(set + i);
        norm_value = norm_value + a * a;
    }
    norm_value = sqrt(norm_value);
    return norm_value;
}
