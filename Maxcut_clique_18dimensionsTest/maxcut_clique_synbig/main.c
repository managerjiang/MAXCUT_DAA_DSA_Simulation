#include "mpi.h"
#include <stdio.h>
#include "time.h"
#include <math.h>
#define AGENT 6   ////���и���������-n 6��ʵ��agent 5��0�����ڷַ����ݣ��ռ�����
#define DIM 18
double norm1_f(double* set, int num);
void main(int argc, char* argv[]) {
    double C1[DIM][DIM] = { 0.0 };
    int i, j;

    C1[0][2] = 1.0 / 3;
    C1[2][0] = 1.0 / 3;
    C1[0][12] = 0.2;
    C1[2][12] = 1.2;
    C1[12][0] = 0.2;
    C1[12][2] = 1.2;

    double C2[DIM][DIM] = { 0.0 };
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
    double C[DIM][DIM] = { 0.0 };
    for (i = 0;i < DIM;i++)
    {
        for (j = 0;j < DIM;j++)
        {
            C[i][j] = C1[i][j] + C2[i][j] + C3[i][j] + C4[i][j] + C5[i][j] + C6[i][j];
        }
    }
    int m4[6] = { 2, 5, 6,10,12,14 };
    int m5[6] = { 2, 7,8,11,12,16 };
    int m2[5] = { 0, 2, 3,12,15 };
    int m2f1[3] = { 0, 2 ,12};
    int m2f2[3] = { 2, 3,15 };
    int m3[4] = { 3, 4,13,17 };
    int m1[5] = { 0, 1, 3,9,13 };
    clock_t start, finish;     //�����һ�ε���CPUʱ�ӵ�λ��ʵ�ʣ���������Ϊ����һ��������
    int myrank, numprocs;
    double Total_time;        //����һ��double���͵ı��������ڴ洢ʱ�䵥λ

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        MPI_Ssend(C1, DIM * DIM, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);
        MPI_Ssend(C2, DIM * DIM, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
        MPI_Ssend(C3, DIM * DIM, MPI_DOUBLE, 3, 3, MPI_COMM_WORLD);
        MPI_Ssend(C4, DIM * DIM, MPI_DOUBLE, 2, 4, MPI_COMM_WORLD);
        MPI_Ssend(C5, DIM * DIM, MPI_DOUBLE, 4, 5, MPI_COMM_WORLD);
        MPI_Ssend(C6, DIM * DIM, MPI_DOUBLE, 5, 6, MPI_COMM_WORLD);
        MPI_Ssend(m1, 5, MPI_INT, 1, 11, MPI_COMM_WORLD);
        MPI_Ssend(m2f1, 3, MPI_INT, 2, 12, MPI_COMM_WORLD);
        MPI_Ssend(m2f2, 3, MPI_INT, 2, 13, MPI_COMM_WORLD);
        MPI_Ssend(m3, 4, MPI_INT, 3, 14, MPI_COMM_WORLD);
        MPI_Ssend(m4, 6, MPI_INT, 4, 15, MPI_COMM_WORLD);
        MPI_Ssend(m5, 6, MPI_INT, 5, 16, MPI_COMM_WORLD);
    }
    double C_local[2][DIM][DIM];
    int m_local[2][6];///��ʼ��
    for (i = 0;i < 2;i++)
    {
        for (j = 0;j < 6;j++)
        {
            m_local[i][j] = -1;
        }
    }
    MPI_Status status_c[2], status_m[2];
    if (myrank == 1)
    {
        MPI_Recv(C_local[0], DIM * DIM, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status_c[0]);
        MPI_Recv(m_local[0], 5, MPI_INT, 0, 11, MPI_COMM_WORLD, &status_m[0]);
    }
    if (myrank == 2)
    {
        MPI_Recv(C_local[0], DIM * DIM, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status_c[0]);
        MPI_Recv(C_local[1], DIM * DIM, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status_c[1]);
        MPI_Recv(m_local[0], 3, MPI_INT, 0, 12, MPI_COMM_WORLD, &status_m[0]);
        MPI_Recv(m_local[1], 3, MPI_INT, 0, 13, MPI_COMM_WORLD, &status_m[1]);
    }
    if (myrank == 3)
    {
        MPI_Recv(C_local[0], DIM * DIM, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status_c[0]);
        MPI_Recv(m_local[0], 4, MPI_INT, 0, 14, MPI_COMM_WORLD, &status_m[0]);
    }
    if (myrank == 4)
    {
        MPI_Recv(C_local[0], DIM * DIM, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status_c[0]);
        MPI_Recv(m_local[0], 6, MPI_INT, 0, 15, MPI_COMM_WORLD, &status_m[0]);
    }
    if (myrank == 5)
    {
        MPI_Recv(C_local[0], DIM * DIM, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &status_c[0]);
        MPI_Recv(m_local[0], 6, MPI_INT, 0, 16, MPI_COMM_WORLD, &status_m[0]);
    }
   
    MPI_Barrier(MPI_COMM_WORLD);
    start = clock();
    double V_variable[DIM][DIM];
    for (i = 0;i < DIM;i++)
    {
        for (j = 0;j < DIM;j++)
        {
            if (i == j)
                V_variable[i][j] = 1.0;
            else
                V_variable[i][j] = 0.0;
        }
    }
    MPI_Status status_v[2], status_seq[4], status_rseq[4],status_cseq[2],status_rcseq[2];
    MPI_Request request_v[4];
    MPI_Request request_rv[4];
    MPI_Request request_c[2], request_rc[2];
    int flag_seq[4], flag_rseq[4];
    int* flag_seq_z = flag_seq;
    int* flag_rseq_z = flag_rseq;
    double alpha = 0.1854;
    int iterations = 390;
    /*         ������ʼ     ͨ��Ϊ�������ͣ���������  */
    for (int ite = 0;ite < iterations;ite++)
    {
        if ((myrank == 4) || (myrank == 5))
        {
            // printf("%d process \n", myrank);
            double send_v[DIM] = { 0.0 };
            double send_v12[DIM] = { 0.0 };
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                if (mm != -1)
                {
                    for (j = 0;j < DIM;j++)
                    {
                        send_v[j] = send_v[j] + C_local[0][2][mm] * V_variable[j][mm];;
                        send_v12[j] = send_v12[j] + C_local[0][12][mm] * V_variable[j][mm];;
                    }

                }
            }
            int tag = 0;
            if (myrank == 4)
                tag = 42;
            if (myrank == 5)
                tag = 52;
             MPI_Ssend(send_v, DIM, MPI_DOUBLE, 2, tag, MPI_COMM_WORLD);///���ϴ���Ϣ
          
            int tag2 = 0;
            if (myrank == 4)
                tag2 = 43;
            if (myrank == 5)
                tag2 = 53;
            MPI_Ssend(send_v12, DIM, MPI_DOUBLE, 2, tag2, MPI_COMM_WORLD); 
       
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                double g_v[DIM] = { 0.0 };
                if ((mm != -1) && (mm != 2)&&(mm!=12))
                {
                    for (j = 0;j < 6;j++)
                    {
                        int linen = m_local[0][j];
                        if (linen != -1)
                        {
                            for (int jj = 0;jj < DIM;jj++)
                            {
                                g_v[jj] = g_v[jj] + C_local[0][mm][linen] * V_variable[jj][linen];
                            }
                        }
                    }
                    for (j = 0;j < DIM;j++)
                    {
                        g_v[j] = V_variable[j][mm] - alpha * g_v[j];
                    }
                    double* vvt = g_v;
                    for (j = 0;j < DIM;j++)
                    {
                        V_variable[j][mm] = g_v[j] * 1.0 / norm1_f(vvt, DIM);
                    }
                
                }
            }
            double v_recv[DIM];
            double v_recv2[DIM];
            int tag_r = 0;
            if (myrank == 4)
                tag_r = 24;
            if (myrank == 5)
                tag_r = 25;
             MPI_Recv(v_recv, DIM, MPI_DOUBLE, 2, tag_r, MPI_COMM_WORLD, &status_v[0]);
         
            for (j = 0;j < DIM;j++)
            {
                V_variable[j][2] = v_recv[j];
            }
            int tag_r2 = 0;
            if (myrank == 4)
                tag_r2 = 34;
            if (myrank == 5)
                tag_r2 = 35;
            MPI_Recv(v_recv2, DIM, MPI_DOUBLE, 2, tag_r2, MPI_COMM_WORLD, &status_v[1]);
           
            for (j = 0;j < DIM;j++)
            {
                V_variable[j][12] = v_recv2[j];
            }
        

        }
        if (myrank == 2)
        {
            double send_v[DIM] = { 0.0 };
            double send_v2[DIM] = { 0.0 };
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                if (mm != -1)
                {
                    for (j = 0;j < DIM;j++)
                    {
                        send_v[j] = send_v[j] + C_local[0][0][mm] * V_variable[j][mm];
                    }

                }
                mm = m_local[1][i];
                if (mm != -1)
                {
                    for (j = 0;j < DIM;j++)
                    {
                        send_v2[j] = send_v2[j] + C_local[1][3][mm] * V_variable[j][mm];;
                    }

                }
            }
            int tag = 21;

             MPI_Ssend(send_v, DIM, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);///���ϴ���Ϣ
             MPI_Ssend(send_v2, DIM, MPI_DOUBLE, 1, 22, MPI_COMM_WORLD);///���ϴ���Ϣ*/
         
            double v_recv[DIM] = { 0.0 };
            double v_recv2[DIM] = { 0.0 };
            MPI_Recv(v_recv, DIM, MPI_DOUBLE, 4, 42, MPI_COMM_WORLD, &status_v[0]);
            MPI_Recv(v_recv2, DIM, MPI_DOUBLE, 5, 52, MPI_COMM_WORLD, &status_v[1]);
       
            MPI_Status status_cv[2];
            double v_rec[DIM] = { 0.0 };
            double v_rec2[DIM] = { 0.0 };
            MPI_Recv(v_rec, DIM, MPI_DOUBLE, 4, 43, MPI_COMM_WORLD, &status_cv[0]);
            MPI_Recv(v_rec2, DIM, MPI_DOUBLE, 5, 53, MPI_COMM_WORLD, &status_cv[1]);
        
            double send_rl[DIM] = { 0.0 };
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                double g_v[DIM] = { 0 };
                if ((mm != -1) && (mm != 0) && (mm != 3))
                {
                    for (j = 0;j < 6;j++)
                    {
                        int linen = m_local[0][j];
                        int linen2 = m_local[1][j];
                        if (mm == 2)
                        {
                            if (linen != -1)
                            {
                                for (int jj = 0;jj < DIM;jj++)
                                {
                                    g_v[jj] = g_v[jj] + C_local[0][mm][linen] * V_variable[jj][linen];
                                }
                            }
                            if (linen2 != -1)
                            {
                                for (int jj = 0;jj < DIM;jj++)
                                {

                                    g_v[jj] = g_v[jj] + C_local[1][mm][linen2] * V_variable[jj][linen2];
                                }
                            }
                        }
                        if (mm == 12)///////������������
                        {
                            if (linen != -1)
                            {
                                for (int jj = 0;jj < DIM;jj++)
                                {
                                    g_v[jj] = g_v[jj] + C_local[0][mm][linen] * V_variable[jj][linen];
                                }
                            }
                        }
                        if (mm == 15)
                        {
                            if (linen2 != -1)
                            {
                                for (int jj = 0;jj < DIM;jj++)
                                {

                                    g_v[jj] = g_v[jj] + C_local[1][mm][linen2] * V_variable[jj][linen2];
                                }
                            }
                        }
                    }
                    if (mm == 2)
                    {
                        for (j = 0;j < DIM;j++)
                        {
                            g_v[j] = g_v[j] + v_recv[j] + v_recv2[j];
                        }
                    }
                    if (mm == 12)
                    {
                        for (j = 0;j < DIM;j++)
                        {
                            g_v[j] = g_v[j] + v_rec[j] + v_rec2[j];
                        }
                    }
                    for (j = 0;j < DIM;j++)
                    {
                        g_v[j] = V_variable[j][mm] - alpha * g_v[j];
                    }
                    double* vvt = g_v;
                    for (j = 0;j < DIM;j++)
                    {
                        V_variable[j][mm] = g_v[j] * 1.0 / norm1_f(vvt, DIM);
                       // send_rl[j] = V_variable[j][mm];
                    }
                    if (mm == 2)
                    {
                        for (j = 0;j < DIM;j++)
                        {
                            send_rl[j] = V_variable[j][mm];
                        }
                    
                        MPI_Ssend(send_rl, DIM, MPI_DOUBLE, 4, 24, MPI_COMM_WORLD);///���´���Ϣ
                        MPI_Ssend(send_rl, DIM, MPI_DOUBLE, 5, 25, MPI_COMM_WORLD);///���´���Ϣ
                        
                    }
                    if (mm == 12)
                    {
                        for (j = 0;j < DIM;j++)
                        {
                            send_rl[j] = V_variable[j][mm];
                        }
       
                        MPI_Ssend(send_rl, DIM, MPI_DOUBLE, 4, 34, MPI_COMM_WORLD);///���´���Ϣ
                        MPI_Ssend(send_rl, DIM, MPI_DOUBLE, 5, 35, MPI_COMM_WORLD);///���´���Ϣ

                    }
                }
            }
            MPI_Status status_v2[2];
            double v_recvp[DIM];
            double v_recvp2[DIM];
            MPI_Recv(v_recvp, DIM, MPI_DOUBLE, 1, 11, MPI_COMM_WORLD, &status_v2[0]);//�������Ը�����Ϣ
            MPI_Recv(v_recvp2, DIM, MPI_DOUBLE, 1, 12, MPI_COMM_WORLD, &status_v2[1]);
       
            for (j = 0;j < DIM;j++)
            {
                V_variable[j][0] = v_recvp[j];
                V_variable[j][3] = v_recvp2[j];
            }
      

        }

        if (myrank == 3)
        {
            double send_v[DIM] = { 0.0 };
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                if (mm != -1)
                {
                    for (j = 0;j < DIM;j++)
                    {
                        send_v[j] = send_v[j] + C_local[0][3][mm] * V_variable[j][mm];
                    }

                }
            }

             MPI_Ssend(send_v, DIM, MPI_DOUBLE, 1, 31, MPI_COMM_WORLD);///���ϴ���Ϣ
          
            double send_v2[DIM] = { 0.0 };
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                if (mm != -1)
                {
                    for (j = 0;j < DIM;j++)
                    {
                        send_v2[j] = send_v2[j] + C_local[0][13][mm] * V_variable[j][mm];
                    }

                }
            }

             MPI_Ssend(send_v2, DIM, MPI_DOUBLE, 1, 32, MPI_COMM_WORLD);///���ϴ���Ϣ
         
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                double g_v[DIM] = { 0.0 };
                if ((mm != -1) && (mm != 3)&&(mm!=13))
                {
                    for (j = 0;j < 6;j++)
                    {
                        int linen = m_local[0][j];
                        if (linen != -1)
                        {
                            for (int jj = 0;jj < DIM;jj++)
                            {
                                g_v[jj] = g_v[jj] + C_local[0][mm][linen] * V_variable[jj][linen];
                            }
                        }
                    }
                    for (j = 0;j < DIM;j++)
                    {
                        g_v[j] = V_variable[j][mm] - alpha * g_v[j];
                    }
                    double* vvt = g_v;
                    for (j = 0;j < DIM;j++)
                    {
                        V_variable[j][mm] = g_v[j] * 1.0 / norm1_f(vvt, DIM);
                    }
            
                }
            }
            double v_recv[DIM];
            int tag_r = 13;
            MPI_Recv(v_recv, DIM, MPI_DOUBLE, 1, tag_r, MPI_COMM_WORLD, &status_v[0]);//�������Ը���Ϣ
        
            for (j = 0;j < DIM;j++)
            {
                V_variable[j][3] = v_recv[j];
            }
            double v_recv2[DIM];
             tag_r = 14;
            MPI_Recv(v_recv2, DIM, MPI_DOUBLE, 1, 14, MPI_COMM_WORLD, &status_v[1]);//�������Ը���Ϣ
        
            for (j = 0;j < DIM;j++)
            {
                V_variable[j][13] = v_recv2[j];
            }
        
        }
        if (myrank == 1)
        {
            double v_recv[DIM];
            double v_recv2[DIM];
            double v_recv3[DIM];
            double v_recv4[DIM];
            MPI_Status status_v3[2];
             MPI_Recv(v_recv, DIM, MPI_DOUBLE, 2, 21, MPI_COMM_WORLD, &status_v[0]);
             MPI_Recv(v_recv2, DIM, MPI_DOUBLE, 2, 22, MPI_COMM_WORLD, &status_v[1]);
             MPI_Recv(v_recv3, DIM, MPI_DOUBLE, 3, 31, MPI_COMM_WORLD, &status_v3[0]);
             MPI_Recv(v_recv4, DIM, MPI_DOUBLE, 3, 32, MPI_COMM_WORLD, &status_v3[1]);
       
            double send_rl[DIM] = { 0 };
            double send_rl4[DIM] = { 0 };
            for (i = 0;i < 6;i++)
            {
                int mm = m_local[0][i];
                double g_v[DIM] = { 0 };
                if ((mm != -1))
                {
                    for (j = 0;j < 6;j++)
                    {
                        int linen = m_local[0][j];
                        if (linen != -1)
                        {
                            for (int jj = 0;jj < DIM;jj++)
                            {
                                g_v[jj] = g_v[jj] + C_local[0][mm][linen] * V_variable[jj][linen];

                            }
                        }
                    }
                    if (mm == 3)
                    {
                        for (j = 0;j < DIM;j++)
                        {
                            g_v[j] = g_v[j] + v_recv2[j] + v_recv3[j];
                        }
                    }
                    if (mm == 0)
                    {
                        for (j = 0;j < DIM;j++)
                        {
                            g_v[j] = g_v[j] + v_recv[j];
                        }
                    }
                    if (mm == 13)
                    {
                        for (j = 0;j < DIM;j++)
                        {
                            g_v[j] = g_v[j] + v_recv4[j];
                        }
                    }
                    for (j = 0;j < DIM;j++)
                    {
                        g_v[j] = V_variable[j][mm] - alpha * g_v[j];
                    }
                    double* vvt = g_v;
                    for (j = 0;j < DIM;j++)
                    {
                        V_variable[j][mm] = g_v[j] * 1.0 / norm1_f(vvt, DIM);

                    }
                
                    if (mm == 0)
                    {
                        for (j = 0;j < DIM;j++)
                            send_rl[j] = V_variable[j][mm];
                        MPI_Ssend(send_rl, DIM, MPI_DOUBLE, 2, 11, MPI_COMM_WORLD);///���´���Ϣ
                      
                    }
                    if (mm == 3)
                    {
                        for (j = 0;j < DIM;j++)
                            send_rl4[j] = V_variable[j][mm];
                         MPI_Ssend(send_rl4, DIM, MPI_DOUBLE, 2, 12, MPI_COMM_WORLD);///���´���Ϣ
                         MPI_Ssend(send_rl4, DIM, MPI_DOUBLE, 3, 13, MPI_COMM_WORLD);///����*/
                     
                    }
                    double sendrl2[DIM];
                    if (mm == 13)
                    {
                        for (j = 0;j < DIM;j++)
                            sendrl2[j] = V_variable[j][mm];
                        MPI_Ssend(sendrl2, DIM, MPI_DOUBLE, 3, 14, MPI_COMM_WORLD);///���´���Ϣ
                     
                    }
                }

            }
        }
        if (ite == iterations - 1)
        {
 
            double Vseq[6][DIM][DIM];
            MPI_Gather(V_variable, DIM * DIM, MPI_DOUBLE, Vseq, DIM * DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            double V_end[DIM][DIM];
            if (myrank == 0)
            {

                for (i = 0;i < DIM;i++)
                {
                    V_end[i][0] = Vseq[1][i][0];
                    V_end[i][1] = Vseq[1][i][1];
                    V_end[i][3] = Vseq[1][i][3];
                    V_end[i][2] = Vseq[2][i][2];
                    V_end[i][4] = Vseq[3][i][4];
                    V_end[i][5] = Vseq[4][i][5];
                    V_end[i][6] = Vseq[4][i][6];
                    V_end[i][7] = Vseq[5][i][7];
                    V_end[i][8] = Vseq[5][i][8];
                    V_end[i][9] = Vseq[1][i][9];
                    V_end[i][10] = Vseq[4][i][10];
                    V_end[i][11] = Vseq[5][i][11];
                    V_end[i][12] = Vseq[2][i][12];
                    V_end[i][13] = Vseq[1][i][13];
                    V_end[i][14] = Vseq[4][i][14];
                    V_end[i][15] = Vseq[2][i][15];
                    V_end[i][16] = Vseq[5][i][16];
                    V_end[i][17] = Vseq[3][i][17];
                }
               
                double Xm[DIM][DIM];
                for (i = 0;i < DIM;i++) {
                    for (j = 0;j < DIM;j++) {
                        double sumn = 0.0;
                        for (int k = 0;k < DIM;k++) {

                            sumn = sumn + V_end[k][i] * V_end[k][j];////* (*(result_matrix + i) + j)
                        }
                        Xm[i][j] = sumn;
                    }
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
                double flim = -36.98024;
                printf("f function value is %f!!!!!!!!!!!\n", fvalue);
                printf("f-f* function value is %f!!!!!!!!!!!\n", fvalue - flim);
            }
        }
    
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    finish = clock();
    Total_time = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("\n��������ʱ�䣺%0.3f���� \n", Total_time); //��ӡС����ĺ���λ������Ϊ��λ

   
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