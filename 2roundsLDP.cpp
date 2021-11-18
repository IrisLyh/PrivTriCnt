#ifndef _BASESMOOTHINGALGORITHM_H_
#define _BASESMOOTHINGALGORITHM_H_
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <typeinfo>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>
#include "data.cpp"
#include "metric.cpp"
#include "utils.cpp"
#include "lf-gdpr.cpp"
#endif

using namespace std;
int main()
{
    bool test_performance = true;
    int datasetNo=5;
    int totalStartTime=clock();
    //重新计算triCnt
    Graph dataset(datasetNo);
    Metric metric(dataset.getDatasetName());
    int userNum = dataset.getUserNum();
    double **oriMatrix = new double*[userNum];
    double *triCnt = new double[userNum];
    double *pathCnt = new double[userNum];
    for(int i=0;i<userNum;i++)
    {
        oriMatrix[i] = new double[userNum];
    }
    oriMatrix = dataset.getMatrix();
    double **ptrMatrix = new double* [userNum];
    for(int i=0;i<userNum;i++)
    {
        ptrMatrix[i] = new double[userNum];
    }
    double epsilon = 1;
    double p = exp(epsilon)/(1+exp(epsilon));
    srand((unsigned)time(NULL));
    double *ptrTime = new double[userNum];
    for(int i=0;i<userNum;i++)
    {
        int ptrStartTime = clock();
        for(int j=i+1;j<userNum;j++)
        {
            double u = (double)rand()/(double)RAND_MAX;
            if(u<=p)
            {
                ptrMatrix[i][j] = oriMatrix[i][j];
                ptrMatrix[j][i] = oriMatrix[i][j];
            }
            else
            {
                ptrMatrix[i][j] = 1-oriMatrix[i][j];
                ptrMatrix[j][i] = 1-oriMatrix[i][j];
            }
        }
        int ptrEndTime = clock();
        ptrTime[i] = ptrEndTime - ptrStartTime;
        cout<<"[INFO]nodeID:"<<i<<"\t perturb time:"<<ptrTime[i]<<endl;
    }
    double* triCntTime = new double[userNum];
    for(int i=0;i<userNum;i++)
    {
        int triCntStartTime = clock();
        for(int j=i+1;j<userNum;j++)
        {
            for(int k=j+1;k<userNum;k++)
            {
                if(oriMatrix[i][j]==1 && oriMatrix[i][k]==1 && ptrMatrix[j][k]==1)
                {
                    triCnt[i] = triCnt[i]+1;
                }
            }
        }
        int triCntEndTime = clock();
        triCntTime[i]=triCntEndTime-triCntStartTime;
        cout<<"[INFO]nodeID:"<<i<<"\t count time:"<<triCntTime[i]<<endl;
    }
    int totalEndTime = clock();
    int totalRunTime = totalEndTime-totalStartTime;
    int* userTime = new int[userNum];
    int totalUserTime = 0;
    for(int i=0;i<userNum;i++)
    {
        userTime[i] = ptrTime[i]+triCntTime[i];
        totalUserTime = totalUserTime + userTime[i];
    }
    double avgTime = (double)totalUserTime/(double)userNum;
    cout<<"**************************AVG TIME*****************************"<<endl;
    cout<<"dataset:"<<dataset.getDatasetName()<<endl;
    cout<<"average time per user for 2Rounds LDP:"<<avgTime<<endl;
    cout<<"total running time:"<<totalRunTime<<endl;
    return 0;
}

        
        