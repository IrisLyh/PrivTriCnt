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
double epsilon_func(double x, double p, double noMax, double delta_prime)
{
    // input x, return y
    double y = 2*x+log(1+p*(exp(x)-1))*sqrt(2*noMax*log(1/delta_prime))+noMax*log(1+p*(exp(x)-1))*p*(exp(x)-1)/(p*(exp(x)-1)+2);
    return y;
}
double get_epsilon(double target_eps, double p, double start, double stop, double threshhold, double noMax, double delta_prime)
{
    //利用二分法找到每个user实际上应该添加多少的噪声
    // p is sample rate
    double res=(double)RAND_MAX;
    double lower = start;
    double upper = stop;
    double middle = (start + stop)/2;
    double low_eps = epsilon_func(lower, p, noMax,delta_prime);
    double up_eps = epsilon_func(upper, p, noMax,delta_prime);
    double middle_eps = 0;
    if(low_eps>target_eps)
    {
        cout<<"ERROR! lower eps is too high, please select another parameter!"<<endl;
        return -1;
    }
    if(up_eps<target_eps)
    {
        cout<<"ERROR! upper eps is too low, please select another parameter!"<<endl;
        return -1;
    }
    while(res>threshhold)
    {
        middle_eps = epsilon_func(middle, p, noMax, delta_prime);
        if(middle_eps>target_eps)
        {
            upper = middle;
        }
        else if (middle_eps<target_eps)
        {
            lower = middle;
        }
        else
        {
            return middle;
        }
        res = abs(middle_eps - target_eps);
        middle = (lower + upper)/2;
    }
    return middle;
}


int main_sign(double num)
{
    if(num>0)
    {
        return 1;
    }
    else if (num<0)
    {
        return -1;
    }
    else if(num==0)
    {
        return 0;
    }
    else
    {
        cout<<"[ERROR]Wrong num for SIGN function!!"<<endl;
        exit(0);
    }
}

bool cmp(pair<int, double>a, pair<int, double>b)
{
    return a.second>b.second;
}

// 这个阶段大家只考虑自己就可以了
// 之后组合的时候整体组合就可以了
// 目前的疑问是，这里是不是受subsample的影响
// 参数传原指针的话在函数中就不能改变原指针的值了，否则会影响最终的结果
double* perturbScale(double eps, double delta, double* scale, int userNum)
{
    
    LaplaceNoise lp_scale(userNum,1,2,eps);
    double** noise = new double*[userNum];
    for(int i=0;i<userNum;i++)
    {
        noise[i] = new double[1];
    }
    noise = lp_scale.getNoise();

    for(int i=0;i<userNum;i++)
    {
        scale[i] = scale[i] + noise[0][i] + (2.0/eps)*log(1/(2*delta));
        if(scale[i]<1)
        {
            scale[i] = 1;
        }
    }
    return scale;
}
// adaptive optimized soltion from CCS'19
double getNoisyMax(double eps, double delta, double* degreeList, double* sensitivityList, int userNum)
{   
    // cout<<"in noisy max!"<<endl;
    double alpha = 0.5;   
    double eps1 = eps * alpha;
    double eps2 = eps * (1 - alpha);
    double* noisyDegree = new double[userNum];
    double* noisyUpperBound = new double[userNum];
    double h_prime = 100.0;
    double h=0;
    double delta_prime = delta / 2*(h_prime+1);
    double max=0;
    srand((unsigned)time(NULL));
    for(int i=0;i<userNum;i++)
    {
        double u =  rand()/(double)(RAND_MAX)-0.5;
        double noise = (-1)*(2.0/eps1)*(double)main_sign(u)*log(1-2*abs(u));
        noisyDegree[i] = 0;
        noisyDegree[i] = degreeList[i] + noise + (2/eps1)*log(1/(2*delta_prime));
    }
    pair<int,double> degreeUpperBound[userNum];
    vector<pair<int, double> >vec;
    for(int i=0;i<userNum;i++)
    {
        degreeUpperBound[i] = make_pair(i,noisyDegree[i]);
        vec.push_back(degreeUpperBound[i]);
    }
    sort(vec.begin(), vec.end(), cmp);
    for(int i=0;i<int(h_prime);i++)
    {
        if(double(i+1)/eps1*log(1/(2*delta_prime)) >= vec[i+1+2].second)
        {
            h = ceil(((double)i+1)/2);
            break;
        }
    }
    if(h==0)
    {
        h=ceil(h_prime/2);
    }
    // cout<<"h:"<<h<<endl;
    for(int i=0;i<h;i++)
    {
        noisyUpperBound[i] = 0;
        noisyUpperBound[i] = sensitivityList[vec[i].first];
        double u =  rand()/(double)(RAND_MAX)-0.5;
        double noise = (-1)*(h/eps2)*(double)main_sign(u)*log(1-2*abs(u));
        noisyUpperBound[i] += noise + (h/eps2)*log(1/(2*delta_prime));
        noisyUpperBound[i] = min(noisyUpperBound[i],vec[i].second);
        if(noisyUpperBound[i]>max)
        {
            max = noisyUpperBound[i];
        }
    }
    return max;
}

int main()
{
    bool test_performance = true;
    int datasetNo=5;
    //重新计算triCnt
    Graph dataset(datasetNo);
    Metric metric(dataset.getDatasetName());
    int userNum = dataset.getUserNum();
    double **oriMatrix = new double*[userNum];
    for(int i=0;i<userNum;i++)
    {
        oriMatrix[i] = new double[userNum];
    }
    oriMatrix = dataset.getMatrix();
    double *senCnt = new double[userNum];
    double *triCnt = new double[userNum];
    double *subsample_count = new double[userNum];
    int *triTime = new int[userNum];
    int *sampleTime = new int[userNum];
    double *degreeList = new double[userNum];
    degreeList = dataset.getDegreeList();
    double *scale = new double[userNum];
    double p=0.01;//sample rate
    double epsilon = 1;//privacy budget
    double global_delta = 1e-4;//invalidation probability
    if(test_performance)
    {
        int TstartTime = clock();
        //计算triangle的时间
        for(int i=0;i<userNum;i++)
        {
            int starttime = clock();
            triCnt[i]=0;
            senCnt[i]=0;
            for(int j=0;j<userNum;j++)
            {
                double tempSenCnt = 0;
                for(int k=j+1;k<userNum;k++)
                {
                    if(i!=j && i!=k && j!=k)
                    {
                        if(oriMatrix[i][j]==1 && oriMatrix[i][k]==1 && oriMatrix[j][k]==1)
                        {
                            triCnt[i] = triCnt[i] + 1;
                        }
                        if(oriMatrix[i][j]==1 && oriMatrix[j][k]==1)
                        {
                            tempSenCnt = tempSenCnt + 1;
                        }
                    }
                }
                if(tempSenCnt>senCnt[i])
                {
                    senCnt[i]=tempSenCnt;
                }
            }
            int endtime = clock();
            triTime[i] = endtime-starttime;
            subsample_count[i] = 0;
            for(int j=0;j<triCnt[i];j++)
            {
                double u = (double)rand()/(double)RAND_MAX;
                if(u<=p)
                {
                    subsample_count[i]++;
                }
            }
            int sampleEndTime=clock();
            sampleTime[i] = sampleEndTime - endtime;
            cout<<"[INFO]nodeID:"<<i<<"\t counting time:"<<triTime[i]<<"\t subsample time:"<<sampleTime[i]<<endl;
        }
        int TendTime = clock();
        int total_time = 0;
        int total_time_DDP = 0;
        for(int user=0;user<userNum;user++)
        {
            total_time = total_time + triTime[user] + sampleTime[user];
            total_time_DDP = total_time_DDP + triTime[user];
        }
        double avgTime = (double)total_time/(double)userNum;
        double avgTime_DDP = (double)total_time_DDP/(double)userNum;
        int Ttime = TendTime - TstartTime;
        cout<<"**************************AVG TIME*****************************"<<endl;
        cout<<"dataset:"<<dataset.getDatasetName()<<endl;
        cout<<"average time per user for edge RLDP:"<<avgTime<<endl;
        cout<<"average time per user for edge DDP:"<<avgTime_DDP<<endl;
        cout<<"total running time:"<<Ttime<<endl;
    }
    return 0;
}