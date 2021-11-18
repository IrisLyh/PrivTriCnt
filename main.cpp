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
    // 初始化所有的参数，以及所需要的数据
    // srand((unsigned)time(NULL));
    // string path ="./eps.txt";
    // ofstream file(path.data());

    // double eps_2 = 2.4;
    // double global_delta = 1e-4;
    // for(int k=1;k<=5000;k++)
    // {
    //     double eps_prime = get_epsilon(eps_2,0.01,0,10,global_delta,k,global_delta);
    //     file<<eps_prime<<endl;
    // }
    // file.close();
    //循环所有数据集
    bool ye = false;
    for(int datasetNo=5;datasetNo<=5;datasetNo++)
    {
        Graph dataset(datasetNo);
        Metric metric(dataset.getDatasetName());
        int userNum = dataset.getUserNum();
        double *senCnt = new double[userNum];
        double **mat = new double*[userNum];
        double global_delta = 1e-4;
         
        for(int i=0;i<userNum;i++)
        {
            mat[i] = new double[userNum];
        }
        mat = dataset.getMatrix();
        int total_edge = 0;
        for(int i=0;i<userNum;i++)
        {
            for(int j=0;j<userNum;j++)
            {
                if(mat[i][j]==1)
                {
                    total_edge++;
                }
            }
        }
        cout<<"total edge:"<<total_edge<<endl;
    }
        /*
        double *degreeList = new double[userNum];
        double *triCnt = new double[userNum];
        degreeList = dataset.getDegreeList();
        string path = "./oriDegree.txt";
        bool flag = false;
        if(flag)
        {
            ofstream outfile_tri(path.data());
            for(int i=0;i<userNum;i++)
            {
                outfile_tri<<degreeList[i]<<endl;
            }
            outfile_tri.close();
        }
        cout<<"counting triangle..."<<endl;
        triCnt = metric.TriangleCount(mat,userNum);
        cout<<"counting sensitivity..."<<endl;
        senCnt = metric.sensitivityCount(mat,userNum);
        double global_sensitivity = 3*(userNum-2); 

        double oriTotal=0;
        double beta = 0.2;
        double sample_rate = 1;
        double oriSenMax = 0;
        for(int i=0;i<userNum;i++)
        {
            if(senCnt[i]>oriSenMax)
            {
                oriSenMax = senCnt[i];
            }
        }
        // cout<<endl<<"***************"<<dataset.getDatasetName()<<"***************"<<endl;
        // cout<<"Max:"<<oriSenMax<<endl;
        // for(int eps=1;eps<=10;eps++)
        // {
        //     sample_rate=0.01;
        //     double target_eps = (double)eps;
        //     double noMax = getNoisyMax(target_eps*beta, global_delta, degreeList, senCnt, userNum);
        //     double user_eps = get_epsilon(target_eps*(1-beta), sample_rate, 0, 10, global_delta, noMax, global_delta);
        //     cout<<"eps:"<<eps<<"\t eps_1:"<<target_eps*beta<<"\t noMax:"<<noMax<<"\t user_eps:"<<user_eps<<endl;
        // }


            // my method
            // subample
        //     double *subsample_count = new double[userNum];
        //     double p=sample_rate;
        //     for(int i=0;i<userNum;i++)
        //     {
        //         subsample_count[i] = 0;
        //         for(int j=0;j<triCnt[i];j++)
        //         {
        //             double u = (double)rand()/(double)RAND_MAX;
        //             if(u<=p)
        //             {
        //                 subsample_count[i]++;
        //             }
        //         }
        //         // cout<<"subsample count:"<<subsample_count[i]<<endl;
        //     }

        //     //add noise to subample result
        //     //compute scale for each user
        //     double* my_scale=new double[userNum];
        //     for(int i=0;i<userNum;i++)
        //     {
        //         my_scale[i] = degreeList[i];
        //     }
        //     my_scale = perturbScale(target_eps*beta/2, global_delta, my_scale,userNum);
        //     char buffer[20];
        //     sprintf(buffer,"%d",eps);
        //     string strEps = buffer;
        //     string res_path = "./upper_sample_"+dataset.getDatasetName()+"_"+strEps+"_degree.txt";
        //     ofstream outfile_res(res_path.data(),ios::app);
        //     for(int userId=0;userId<userNum;userId++)
        //     {
        //         outfile_res<<my_scale[userId]<<endl;
        //     }
        //     outfile_res.close();
        // }
 
        for(int i=0;i<userNum;i++)
        {
            oriTotal += triCnt[i];
        }
        if(!ye)
        {
            int rateLen = 8;
            double* sampleRate = new double[rateLen];
            string rateFilePath = "./rate.txt";
            ifstream rateFile(rateFilePath.data()); 
            for(int rateNo=0;rateNo<rateLen;rateNo++)
            {
                rateFile>>sampleRate[rateNo];
            }
            //循环所有sample rate
            for(int rateTrans=0;rateTrans<rateLen;rateTrans++)
            {
                sample_rate = sampleRate[rateTrans];
                double target_eps = 0.5;
                //循环所有epsilon
                for(int eps_trans = 0;eps_trans<19;eps_trans++)
                {
                    // double *outScale = new double[userNum]; //记录每个user在不同的参数下的sensitivity
                    // double noisyNumber = 0; //记录每个数据集noisy的边最大重复使用次数
                    double my_MSE=0;
                    double baseline_MSE = 0;
                    double my_baseline_MSE = 0;
                    double DDP_MSE = 0;
                    double sample_MSE = 0;
                    target_eps += 0.5;
                    // 多次实验只获取一次noisy max
                    cout<<"[INFO] Total triangle count:"<<oriTotal<<endl;
                    // pre-try for subsampling and perturb
                    int trials = 1000;
                    // 实验做trials这么多次取最后结果的平均值
                    for(int times=0;times<trials;times++)
                    {
                        // 单次变量初始化
                        double noMax = getNoisyMax(target_eps*beta, global_delta, degreeList, senCnt, userNum);
                        double user_eps = get_epsilon(target_eps*(1-beta), sample_rate, 0, 10, global_delta, noMax, global_delta);
                        

                        double my_count=0;
                        double my_subsample_count = 0;
                        double baseline_count = 0;
                        double DDP_count = 0;
                        double my_baseline_count = 0;

                        double my_eps = user_eps;
                        double baseline_eps = target_eps;
                        double DDP_eps = target_eps*(1-beta);
                        double my_baseline_eps = target_eps*(1-beta)/noMax;

                        double my_subsample_error = 0;
                        double my_error = 0;
                        double baseline_error = 0;
                        double DDP_error = 0;
                        double my_baseline_error = 0;

                        
                        // my method
                        // subample
                        double *subsample_count = new double[userNum];
                        double p=sample_rate;
                        for(int i=0;i<userNum;i++)
                        {
                            subsample_count[i] = 0;
                            for(int j=0;j<triCnt[i];j++)
                            {
                                double u = (double)rand()/(double)RAND_MAX;
                                if(u<=p)
                                {
                                    subsample_count[i]++;
                                }
                            }
                            // cout<<"subsample count:"<<subsample_count[i]<<endl;
                        }
                        // 测试subsample的error
                        for(int i=0;i<userNum;i++)
                        {
                            my_subsample_count += subsample_count[i];
                        }
                        my_subsample_count = my_subsample_count / p;
                        // cout<<"count:"<<my_subsample_count<<endl;
                        my_subsample_error=abs(my_subsample_count-oriTotal)/my_subsample_count*100;
                        sample_MSE += my_subsample_error;
                        // cout<<"[INFO]Subsampling error:"<<my_subsample_error<<"%"<<endl;
                        
                        
                        //add noise to subample result
                        //compute scale for each user
                        double* my_scale=new double[userNum];
                        double* my_baseline_scale = new double[userNum];
                        for(int i=0;i<userNum;i++)
                        {
                            my_scale[i] = degreeList[i];
                        }
                        my_scale = perturbScale(target_eps*beta/2, global_delta, my_scale,userNum);


                        for(int i=0;i<userNum;i++)
                        {
                            my_baseline_scale[i] = my_scale[i];
                            my_scale[i] = my_scale[i] * p + 4*sqrt(p*(1-p)*my_scale[i]);
                            if(my_scale[i]<1)
                            {
                                my_scale[i]=1;
                            }
                            if(my_baseline_scale[i]<1)
                            {
                                my_baseline_scale[i]=1;
                            }
                            my_scale[i] = my_scale[i] / my_eps;
                            my_baseline_scale[i] = my_baseline_scale[i] / my_baseline_eps;
                            double u = rand()/(double)(RAND_MAX)-0.5;
                            double noise = (-1)*my_scale[i]*(double)main_sign(u)*log(1-2*abs(u));
                            subsample_count[i] = subsample_count[i] + noise;
                            my_count += subsample_count[i];

                            double baseline_noise = (-1)*(3*global_sensitivity/baseline_eps)*(double)main_sign(u)*log(1-2*abs(u));
                            double tmp_cnt = triCnt[i] + baseline_noise;
                            baseline_count += tmp_cnt;


                            double DDP_noise = (-1)*(3*noMax/DDP_eps)*(double)main_sign(u)*log(1-2*abs(u));
                            tmp_cnt = triCnt[i] + DDP_noise;
                            DDP_count += tmp_cnt;

                            double my_baseline_noise = (-1)*(my_baseline_scale[i])*(double)main_sign(u)*log(1-2*abs(u));
                            tmp_cnt = triCnt[i] + my_baseline_noise;
                            my_baseline_count += tmp_cnt;
                            
                        }
                        my_count = my_count / p;
                        my_error= abs(my_count - oriTotal)/my_count*100;
                        my_MSE += my_error;

                        baseline_error = abs(baseline_count - oriTotal)/oriTotal*100;
                        baseline_MSE += baseline_error;

                        DDP_error = abs(DDP_count-oriTotal)/oriTotal*100; 
                        DDP_MSE += DDP_error;

                        my_baseline_error = abs(my_baseline_count - oriTotal)/oriTotal*100;
                        my_baseline_MSE += my_baseline_error;
                    }
                    // for(int ses=0;ses<userNum;ses++)
                    // {
                    //     outScale[ses] /= trials;
                    // }
                    sample_MSE = sample_MSE/trials;
                    cout<<"sample rate:"<<sample_rate<<endl;
                    cout<<"sample error:"<<sample_MSE<<"%"<<endl;
                    cout<<"my MRE error:"<<my_MSE/trials<<"%"<<endl;
                    cout<<"baseline error:"<<baseline_MSE/trials<<"%"<<endl;
                    cout<<"DDP error:"<<DDP_MSE/trials<<"%"<<endl;
                    cout<<"my baseline error:"<<my_baseline_MSE/trials<<"%"<<endl;
                    // cout<<"noisy max:"<<noMax<<endl;
                    char buffer[20];
                    sprintf(buffer,"%f",sample_rate);
                    string strRate = buffer;
                    string res_path = "./result_"+dataset.getDatasetName()+"_"+strRate+"_degree.txt";
                    ofstream outfile_res(res_path.data(),ios::app);
                    outfile_res<<target_eps<<":"<<my_MSE/trials/100<<";"<<DDP_MSE/trials/100<<";"<<my_baseline_MSE/trials/100<<";"<<sample_MSE/100<<endl;
                    outfile_res.close();
                    // string compare_file_path = "./" + dataset.getDatasetName() + "_" + strRate + "_noisySensitivity.txt";
                    // ofstream compare_file(compare_file_path.data(),ios::app);
                    // for(int ses=0;ses<userNum;ses++)
                    // {
                    //     compare_file<<target_eps<<":"<<senCnt[ses]<<";"<<degreeList[ses]<<";"<<outScale[ses]<<endl;
                    // }
                    // compare_file.close();
                }
            }
        }
        else
        {
            double target_eps = 0;
            //循环所有epsilon
            string datasetName = dataset.getDatasetName();
            string ye_path ="./"+ datasetName +"_ye_result.txt";
            ofstream ye_file(ye_path.data(),ios::app);
            for(int eps_trans = 0;eps_trans<10;eps_trans++)
            {
                target_eps += 1;
                lf_gdpr lf(datasetNo,mat,target_eps,degreeList);
                double noiTotal = lf.getTCnt();
                double error = abs(3*noiTotal-oriTotal)/oriTotal*100;
                ye_file<<target_eps<<":"<<error<<"%"<<endl;
            }
            ye_file.close();
        }
    }*/
    return 0;
}