#ifndef _BASESMOOTHINGALGORITHM_H_
#define _BASESMOOTHINGALGORITHM_H_
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>
#endif
using namespace std;
int file_check(string filename)
{
    fstream filecheck;
    filecheck.open(filename.data(),ios::in);
    if(!filecheck)
    {
        cout<<filename<<" not exists!";
        return -1;
    }
    else
    {
        cout<<filename<<" exists! Reading file...";
        filecheck.close();
        return 1;
    }
}
class LaplaceNoise
{
    public:
        void setLength(int len) //设置单个向量长度
        {
            length = len;
        }
        void setWidth(int wid) //设置向量个数
        {
            width = wid;
        }
        void setScale(int sensitivity, double eps)//根据敏感度设置Laplace分布的参数,设置隐私预算
        {
            epsilon = eps;
            scale = (double)sensitivity/epsilon;
        }
        int sign(double num)
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
        double** setNoiseMatrix()
        {
            if(length<=0 || width<=0)
            {
                cout<<"[ERROR]Length and width must be larger that 0!!"<<endl;
                exit(0);
            }
            double** LaplaceDis = new double* [width];
            for(int i=0;i<width;i++)
            {
                LaplaceDis[i]=new double[length];
            }
            return LaplaceDis;
        }
        double ** getNoise()
        {
            double ** Noise = setNoiseMatrix();
            srand((unsigned)time(NULL));
            for(int i=0;i<width;i++)
            {
                for(int j=0;j<length;j++)
                {
                    double u = rand()/(double)(RAND_MAX)-0.5;
                    Noise[i][j] = (-1)*scale*(double)sign(u)*log(1-2*abs(u));
                }
            }
            return Noise;
        }
        LaplaceNoise(int len, int wid, int sensitivity, double eps)
        {
            setWidth(wid);
            setLength(len);
            setScale(sensitivity,eps);
        }
    protected:
        int length;
        int width;
        double scale;
        double epsilon;
};


// int main()
// {
//     ofstream outfile("noise_test.txt",ios::app);
//     LaplaceNoise LapDis = LaplaceNoise(10000,1,1,1);
//     double** noise = LapDis.getNoise();
//     return 0;
// }