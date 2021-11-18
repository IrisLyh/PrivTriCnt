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
#include "perturbation.cpp"
#include "utils.cpp"
#endif
using namespace std;

class Metric
{
    private:
        void setDatasetName(string dName)
        {
            dataset = dName;
        }
    public:
        double* LocalClusteringCoefficient(double* degreeList, double* triangle, int userNum)
        {
            double* LCC = new double[userNum];
            for(int i=0;i<userNum;i++)
            {
                if(degreeList[i]==1||degreeList[i]==0)
                {
                    LCC[i]=0;
                }
                else
                {
                    LCC[i] = 2.0*triangle[i]/(degreeList[i]*(degreeList[i]-1));
                }
            }
            return LCC;
        }
        double* TriangleCount(double** matrix, int userNum)
        {
            cout<<"[INFO] In this func for triangle counting!"<<endl;
            double *triCnt = new double[userNum];
            string path = "./"+dataset+"_oriCount.txt";
            if(file_check(path)==1)
            {
                ifstream infile_tri(path.data());
                for(int i=0;i<userNum;i++)
                {
                    infile_tri>>triCnt[i];
                }
                infile_tri.close();
            }
            else
            {
                cout<<"[INFO]File "+path+" not exists, creating..."<<endl;
                for(int i=0;i<userNum;i++)
                {
                    triCnt[i]=0;
                    for(int j=0;j<userNum;j++)
                    {
                        for(int k=j+1;k<userNum;k++)
                        {
                            if(j!=i && j!=k && i!=k)
                            {
                                if(matrix[i][j]==1 && matrix[i][k]==1 && matrix[j][k]==1)
                                {
                                    triCnt[i]++;
                                }
                            }
                        }
                    }
                    cout<<"[INFO]Triangel counting for node "<<i<<": "<<triCnt[i]<<endl;
                }
                ofstream outfile_tri(path.data());
                for(int i=0;i<userNum;i++)
                {
                    outfile_tri<<triCnt[i]<<endl;
                }
                outfile_tri.close();
            }
            return triCnt;
        }
        double* improvedTriCnt(double** matrix, int userNum)
        {
            cout<<"testMatrix:"<<matrix[0][1]<<endl;
            cout<<"[INFO] In this func for triangle counting!"<<endl;
            double *triCnt = new double[userNum];
            string path = "./"+dataset+"_impCount.txt";
            if(file_check(path)==1)
            {
                ifstream infile_tri(path.data());
                for(int i=0;i<userNum;i++)
                {
                    infile_tri>>triCnt[i];
                }
                infile_tri.close();
            }
            else
            {
                cout<<"[INFO]File "+path+" not exists, creating..."<<endl;
                for(int i=0;i<userNum;i++)
                {
                    triCnt[i]=0;
                    for(int j=i+1;j<userNum;j++)
                    {
                        for(int k=j+1;k<userNum;k++)
                        {
                            if(matrix[i][j]==1 && matrix[i][k]==1 && matrix[j][k]==1)
                            {
                                triCnt[i]=triCnt[i]+1;
                            }
                        }
                    }
                    cout<<"[INFO]Triangle number for node "<<i<<":"<<triCnt[i]<<endl;
                }
            }
            double total=0;
            for(int i=0;i<userNum;i++)
            {
                total = total+triCnt[i];
            }
            cout<<"Total:"<<total<<endl;
            ofstream outfile_tri(path.data());
            for(int i=0;i<userNum;i++)
            {
                outfile_tri<<triCnt[i]<<endl;
            }
            outfile_tri.close();
            return triCnt;
        }
        double* sensitivityCount(double** matrix, int userNum)
        {
            double *senCnt = new double[userNum];
            string path = "./"+dataset+"_oriSensitivity.txt";
            if(file_check(path)==1)
            {
                ifstream infile_sen(path.data());
                for(int i=0;i<userNum;i++)
                {
                    infile_sen>>senCnt[i];
                }
                infile_sen.close();
            }
            else
            {
                for(int i=0;i<userNum;i++)
                {
                    senCnt[i] = 0;
                    for(int j=0;j<userNum;j++)
                    {
                        int temp_sensitivity = 0;
                        for(int k=0;k<userNum;k++)
                        {
                            if(j!=i && j!=k && i!=k)
                            {
                                if(matrix[i][k] == 1 && matrix[j][k] == 1)
                                {
                                    temp_sensitivity++;
                                }
                            }
                        }
                        if(temp_sensitivity > senCnt[i])
                        {
                            senCnt[i] = temp_sensitivity;
                        }
                    }
                }
                ofstream outfile_sen(path.data());
                for(int i=0;i<userNum;i++)
                {
                    outfile_sen<<senCnt[i]<<endl;
                }
                outfile_sen.close();
            }
            return senCnt;
        }
        
        double* improvedSenCnt(double** matrix, int userNum)
        {
            double *senCnt = new double[userNum];
            string path = "./"+dataset+"_impSensitivity.txt";
            if(file_check(path)==1)
            {
                ifstream infile_sen(path.data());
                for(int i=0;i<userNum;i++)
                {
                    infile_sen>>senCnt[i];
                }
                infile_sen.close();
            }
            else
            {
                for(int i=0;i<userNum;i++)
                {
                    senCnt[i] = 0;
                    for(int j=i+1;j<userNum;j++)
                    {
                        int temp_sensitivity = 0;
                        for(int k=j+1;k<userNum;k++)
                        {
                            if(j!=i && j!=k && i!=k)
                            {
                                if(matrix[i][k] == 1 && matrix[j][k] == 1)
                                {
                                    temp_sensitivity++;
                                }
                            }
                        }
                        if(temp_sensitivity > senCnt[i])
                        {
                            senCnt[i] = temp_sensitivity;
                        }
                    }
                    cout<<"[INFO]Sensitivity number for node "<<i<<":"<<senCnt[i]<<endl;
                }
                ofstream outfile_sen(path.data());
                for(int i=0;i<userNum;i++)
                {
                    outfile_sen<<senCnt[i]<<endl;
                }
                outfile_sen.close();
            }
            return senCnt;
        }

        double* sensitivityAfterSample(double** matrix, int userNum, double sample_rate)
        {
            double *senCnt = new double[userNum];
            char buffer[20];
            sprintf(buffer,"%f",sample_rate);
            string strRate = buffer;
            string path = "./"+dataset+"_"+strRate+"_SampleSensitivity.txt";
            srand((unsigned)time(NULL));
            if(file_check(path)==1)
            {
                ifstream infile_sen(path.data());
                for(int i=0;i<userNum;i++)
                {
                    infile_sen>>senCnt[i];
                }
                infile_sen.close();
            }
            else
            {
                for(int i=0;i<userNum;i++)
                {
                    senCnt[i] = 0;
                    for(int j=0;j<userNum;j++)
                    {
                        int temp_sensitivity = 0;
                        int sample_sensitivity = 0;
                        for(int k=0;k<userNum;k++)
                        {
                            if(j!=i && j!=k && i!=k)
                            {
                                if(matrix[i][k] == 1 && matrix[j][k] == 1)
                                {
                                    temp_sensitivity++;
                                }
                            }
                        }
                        for(int rr=0;rr<temp_sensitivity;rr++)
                        {
                            double u = (double)rand()/(double)RAND_MAX;
                            if(u<=sample_rate)
                            {
                                sample_sensitivity++;
                            }
                        }
                        if(sample_sensitivity > senCnt[i])
                        {
                            senCnt[i] = sample_sensitivity;
                        }
                    }
                }
                ofstream outfile_sen(path.data());
                for(int i=0;i<userNum;i++)
                {
                    outfile_sen<<senCnt[i]<<endl;
                }
                outfile_sen.close();
            }
            return senCnt;
        }


        double* lapTriangleCnt(double** matrix, int userNum)
        {
            double* triCnt= new double [userNum];
            for(int i=0;i<userNum;i++)
            {
                triCnt[i]=0;
                for(int j=0;j<userNum;j++)
                {
                    for(int k=j+1;k<userNum;k++)
                    {
                        if(i!=j && i!=k && j!=k)
                        {
                            triCnt[i]=triCnt[i]+(matrix[i][j]+matrix[i][k]+matrix[j][k])/3;
                        }
                    }
                }
            }
            return triCnt;
        }
        double* lccCalibration(double* triCnt, double noisyDensity, double* noisyDegree, double eps, double alp, int userNum)
        {
            double p=0.075;
            double q=0.01;
            double* caliLcc=new double[userNum];
            double* first_second = new double[userNum];
		    double* third_fourth = new double[userNum];
		    double* fifth_sixth = new double[userNum];
            for(int i=0; i<userNum; i++)
            {
                // third_fourth[i] = noisyDegree[i]*(userNum-noisyDegree[i]-1)*p*p*(1-p)*(1-p)*noisyDensity;
                // fifth_sixth[i] = 0.5*(userNum-noisyDegree[i]-1)*(userNum-noisyDegree[i]-2)*pow(1-p, 4)*noisyDensity;
                third_fourth[i] = noisyDegree[i]*(userNum-noisyDegree[i]-1)*p*(1-p)*noisyDensity;
                fifth_sixth[i] = 0.5*(userNum-noisyDegree[i]-1)*(userNum-noisyDegree[i]-2)*pow(1-p, 2)*noisyDensity;
		    }
            for(int i=0; i<userNum; i++)
            {
			    if(third_fourth[i] + fifth_sixth[i] < triCnt[i])
                {
                    first_second[i] = triCnt[i] - third_fourth[i] - fifth_sixth[i];
                }
		    }
            for(int i=0; i<userNum; i++)
            {
			    // double estimatedTriangle = (2.0*first_second[i]-noisyDegree[i]*(noisyDegree[i]-1)*pow(p, 4)*pow(1-p, 2) ) / (2.0*pow(p, 4)*(2*p-1.0));
                double estimatedTriangle = (2.0*first_second[i]-noisyDegree[i]*(noisyDegree[i]-1)*pow(p, 2)*pow(1-p, 1) ) / (2.0*pow(p, 2)*(2*p-1.0));
			    estimatedTriangle = fmax(estimatedTriangle, 0);
			    if(noisyDegree[i]>1)
                {				
                    double onelcc = 2.0*estimatedTriangle / (noisyDegree[i]*(noisyDegree[i]-1));
                    if(onelcc>1)
                    {
                        caliLcc[i]=0.5;
                    }
                    else
                    {
                        caliLcc[i]=onelcc;
                    }
                    						
			    }
		    }
            return caliLcc;
        }
        double* adaptiveCali(double* triCnt, double noisyDensity, double* noisyDegree, double* p, double* q, int userNum)
        {
            for(int i=0;i<userNum;i++)
            {
                double term1 = noisyDegree[i]*((double)userNum-noisyDegree[i]-1);
            }
        }
        double *getTriCnt()
        {
            return triCnt;
        }
        double *getSenCnt()
        {
            return senCnt;
        }
        Metric(string dName)
        {
            setDatasetName(dName);
        }
        protected:
            string dataset;
            double *triCnt;
            double *senCnt;
};

// int main()
// { 
//     Graph rg(1);
//     PerturbGraph perGra(rg,8.0,0.95,3);
//     Metric metric;
//     double **perMatrix=perGra.getNoisyMatrix();
//     // double* p = perGra.getP();
//     // double* q = perGra.getQ();
//     double *triCnt = metric.TriangleCount(perMatrix,rg.getUserNum());
//     // double *triCnt = metric.lapTriangleCnt(perMatrix,rg.getUserNum());
//     // double *estLcc = metric.lccCalibration(triCnt,perGra.getNoisyDensity(),perGra.getNoisyDegreeList(),4,0.8923,rg.getUserNum());
//     // ofstream outfile("./cali_lcc",ios::app);
//     // for(int i=0;i<rg.getUserNum();i++)
//     // {
//     //     outfile<<estLcc[i]<<endl;
//     // }
//     // outfile.close();
//     // outfile.clear();
//     ofstream outfile;
//     outfile.open("./noisyTri_lap.txt",ios::app);
//     for(int i=0;i<rg.getUserNum();i++)
//     {
//         outfile<<triCnt[i]<<endl;
//     }
//     // double **oriMatrix = rg.getMatrix();
//     // double *oriCnt = metric.TriangleCount(oriMatrix,rg.getUserNum());
//     // outfile.close();
//     // outfile.clear();
//     // // double *Lcc = metric.LocalClusteringCoefficient(rg.getDegreeList(),oriCnt,rg.getUserNum());
//     // outfile.open("./oriCount.txt",ios::app);
//     // for(int i=0;i<rg.getUserNum();i++)
//     // {
//     //     outfile<<oriCnt[i]<<endl;
//     // }
//     outfile.close();
//     return 0;
// }