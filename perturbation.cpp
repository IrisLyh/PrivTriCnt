#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <typeinfo>
#include <vector>
#include<algorithm>
#include "data.cpp"
using namespace std;
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
bool cmp(pair<int, double>a, pair<int, double>b)
{
    return a.second > b.second;
}
class PerturbGraph
{
     public:
        double** perturbMatrixRR(double eps)
            {
                double **noiMatrix = new double* [userNum];
                double p = exp(eps)/(exp(eps)+1);
                for(int i=0;i<userNum;i++)
                {
                    noiMatrix[i]=new double[userNum];
                }
                srand((unsigned)time(NULL));
                for(int i=0;i<userNum;i++)
                {
                    noisyEdgeCount=0;
                    for(int j=i+1;j<userNum;j++)
                    {
                        double u = rand()/(double)RAND_MAX;
                        if(u<=p)
                        {
                            noiMatrix[i][j]=originMatrix[i][j];
                            noiMatrix[j][i]=originMatrix[i][j];
                        }
                        else
                        {
                            noiMatrix[i][j]=1-originMatrix[i][j];
                            noiMatrix[j][i]=1-originMatrix[i][j];
                        }
                    }
                    for(int j=0;j<userNum;j++)
                    {
                        if(noiMatrix[i][j]==1)
                        {
                            noisyEdgeCount++;
                        }
                    }
                }
                noisyDensity=(double)noisyEdgeCount/(double)(userNum*(userNum-1));
                return noiMatrix;
            }

            double** perturbMatrixAdaptive(double* noiDegree, double eps)
            {
                double **noiMatrix = new double* [userNum];
                for(int i=0;i<userNum;i++)
                {
                    noiMatrix[i]=new double[userNum];
                }
                p=new double[userNum];
                q=new double[userNum];
                LaplaceNoise lapMatrix(userNum,userNum,1,eps);
                double **LapNoise = lapMatrix.getNoise();
                double noiEdgeCount = 0;
                double threshold = 1.0/(exp(eps)+1.0);
                for(int i=0;i<userNum;i++)
                {
                    double rho = noiDegree[i]/(double)userNum;
                    double theta = 0;
                    if(rho<threshold)
                    {
                        theta = (1.0/eps)*log((1.0-rho+rho*exp(eps))/(2*rho));
                        // cout<<"i:"<<i<<" threshold:"<<threshold<<" theta:"<<theta<<" rho:"<<rho<<" degree:"<<noiDegree[i]<<endl;
                        if(theta>=1.0)
                        {
                            p[i]=0.5*exp((-1)*eps*(theta-1));
                            q[i]=0.5*exp((-1)*eps*theta);
                        }
                        else
                        {
                            cout<<"[ERROR]Wrong theta for theta > 1"<<endl;
                            exit(0);
                        }
                    }
                    else if(rho>=threshold)
                    {
                        theta = (1/(2*eps))*log((1-rho)/rho)+0.5;
                        if(theta<=1)
                        {
                            p[i]=1-0.5*exp(eps*(theta-1));
                            q[i]=0.5*exp((-1)*eps*theta);
                        }
                        else
                        {
                            cout<<"[ERROR]Wrong theta for theta < 1"<<endl;
                            exit(0);
                        }
                    }
                    for(int j=0;j<userNum;j++)
                    {
                        noiMatrix[i][j] = originMatrix[i][j] + LapNoise[i][j];
                        if(noiMatrix[i][j]>=theta)
                        {
                            noiMatrix[i][j]=1;
                            noiEdgeCount++;
                        }
                        else
                        {
                            noiMatrix[i][j]=0;
                        }
                    }
                }
                noisyDensity = noiEdgeCount/(userNum*(userNum-1));             
                return noiMatrix;
            }
            
            double** perturbLapOnly(double eps)
            {
                double **noiMatrix = new double* [userNum];
                for(int i=0;i<userNum;i++)
                {
                    noiMatrix[i]=new double[userNum];
                }
                p=new double[userNum];
                q=new double[userNum];
                LaplaceNoise lapMatrix(userNum,userNum,1,eps);
                double **LapNoise = lapMatrix.getNoise();
                srand(time(NULL));
                for(int i=0;i<userNum;i++)
                {
                    for(int j=i+1;j<userNum;j++)
                    {
                        noiMatrix[i][j]=originMatrix[i][j]+LapNoise[i][j];
                        noiMatrix[i][i]=0;
                        noiMatrix[j][i]=noiMatrix[i][j];
                    }
                }
                for(int i=0;i<userNum;i++)
                {
                    double max = -10;
                    double min = 10;
                    for(int j=0;j<userNum;j++)
                    {
                        if(noiMatrix[i][j]>max)
                        {
                            max=noiMatrix[i][j];
                        }
                        if(noiMatrix[i][j]<min)
                        {
                            min=noiMatrix[i][j];
                        }
                    }
                    for(int k=0;k<userNum;k++)
                    {
                        noiMatrix[i][k]=(noiMatrix[i][k]-min)/(max-min);
                        double ins = rand()/(double)RAND_MAX;
                        if(ins<noiMatrix[i][k])
                        {
                            noiMatrix[i][k]=1;
                        }
                        else
                        {
                            noiMatrix[i][k]=0;
                        }
                    }
                }
                return noiMatrix;
            }

            double *perturDegree(double eps)
            {
                LaplaceNoise LapDis = LaplaceNoise(userNum,1,2,eps);
                double * noisyList = new double[userNum];
                double ** LapNosie = LapDis.getNoise();
                for(int i=0;i<userNum;i++)
                {
                    noisyList[i]=originDegreeList[i]+LapNosie[0][i];
                    noisyList[i] = round(noisyList[i]);
                    if(noisyList[i]<=0)
                    {
                        noisyList[i]=1;
                    }
                }
                return noisyList;
            }
            void setNoisyMatrixDegreeList()
            {
                noisyMatrixDegreeList = new double[userNum];
                for(int i=0;i<userNum;i++)
                {
                    noisyMatrixDegreeList[i]=0;
                    for(int j=0;j<userNum;j++)
                    {
                        if(noisyMatrix[i][j]==1)
                        {
                            noisyMatrixDegreeList[i]++;
                        }
                    }
                }
            }
            double ** getNoisyMatrix()
            {
                return noisyMatrix;
            }
            double getEpsilon()
            {
                return epsilon;
            }
            double getAlpha()
            {
                return alpha;
            }
            double getNoisyDensity()
            {
                return noisyDensity;
            }
            int getNoisyEdgeCount()
            {
                return noisyEdgeCount;
            }
            double* getNoisyDegreeList()
            {
                return noisyDegreeList;
            }
            double* getNoisyMatrixDegreeList()
            {
                return noisyMatrixDegreeList;
            }
            double* getP()
            {
                if(perturbMethod!=2)
                {
                    cout<<"[ERROR]current method has no q"<<endl;
                    exit(0);
                }
                else
                {
                    return p;
                }
            }
            double* getQ()
            {
                if(perturbMethod!=2)
                {
                    cout<<"[ERROR]current method has no q"<<endl;
                    exit(0);
                }
                else
                {
                    return q;
                }
            }
            void perturbRR(double eps1,double eps2)
            {
                noisyDegreeList = perturDegree(eps1);
                noisyMatrix = perturbMatrixRR(eps2);
                setNoisyMatrixDegreeList();
            }
            void perturbAdaptive(double eps1, double eps2)
            {
                noisyDegreeList = perturDegree(eps1);
                noisyMatrix = perturbMatrixAdaptive(noisyDegreeList,eps2);
            }
            void perturbLap(double eps1, double eps2)
            {
                noisyDegreeList = perturDegree(eps1);
                noisyMatrix = perturbLapOnly(eps2);
            }
            PerturbGraph(Graph oriGraph,double eps,double alp, int method)
            {
                perturbMethod = method;
                userNum = oriGraph.getUserNum();
                originMatrix = oriGraph.getMatrix();
                originDegreeList = oriGraph.getDegreeList();
                alpha = alp;
                double eps1 = eps*(1-alpha);
                double eps2 = eps*alpha;
                switch (method)
                {
                case 1:
                    perturbRR(eps1,eps2);
                    break;
                case 2:
                    perturbAdaptive(eps1,eps2);
                case 3:
                    perturbLap(eps1,eps2);
                default:
                    break;
                }
                setNoisyMatrixDegreeList();
            }
        protected:
            double** noisyMatrix;
            double** originMatrix;
            double epsilon;
            double* noisyDegreeList;
            double* noisyMatrixDegreeList;
            double* originDegreeList;
            double noisyDensity;
            double alpha;//for epsilon alignment
            double* theta;  //for adaptive pertubation
            double* p;  //for adaptive perturbation
            double* q;  //for adaptive perturbation
            int noisyEdgeCount;
            int userNum;
            int perturbMethod;
};