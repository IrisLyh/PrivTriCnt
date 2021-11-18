#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <string>
#include <sstream>
using namespace std;
class lf_gdpr
{
    public:
    double* perturb_degree(double epsilon, double* oriDegree, int userNum)
    {
        cout<<"userNum:"<<userNum<<endl;
        double* ptrDegree = new double[userNum];
        for(int i=0;i<userNum;i++)
        {
            double u = rand()/(double)(RAND_MAX)-0.5;
            int sign;
            if(u>0)
            {
                sign = 1;
            }
            else if(u<0)
            {
                sign = -1;
            }
            else
            {
                sign = 0;
            }
            double noise = (-1)*(2/epsilon)*(double)sign*log(1-2*abs(u));
            ptrDegree[i] = round(oriDegree[i] + noise);
            if(ptrDegree[i]<1)
            {
                ptrDegree[i]=1;
            }
        }
        double max=-1;
        for(int i=0;i<userNum;i++)
        {
            if(ptrDegree[i]>max)
            {
                max = ptrDegree[i];
            }
        }
        maxDegree = max;
        return ptrDegree;
    }
    double **perturb_matrix(string name,double epsilon, double** oriMatix, int userNum, double d_max)
    {
        double **ptrMatrix = new double* [userNum];
        double ptrEdgeCnt=0;
        for(int i=0;i<userNum;i++)
        {
            ptrMatrix[i] = new double[userNum];
        }
        double p = exp(epsilon)/(1+exp(epsilon));
        
        cout<<"p for perturb matrix:"<<p<<endl;
        srand((unsigned)time(NULL));
        for(int i=0;i<userNum;i++)
        {
            for(int j=i+1;j<userNum;j++)
            {
                double u = (double)rand()/(double)RAND_MAX;
                if(u<=p)
                {
                    ptrMatrix[i][j] = oriMatix[i][j];
                    ptrMatrix[j][i] = oriMatix[i][j];
                }
                else
                {
                    ptrMatrix[i][j] = 1-oriMatix[i][j];
                    ptrMatrix[j][i] = 1-oriMatix[i][j];
                }
            }
        }
        for(int i=0;i<userNum;i++)
        {
            for(int j=0;j<userNum;j++)
            {
                if(ptrMatrix[i][j]==1)
                {
                    ptrEdgeCnt = ptrEdgeCnt + 1;
                }
            }
        }
        cout<<"ptrEdge:"<<ptrEdgeCnt<<endl;
        ptrDensity = ptrEdgeCnt / (double)userNum / ((double)userNum -1);
        // cout<<"Creating fil: " << name <<endl; 
        // ofstream ofile(name.data(),ios::app);
        // for(int i=0;i<userNum;i++)
        // {
        //     for(int j=0;j<userNum;j++)
        //     {
        //         ofile<<ptrMatrix[i][j];
        //         if(j<userNum-1)
        //         {
        //             ofile<<" ";
        //         }
        //         if(j == userNum-1)
        //         {
        //             ofile<<endl;
        //         }
        //     }
        // }
        // ofile.close();
        // cout<<"File written!"<<endl;
        return ptrMatrix;
    }
    double* noisyTriCnt(double** ptrMatrix,double** oriMatrix, int userNum, double* ptrDegree)
    {
        double* cnt = new double[userNum];
        double totalCnt = 0;
        for(int i=0;i<userNum;i++)
        {
            cnt[i] = 0;
            for(int j=i+1;j<userNum;j++)
            {
                for(int k=j+1;k<userNum;k++)
                {
                    
                    if(j!=i && j!=k && i!=k)
                    {
                        if(oriMatrix[i][j]==1 && oriMatrix[i][k]==1 && ptrMatrix[j][k]==1)
                        {
                            cnt[i]++;
                        }
                    }
                }
            }
        }
        double noTr=0;
        for(int i=0;i<userNum;i++)
        {
            noTr = noTr+cnt[i];
        }
        cout<<"noisy triangle cnt:"<<noTr<<endl;
        return cnt;
    }
    double cali_cnt(double* cnt, double epsilon,double* ptrDegree, int userNum, double ptrDensity)
    {
        // string path="./ptrCnt.txt";
        // ofstream my_file(path.data(),ios::app);
        // my_file<<"epsilon:"<<epsilon<<":"<<endl;
        double density = ptrDensity;
        double p = exp(epsilon)/(exp(epsilon)+1);

        double total_cnt = 0;
        double* first_second = new double[userNum];
        double* third_fourth = new double[userNum];
        double* fifth_sixth = new double[userNum];
        for(int i=0;i<userNum;i++)
        {
            third_fourth[i] = ptrDegree[i]*(userNum-ptrDegree[i]-1)*p*(1-p)*density;
            fifth_sixth[i] = 0.5*(userNum-ptrDegree[i]-1)*(userNum-ptrDegree[i]-2)*pow(1-p,2)*density;
        }
        for(int i=0;i<userNum;i++)
        {
            if(third_fourth[i]+fifth_sixth[i]<cnt[i])
            {
                first_second[i]=cnt[i]-third_fourth[i]-fifth_sixth[i];
            }
        }
        for(int i=0; i<userNum; i++)
        {
            double tmp= (first_second[i]-0.5*ptrDegree[i]*(ptrDegree[i]-1)*pow(p, 2)*pow(1-p, 1) ) / (pow(p, 2)*(2*p-1.0));
            if(tmp<0)
            {
                tmp = 0;
            }
            // my_file<<tmp<<endl;
            total_cnt=total_cnt+tmp;
        }
        // my_file.close();
        return total_cnt;
    }

    double newCaliCnt(double *cnt, double eps1,double eps2,double ** oriMatrix, int userNum,double dmax)
    {
        srand((int)time(0));
        double *pathCnt= new double[userNum];
        for(int i=0;i<userNum;i++)
        {
            pathCnt[i]=0;
            for(int j=i+1;j<userNum;j++)
            {
                for(int k=j+1;k<userNum;k++)
                {
                    if(i!=j && i!=k && j!=k)
                    {
                        if(oriMatrix[i][j]==1 && oriMatrix[i][k]==1)
                        {
                            pathCnt[i]=pathCnt[i]+1;
                        }
                    }
                }
            }
        }
        double p=1/(exp(eps1)+1);
        cout<<"probability:"<<p<<endl;
        double totalCnt=0;
        for(int i=0;i<userNum;i++)
        {
            double u = rand()/(double)(RAND_MAX)-0.5;
            int sign;
            if(u>0)
            {
                sign = 1;
            }
            else if(u<0)
            {
                sign = -1;
            }
            else
            {
                sign = 0;
            }
            double noise = (-1)*(dmax/eps2)*(double)sign*log(1-2*abs(u));
            totalCnt = totalCnt + (cnt[i]-p*pathCnt[i]+noise)/(1-2*p);
        }
        return totalCnt;
    }

    double getTCnt()
    {
        return tCnt;
    }

    double getMaxPDegree()
    {
        return maxDegree;
    }

    lf_gdpr(int datasetNo, double** oriMatrix, double epsilon, double* oriDegree)
    {
        double epsilon1 = 0;
        double epsilon2 = 0;
        double alpha = 0;
        int userNum = 0;
        // int index = index = ceil(epsilon - 1);
        // double fb[10] = {0.7071,0.8098,0.8614,0.8923,0.9126,0.9268,0.9372,0.9449,0.95,0.95};
        // double ee[10] = {0.6929,0.7927,0.8460,0.8792,0.9016,0.9175,0.9292,0.9379,0.94,0.94};
        // double ap[10] = {0.6929,0.7927,0.8460,0.8792,0.9016,0.9175,0.9291,0.9379,0.94,0.94};
         switch (datasetNo)
        {
        case 1:
            userNum = 4039;
            break;
        case 2:
            userNum = 36692;
            break;
        case 3:
            userNum = 18772;
            break;
        case 4:
            userNum = 7115;
            break;
        case 5:
            userNum = 27770;
            break;
        default:
            break;
        }
        cout<<"[In]***************"<<endl;
        double epsDegree = epsilon*0.2;
        epsilon1= 0.4*epsilon;
        epsilon2=0.4*epsilon;
        double **ptrMatrix;
        double *ptrDegree;
        ptrMatrix = new double*[userNum];
        for(int i=0;i<userNum;i++)
        {
            ptrMatrix[i] = new double[userNum];
        }
        ptrDegree = new double[userNum];
        cout<<"Perturbing degree..."<<endl;
        ptrDegree = perturb_degree(epsDegree,oriDegree,userNum);
        cout<<"Perturbing matrix..."<<endl;
        std::ostringstream os,nm; 
        os << epsilon; 
        nm << datasetNo;
        std::string str = os.str();
        string str1 = nm.str();
        string ptrMat_path ="./" + str1 + "_" + str + "_ptrMat.txt";
        cout<<"perturb matrix:"<<epsilon1<<endl;
        ptrMatrix = perturb_matrix(ptrMat_path,epsilon1,oriMatrix,userNum,maxDegree);
        double* ptrCnt = new double[userNum];
        cout<<"Calculating triangle..."<<endl;
        cout<<"Noisy counting:"<<epsilon2<<endl;
        ptrCnt = noisyTriCnt(ptrMatrix,oriMatrix,userNum,ptrDegree);
        cout<<"Calibrating triangle..."<<endl;
        tCnt = newCaliCnt(ptrCnt, epsilon1,epsilon2,oriMatrix,userNum,maxDegree);
    }
    protected:
        double tCnt;
        double ptrDensity;
        double maxDegree;
};