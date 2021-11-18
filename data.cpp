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
#endif
#include "utils.h"
using namespace std;
class Graph
{
    public:
        void setDataSet(int data)
        {
            switch (data)
            {
            case 1:
                dataset = "facebook";
                userNum = 4039;
                break;
            case 2:
                dataset = "email-Enron";
                userNum = 36692;
                break;
            case 3:
                dataset = "ca-AstroPh";
                userNum = 18772;
                break;
            case 4:
                dataset = "wiki-Vote";
                userNum = 7115;
                break;
            case 5:
                dataset = "cit-HepTh";
                userNum = 27770;
                break;
            default:
                break;
            }
        }
        double **getMatrix()
        {
            return matrix;
        }
        int getUserNum()
        {
            return userNum;
        }
        int getEdgeCount()
        {
            return edgeCount;
        }
        double getDensity()
        {
            return density;
        }
        string getDatasetName()
        {
            return dataset;
        }
        double * getDegreeList()
        {
            double* degreeList=new double[userNum];
            string path = "./"+dataset+"_oriDegree.txt";
            if(file_check(path)==1)
            {
                ifstream infile_tri(path.data());
                for(int i=0;i<userNum;i++)
                {
                    infile_tri>>degreeList[i];
                }
                infile_tri.close();
            }
            else
            {
                for(int i=0;i<userNum;i++)
                {
                    degreeList[i]=0;
                    for(int j=0;j<userNum;j++)
                    {
                        if(matrix[i][j]==1.0)
                        {
                            degreeList[i]++;
                        }
                    }
                }
                ofstream outfile_deg(path.data());
                for(int i=0;i<userNum;i++)
                {
                    outfile_deg<<degreeList[i]<<endl;
                }
                outfile_deg.close();
            }
            return degreeList;
        }
        double getDegreeOne(int user)
        {
            if(user>=0 && user<userNum)
            {
                return degreeList[user];
            }
            else
            {
                cout<<"[ERROR]Wrong userID!!"<<endl;
                exit(0);
            }
            
        }
        double ** readFile(string dataset) //目前主要针对的是facebook
        {
            cout<<"[INFO] Reading graph data file! Dataset name:"<<dataset<<endl;
            string path = "./dataset/"+dataset+".txt";
            ifstream infile(path.data());
            double **mat = new double*[userNum];
            char pair[100];
            int check_number = 0;
            for(int i=0;i<userNum;i++)
            {
                mat[i]=new double[userNum];
            }
            edgeCount=0;
            while(!infile.eof())
            {
                check_number++;
                infile.getline(pair,99);
                int indicatior=0,nodeA=0,nodeB=0,tmp=0;
                while(indicatior<100)
                {
                    if(pair[indicatior]==13)
                    {
                        pair[indicatior]='\0';
                    }
                    if(pair[indicatior]=='\0')
                    {
                        if(indicatior==0)
                        {
                            break;
                        }
                        else
                        {                         
                            nodeB=tmp;
                            break;
                        }
                    }
                    else if(pair[indicatior]==' ')
                    {
                        nodeA=tmp;
                        tmp=0;
                    }
                    else
                    {
                        tmp=tmp*10+pair[indicatior]-48;
                    }
                    indicatior++;
                }
                if(nodeA!=nodeB)
                {
                    mat[nodeA][nodeB]=1;
                    mat[nodeB][nodeA]=1;
                    edgeCount++;
                }
                else
                {
                    mat[nodeA][nodeB]=0;
                }
            }
            density=(double)edgeCount/(double)(userNum*(userNum-1));
            cout<<"edgeCount:"<<edgeCount<<endl;
            infile.close();
            return mat;
        }
        int** setELV(double** matrix, int userNum)
        {
            int** ELV = new int* [userNum];
            for(int i=0;i<userNum;i++)
            {
                ELV[i] = new int[userNum];
                for(int j=0;j<userNum;j++)
                {
                    ELV[i][j]=-1;
                }
            }
            for(int i=0;i<userNum;i++)
            {
                for(int j=0;j<userNum;j++)
                {
                    if(matrix[i][j]==1)
                    {
                        ELV[i][j]=1;
                        for(int k=0;k<userNum;k++)
                        {
                            if(k!=i && matrix[i][k]==0 && matrix[j][k]==1)
                            {
                                ELV[i][k]=2;
                            }
                        }
                    }
                }
            }
            return ELV;
        }
        int** getELV()
        {
            return ELV;
        }
        Graph(int data)
        {
            setDataSet(data);
            matrix = readFile(dataset);
            degreeList = getDegreeList();
        }

    protected:
        string dataset;
        int userNum;
        int edgeCount;
        double density;
        double **matrix;
        double *degreeList;
        int ** ELV;
};