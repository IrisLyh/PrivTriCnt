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
    int datasetNo=5;
    Graph dataset(datasetNo);
    Metric metric(dataset.getDatasetName());
    int userNum = dataset.getUserNum();
    double *senCnt = new double[userNum];
    double **mat = new double*[userNum];
    for(int i=0;i<userNum;i++)
    {
        mat[i] = new double[userNum];
    }
    mat = dataset.getMatrix();
    double *degreeList = new double[userNum];
    double *triCnt = new double[userNum];
    degreeList = dataset.getDegreeList();
    // string path = "./oriDegree.txt";
    // bool flag = true;
    // if(flag)
    // {
    //     ofstream outfile_tri(path.data());
    //     for(int i=0;i<userNum;i++)
    //     {
    //         outfile_tri<<degreeList[i]<<endl;
    //     }
    //     outfile_tri.close();
    // }
    // cout<<"counting triangle..."<<endl;
    // triCnt = metric.TriangleCount(mat,userNum);
    cout<<"counting sensitivity..."<<endl;
    senCnt = metric.sensitivityCount(mat,userNum);
    return 0;
}