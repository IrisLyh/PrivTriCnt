
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

using namespace std;

int main()
{
    for(int datasetNo=4;datasetNo<=4;datasetNo++)
    {
        Graph dataset(datasetNo);
        Metric metric(dataset.getDatasetName());
        int userNum = dataset.getUserNum();
        double **mat = new double*[userNum];
        for(int i=0;i<userNum;i++)
        {
            mat[i] = new double[userNum];
        }
        double sample_rate[12]= {0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2};
        mat = dataset.getMatrix();
        for(int rt=0;rt<12;rt++)
        {
            double *sample_sensitivity = new double[userNum];
            double temp_rate = sample_rate[rt];
            sample_sensitivity = metric.sensitivityAfterSample(mat,userNum,temp_rate);
        } 
    }
    return 0;
}