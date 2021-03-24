#ifndef AGGLOMERATIVECLUSTERING1_H
#define AGGLOMERATIVECLUSTERING1_H

#include<string.h>
#include<Eigen/Dense>
#include<cmath>
#include<vector>
#include<iostream>

using namespace std;

class agglomerativeClustering1
{
    public:
        agglomerativeClustering1();
        string method;
        string affinity;
        unsigned int numOfClusters;
        Eigen::MatrixXd fit(Eigen::MatrixXd data_,unsigned int numOfClusters_,string method_,string affinity_);
        Eigen::MatrixXd trackerMatrix;

    private:
        Eigen::MatrixXd data;
        int rows;
        int columns;
        int distDim;
        Eigen::VectorXd data1;
        Eigen::VectorXd data2;
        Eigen::VectorXd minVec;
        Eigen::MatrixXd distMatrix;
        Eigen::MatrixXd createArray_dist(Eigen::MatrixXd data);
        double findDistance(Eigen::VectorXd data1,Eigen::VectorXd data2,string affinity);
        void updateDistance(Eigen::MatrixXd &distMatrix,string method);
        Eigen::VectorXd findmin(Eigen::MatrixXd distMatrix);
        void removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove);
        void removeColumn(Eigen::MatrixXd &matrix, unsigned int colToRemove);
        void centroid(int loc1,int loc2,Eigen::MatrixXd &centTracker);
        void tracker(int loc1,int loc2);

};

#endif // AGGLOMERATIVECLUSTERING1_H
