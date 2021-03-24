#include <iostream>
#include "agglomerativeClustering1.h"
#include<Eigen/Dense>
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>

using namespace std;
using namespace Eigen;

void saveData(string fileName, Eigen::MatrixXd  matrix)
{
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

Eigen::MatrixXd openData(string fileToOpen)
{

    // later on, this vector is mapped into the Eigen matrix format
    vector<double> matrixEntries;

    // in this object we store the data from the matrix
    ifstream matrixDataFile(fileToOpen);

    // this variable is used to store the row of the matrix that contains commas
    string matrixRowString;

    // this variable is used to store the matrix entry;
    string matrixEntry;

    // this variable is used to track the number of rows
    int matrixRowNumber = 0;


    while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
    {
        stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

        while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
        {
            matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
        }
        matrixRowNumber++; //update the column numbers
    }

    return Eigen::Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}
int main()
{
    Eigen::MatrixXd data;
    Eigen::MatrixXd result;
    data = openData("sample.csv");
    agglomerativeClustering1 ahc;
    result = ahc.fit(data,3,"centroid","euclidean");
    cout<<result<<endl;
    saveData("opclsdtscen.csv",result);
    cout << "Data Saved" << endl;
    return 0;
}
