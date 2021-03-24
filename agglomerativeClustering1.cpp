#include "agglomerativeClustering1.h"

agglomerativeClustering1::agglomerativeClustering1()
{
    //ctor
    cout<<"Object Created"<<endl;
}

void agglomerativeClustering1::centroid(int loc1,int loc2,Eigen::MatrixXd &centTracker){

int centroidRows = 0;

for(int i=0;i<rows;i++){

    if(centTracker(i,1)== loc1 || centTracker(i,1)==loc2){

        centroidRows++;
    }
}

Eigen::MatrixXd cent(centroidRows,columns);
int index=0;

for(int i=0;i<rows;i++){

    if(centTracker(i,1)== loc1 || centTracker(i,1)==loc2){

        cent.row(index) = data.row(i);
        index++;

    }
}

index=index-1;

for(int i=0;i<centroidRows-1;i++){

    cent.row(i+1)=cent.row(i)+cent.row(i+1);
}

cent.row(index)=cent.row(index)/(index+1);


for(int i=0;i<rows;i++){

    if(centTracker(i,1)==loc2){

        centTracker(i,1)=loc1;
    }
    if(centTracker(i,1)>loc2){

        centTracker(i,1)=centTracker(i,1)-1;
    }
}
for(int i=0;i<rows;i++){

    if(centTracker(i,1)== loc1 ){

        centTracker.block(i,2,1,columns)=cent.row(centroidRows-1);
}

}
Eigen::VectorXd data1_ = cent.row(index);

for(int i=loc2;i<distDim;i++){
    if(i!=loc1 || i!=loc2){
        Eigen::VectorXd data2_ = centTracker.block(i,2,1,columns).transpose();

        distMatrix(loc1,i)=findDistance(data1_,data2_,affinity);

    }

    }

}

void agglomerativeClustering1::removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();
    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
    matrix.conservativeResize(numRows,numCols);
}

void agglomerativeClustering1::removeColumn(Eigen::MatrixXd &matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;
    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
    matrix.conservativeResize(numRows,numCols);
}

Eigen::VectorXd agglomerativeClustering1::findmin(Eigen::MatrixXd distMatrix){

  double mini = distMatrix.minCoeff() ;

  int miniIndex,minjindex;

  int flag=0;

  for(int i=0;i<distDim-1;i++){

    if(flag == 1){
        break;
    }

    for(int j=i+1;j<distDim;j++){

        if (distMatrix(i,j) == mini){

            miniIndex = i;
            minjindex = j;
            flag = 1;
            break;
        }
    }
  }

  Eigen::VectorXd a(3);
  a << mini,miniIndex,minjindex;
  return a;
}

void agglomerativeClustering1::tracker(int loc1,int loc2){

for (int i=0;i<rows;i++){

    if(trackerMatrix(i,1)== loc2){

        trackerMatrix(i,1)=loc1;
    }
    if(trackerMatrix(i,1)>loc2){

        trackerMatrix(i,1)=trackerMatrix(i,1)-1;
    }
}
}

void agglomerativeClustering1::updateDistance(Eigen::MatrixXd &distMatrix,string method){

// Single - Link
if (method == "singleLink"){
    minVec = findmin(distMatrix);

    for(int i=minVec(2);i<distDim;i++){

    if(distMatrix(minVec(1),i) > distMatrix(minVec(2),i)){

        distMatrix(minVec(1),i) = distMatrix(minVec(2),i);
    }
}
}
else if (method == "completeLink"){
        cout<<"chk6"<<endl;
minVec = findmin(distMatrix);
cout<<minVec<<endl;
    for(int i=minVec(2);i<distDim;i++){

    if(distMatrix(minVec(1),i) < distMatrix(minVec(2),i)){

        distMatrix(minVec(1),i) = distMatrix(minVec(2),i);
    }

}
   cout<<"chk8"<<endl;
}
else if (method == "averageLink"){
    minVec = findmin(distMatrix);
    for(int i=minVec(2);i<distDim;i++){

        distMatrix(minVec(1),i) = (distMatrix(minVec(1),i)+distMatrix(minVec(2),i))/2.0;
    }
}

else if (method == "centroid"){

    int matCols = 2+columns;
    Eigen::MatrixXd centTracker(rows,matCols);
    centTracker.col(0) = trackerMatrix.col(0);
    centTracker.col(1) = trackerMatrix.col(1);
    centTracker.block(0,2,rows,columns)=data;
cout<<distMatrix<<endl;
for(int i=0;i< (rows-numOfClusters) ;i++){

    minVec = findmin(distMatrix);
    cout<<distMatrix.rows()<<distMatrix.cols()<<endl;
    centroid(minVec(1),minVec(2),centTracker);
    removeRow(distMatrix,minVec(2));
    removeColumn(distMatrix,minVec(2));
    distDim--;
}
trackerMatrix.col(0) = centTracker.col(0);
trackerMatrix.col(1) = centTracker.col(1);

}
else{

    cout<<"Couldn't find the [method] parameter to proceed"<<endl;
    exit(0);
}

if(method !="centroid"){
cout<<"chk12"<<endl;
cout<<minVec<<endl;
removeColumn(distMatrix,minVec(2));
cout<<"chk10"<<endl;
removeRow(distMatrix,minVec(2));
cout<<"chk9"<<endl;

distDim--;



tracker(minVec(1),minVec(2));
}
}

double agglomerativeClustering1::findDistance(Eigen::VectorXd data1,Eigen::VectorXd data2,string affinity){

    if(affinity == "manhattan"){
    return (data1-data2).array().abs().sum();
    }

    if(affinity == "euclidean"){
    return sqrt((data1-data2).array().square().sum());
    }
    return 2;
}

Eigen::MatrixXd agglomerativeClustering1::createArray_dist(Eigen::MatrixXd data){

   for(int i = 0; i < rows-1;i++){
    data1 = data.row(i);
    for (int j = i+1 ; j < rows ; j++){
        data2 = data.row(j);
        distMatrix(i,j) = findDistance(data1,data2,affinity);
    }
   }
   return distMatrix;
}

Eigen::MatrixXd agglomerativeClustering1::fit(Eigen::MatrixXd data_,unsigned int numOfClusters_,string method_,string affinity_){

cout<<"chk1"<<endl;
method = method_;
affinity = affinity_;
numOfClusters = numOfClusters_;
data = data_;
rows=data.rows();
columns=data.cols();
distDim=data.rows();
distMatrix =Eigen::MatrixXd::Constant(rows,rows,10000);
cout<<"chk2"<<endl;
distMatrix = createArray_dist(data);
cout<<"chk3"<<endl;
trackerMatrix = Eigen::MatrixXd(rows,2);
cout<<"chk4"<<endl;
for(int i=0;i<rows;i++){

    trackerMatrix(i,0) = i;
    trackerMatrix(i,1) = i;

}
cout<<trackerMatrix<<endl;
cout<<"chk5"<<endl;

if(method!="centroid"){

for (unsigned int i = 0 ; i<(rows-numOfClusters) ; i++){

    cout<<"iter"<<i<<endl;
    updateDistance(distMatrix,method);



}
}
else if(method=="centroid"){
    cout<<"in fit"<<endl;
    updateDistance(distMatrix,method);
}
return trackerMatrix;
}

