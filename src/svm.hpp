//SVM Header File

#include <stdio.h>
#include <algorithm>
#include <random>
#include <vector>

#include "data.hpp"

using namespace std;

class svm
{
public:
  svm(datasets *info);
  

  private:
  void bestHyperParameters();
  void testing();
  void runAlgorithm(double hp1, double hp2, int trainingData);

  vector<double> checkPrecision(vector<int> predictions, vector<int> answers);
  vector<int> getLabels(int dataSet);
  vector<int> getPredictions(int dataset);
  double objective(vector<pair<int,int>> x);


  default_random_engine generator;
  datasets *info;
  vector<double> weights;
  double bestHP1;
  double bestHP2;

  
};