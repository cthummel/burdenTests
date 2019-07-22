//SVM 

#include <iostream>
#include "svm.hpp"

svm::svm(datasets *inputData)
{
    generator.seed(0);
    info = inputData;
    for(int i = 0; i < 250; i++)
    {
        weights.push_back(0);
    }
    bestHyperParameters();
    testing();
}




void svm::runAlgorithm(double hp1, double hp2, int dataset)
{
    //Training Data
    if(dataset == -1)
    {
        shuffle(info->trainData.begin(), info->trainData.end(), generator);
        for(int i = 0; i < info->trainData.size(); i++)
        {
            double sign = 0;
            int label = info->trainData[i][0].second;
            for(int j = 1; j < info->trainData[i].size(); j++)
            {
                double x = info->trainData[i][j].second;
                sign += x * weights[info->trainData[i][j].first];
            }
            sign = sign * label;

            //No matter the sign we update the weights a base amount.
            for(int j = 0; j < weights.size(); j++)
            {
                weights[j] = (1 - hp1) * weights[j];
            }

            if(sign <= 1)
            {
                for (int j = 1; j < info->trainData[i].size(); j++)
                {
                    weights[info->trainData[i][j].first] += hp1 * hp2 * label * info->trainData[i][j].second;
                }
            }
        }
    }
    //Doing cross validation
    else
    {
        shuffle(info->combinedData[dataset].begin(), info->combinedData[dataset].end(), generator);
        for(int i = 0; i < info->combinedData[dataset].size(); i++)
        {
            double sign = 0;
            int label = info->combinedData[dataset][i][0].second;
            for(int j = 1; j < info->combinedData[dataset][i].size(); j++)
            {
                double x = info->combinedData[dataset][i][j].second;
                sign += x * weights[info->combinedData[dataset][i][j].first];
            }
            sign = sign * label;

            //No matter the sign we update the weights a base amount.
            for(int j = 0; j < weights.size(); j++)
            {
                weights[j] = (1 - hp1) * weights[j];
            }

            if(sign <= 1)
            {
                for (int j = 1; j < info->combinedData[dataset][i].size(); j++)
                {
                    weights[info->combinedData[dataset][i][j].first] += hp1 * hp2 * label * info->combinedData[dataset][i][j].second;
                }
            }
        }
    }    
}

void svm::bestHyperParameters()
{
    vector<double> hp1 = {1, .1, .01, .001};
    vector<double> hp2 = {1000, 100, 10, 1, .1, .01};
    vector<double> precisionAccuracy;
    vector<double> recallAccuracy;
    vector<double> F1Accuracy;
    vector<pair<double,double> > parameters;

    for(int i = 0; i < hp1.size(); i++)
    {
        for(int j = 0; j < hp2.size(); j++)
        {
            double avgPrecision = 0;
            double avgRecall = 0;
            double avgF1 = 0;
            
            for(int k = 0; k < info->combinedData.size(); k++)
            {
                double change = 20;
                for(int m = 0; m < weights.size(); m++)
                {
                    weights[m] = 0;
                }
                //Epochs
                for(int m = 0; m < change; m++)
                {
                    //cout << "lr = " << (1 - hp1[i]/(1 + m)) << endl;
                    runAlgorithm(hp1[i]/(1 + m), hp2[j], k);
                }
                vector<int> predictions = getPredictions(k);
                vector<double> results = checkPrecision(predictions, getLabels(k));
                avgPrecision += results[0];
                avgRecall += results[1];
                avgF1 += results[2];
            }
            avgPrecision = avgPrecision / info->combinedData.size();
            avgRecall = avgRecall / info->combinedData.size();
            avgF1 = avgF1 / info->combinedData.size();
            precisionAccuracy.push_back(avgPrecision);
            recallAccuracy.push_back(avgRecall);
            F1Accuracy.push_back(avgF1);
            parameters.push_back(pair<double,double>(hp1[i], hp2[j]));
            cout << "Average Values for (" << hp1[i] << ", " << hp2[j] << ") = [" 
            << avgPrecision << ", " << avgRecall << ", " << avgF1 << "]" << endl;
        }
    }

    double max = 0;
    int index = 0;
    for(int i = 0; i < F1Accuracy.size(); i++)
    {
        if(F1Accuracy[i] > max)
        {
            max = F1Accuracy[i];
            bestHP1 = parameters[i].first;
            bestHP2 = parameters[i].second;
            index = i;
        }
    }
    cout << "Best hyper parameter values are (" << bestHP1 << ", " << bestHP2 << ")" << endl; 
    cout << "Values were (" << precisionAccuracy[index] << ", " << recallAccuracy[index] << ", " << F1Accuracy[index] << ")" << endl;
}


void svm::testing()
{
    for(int i = 0; i < weights.size(); i++)
    {
        weights[i] = 0;
    }
    //20 epochs
    for(int i = 0; i < 20; i++)
    {
        runAlgorithm(bestHP1, bestHP2, -1);
    }
    vector<int> predictions = getPredictions(-1);
    vector<double> results = checkPrecision(predictions, getLabels(-1));
    double precision = results[0];
    double recall = results[1];
    double fvalue = results[2];
    cout << "Values for run with training data on test data: (" << precision << ", " << recall << ", " << fvalue << ")"<< endl;
}


vector<int> svm::getPredictions(int dataset)
{
    vector<int> results;

    if(dataset == -1)
    {
        for(int i = 0; i < info->testData.size(); i++)
        {
            double predictedLabel = 0;
            for(int j = 1; j < info->testData[i].size(); j++)
            {
                predictedLabel += weights[info->testData[i][j].first] * info->testData[i][j].second;
            }
            
            if(predictedLabel >= 0)
            {
                results.push_back(1);
            }
            else
            {
                results.push_back(-1);
            }
        }
    }
    else
    {
        for(int i = 0; i < info->hyperData[dataset].size(); i++)
        {
            double predictedLabel = 0;
            for(int j = 1; j < info->hyperData[dataset][i].size(); j++)
            {
                predictedLabel += weights[info->hyperData[dataset][i][j].first] * info->hyperData[dataset][i][j].second;
            }
            if(predictedLabel >= 0)
            {
                results.push_back(1);
            }
            else
            {
                results.push_back(-1);
            }
        }
    }

    return results;
}

//Calculates precision, recall, and f-value for a given set of predictions and answers.
vector<double> svm::checkPrecision(vector<int> predictions, vector<int> answers)
{
    vector<double> results;
    double TP = 0;
    double FP = 0;
    double FN = 0;
    double TN = 0;
    double precision;
    double recall;
    double fvalue;
    
    for(int i = 0; i < predictions.size(); i++)
    {
        if(predictions[i] == answers[i])
        {
            if(predictions[i] == 1)
            {
                TP++;
            }
            else
            {
                TN++;
            }
        }
        else
        {
            if(predictions[i] == 1)
            {
                FP++;
            }
            else
            {
                FN++;
            }
        }
    }
    if(TP == 0 && FP == 0)
    {
        //cout << predictions.size() << " = " << FN << " + " <<  TN << " = "<< (TN + FN) << " ?" << endl;
        for(int i = 0; i < predictions.size(); i++)
        {
            if(predictions[i] == 1)
            {
                cout << "Yippee" << endl;
            }
            //cout << predictions[i] << ", " << answers[i] << endl;
        }
        //cout << "We making nans lads" << endl;
    }
    precision = TP / (TP + FP);
    recall = TP / (TP + FN);
    fvalue = (2 * (precision * recall)) / (precision + recall);


    results.push_back(precision);
    results.push_back(recall);
    results.push_back(fvalue);

    return results;
}

//Returns the labels for test data.
//dataSet < 0 means test data otherwise we send back hyperdata labels.
vector<int> svm::getLabels(int dataSet)
{
    vector<int> results;
    if(dataSet < 0)
    {
        for (int i = 0; i < info->testData.size(); i++)
        {
            results.push_back(info->testData[i][0].second);
        }
    }
    else
    {
        for (int i = 0; i < info->hyperData[dataSet].size(); i++)
        {
            results.push_back(info->hyperData[dataSet][i][0].second);
        }
    }
    return results;
}

