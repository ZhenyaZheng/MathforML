#include "mymath.h"
#include "metricforml.h"
#include <iostream>
int main()
{
	vector<size_t> labels{0,1,0,0,1,0,1,1,0};
	vector<vector<double>> predictions{ {0.9,0.2,0.8,0.8,0.2,0.6,0.1,0.1,0.9}, {0.1,0.8,0.2,0.2,0.8,0.4,0.9,0.9,0.1} };
	cout << "test1 auc = " << MyMath::AUC<MyMath::Binary>::getAUC(labels, predictions) << endl;
	cout << "test1 logloss = " << MyMath::LogLoss<MyMath::Binary>::getLogLoss(labels, predictions) << endl;
	vector<size_t> labels2{0,1,2,1,2,2,0,1,0};
	vector<vector<double>> predictions2{ 
		{0.8,0.2,0.6,0.6,0.2,0.4,0.1,0.6,0.3}, 
		{0.1,0.6,0.2,0.2,0.2,0.2,0.6,0.3,0.1}, 
		{0.1,0.2,0.2,0.2,0.6,0.4,0.3,0.1,0.6} };
	cout << "test2 mirco auc = " << MyMath::AUC<MyMath::Micro>::getAUC(labels2, predictions2, 3) << endl;
	cout << "test2 marco auc = " << MyMath::AUC<MyMath::Macro>::getAUC(labels2, predictions2, 3) << endl;
	cout << "test2 marco logloss = " << MyMath::LogLoss<MyMath::Binary>::getLogLoss(labels2, predictions2, 3) << endl;
	return 0;
}