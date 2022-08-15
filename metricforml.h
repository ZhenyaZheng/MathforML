#ifndef METRICFORML_H
#define METRICFORML_H
#include "MyExcept.h"
#include <algorithm>
namespace MyMath
{
    enum AverageStrategy
    {
        Binary,
        Micro,
        Macro
    };
    
    template<AverageStrategy AS, size_t PositiveClass = 1>
    class AUC
    {
    public:
        template<typename DataType>
        static double getAUC(const vector<size_t>& labels, const DataType & predictions, const int& n=2)
        {
            if (AS == Macro)
            {
                double auc = 0.0;
                for (int i = 0; i < n; ++i)
                {
                    vector<size_t>labeli(labels.size());
                    vector<vector<double>>predictioni(2, vector<double>(n * labels.size()));
                    for (int j = 0; j < labels.size(); ++j)
                    {
                        labeli[j] = (labels.at(j) == static_cast<size_t>(i) ? 1 : 0);
                        predictioni[1][j] = predictions[i][j];
                    }
                    auto tempauc = _getAUC(labeli, predictioni);
                    auc += tempauc;
                    predictioni.clear();
                }
                return auc / n;
            }
            else if (AS == Micro)
            {
                vector<size_t>labeli(n * labels.size());
                vector<vector<double>>predictioni(2, vector<double>(n * labels.size()));
                for (int i = 0; i < labels.size(); ++i)//ÐÐ
                {
                    for (int j = 0; j < n; ++j)//ÁÐ
                    {
                        labeli[i * n + j] = labels[i] == j ? 1 : 0;
                        predictioni[1][i * n + j] = predictions[j][i];
                    }
                }
                return _getAUC(labeli, predictioni);
            }
            else if (AS == Binary)
                return _getAUC(labels, predictions);
            else
                MyExcept("AUC getAUC Error : typeerror !");
        }    

    private:
       
        template<typename DataType>
        static double _getAUC(const vector<size_t>& labels, const DataType& predictions)
        {
            std::vector<std::pair<int, double>> indexval;
            int numpost = 0;
            for (int i = 0; i < labels.size(); ++i)
                indexval.emplace_back(i, predictions[1][i]), numpost += labels[i] == PositiveClass ? 1 : 0;
            if (numpost == 0) return 0.0;
            sort(indexval.begin(), indexval.end(), [](const std::pair<int, double>& a, const  std::pair<int, double>& b) {return a.second < b.second; });
            double sumnum = 0.0;
            //int testnum = 0;
            for (int i = 0; i < labels.size(); ++i)
            {
                int numipos = 0;
                if (labels[indexval[i].first] > 0) numipos++;
                int posti = i;
                int sumi = i + 1;
                int j = i + 1;
                for (; j < labels.size(); ++j)
                {
                    if (indexval[j].second == indexval[posti].second)
                    {
                        if (labels[indexval[j].first] > 0)
                            numipos++;
                        sumi += j + 1;
                    }
                    else { i = j - 1; break; }
                }
                if (j == labels.size()) i = j - 1;
                //testnum += numipos;
                sumnum += (double)(sumi * numipos) / (j - posti);
            }
            return  (double)(sumnum - numpost * (numpost + static_cast<double>(1)) / 2) / (numpost * (labels.size() - numpost));
        }

    };

    template<AverageStrategy AS, size_t PositiveClass = 1>
    class LogLoss
    {
    public:
        template<typename DataType>
        static double getLogLoss(const vector<size_t>& labels, const DataType& predictions, const int n=2)
        {
            if (AS == Binary)
            {
                double logloss = 0.0;
                for (int i = 0; i < labels.size(); ++i)
                {
                    double pred_actual = predictions[1][i];
                    pred_actual = std::max(std::min(pred_actual, 1 - 1E-15), 1E-15);
                    logloss += labels[i] == 0 ? (1 - labels[i]) * log(1 - pred_actual) : labels[i] * log(pred_actual);
                }
                return -logloss / labels.size();
            }
            else 
            {
                double logloss = 0.0;
                for (int i = 0; i < labels.size(); ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        double pred_actual = predictions[j][i];
                        pred_actual = std::max(std::min(pred_actual, 1 - 1E-15), 1E-15);
                        logloss += (labels[i] == j ? 1 : 0) * log(pred_actual);
                    }
                }
                return -logloss / labels.size();
            }
        }

    private:
    };
}//MyMath
#endif