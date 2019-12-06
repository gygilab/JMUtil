#ifndef _JMUTIL_HPP
#define	_JMUTIL_HPP

#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <numeric>
#include <limits>
#include <algorithm>
#include <functional>

#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_statistics.h>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/students_t.hpp>

using namespace std;

namespace JMUtil {
    
    inline double roundTo(double x, double places) {
        double factor = pow(10, places);
        return round(x * factor) / factor;
    }

    template<class T> inline double Lerp(T v0, T v1, T t) {
        return (1 - t)*v0 + t*v1;
    }

    template<class T> inline T quantile(const vector<T>& inData, const double& prob) {
        if (inData.empty()) {
            cerr << "ERROR: cannot determine quantile of empty distribution" << endl;
            exit(EXIT_FAILURE);
        }

        if (1 == inData.size()) {
            return inData[0];
        }

        vector<T> data = inData;
        if (! is_sorted(data.begin(), data.end())) {
            sort(data.begin(), data.end());
        }
        
        T poi = Lerp<T>(-0.5, data.size() - 0.5, prob);

        size_t left = max(int64_t(floor(poi)), int64_t(0));
        size_t right = min(int64_t(ceil(poi)), int64_t(data.size() - 1));

        T datLeft = data.at(left);
        T datRight = data.at(right);

        T quantile = Lerp<T>(datLeft, datRight, poi - left);

        return quantile;
    }
    
    template<class T> inline string toStr(T x, int prec = 5) {
        ostringstream temp;
        temp << fixed << setprecision(prec);
        temp << x;
        return temp.str();
    }
    
    template<class T> inline string toStrE(T x, int prec = 1) {
        ostringstream temp;
        temp << scientific << setprecision(prec);
        temp << x;
        return temp.str();
    }

    template<class T> inline string toStr(vector<T> vec, string delim, int prec = 5) {
        ostringstream oss;
        oss << fixed << setprecision(prec);
        if (!vec.empty()) {
            // Convert all but the last element to avoid a trailing ","
            copy(vec.begin(), vec.end() - 1, ostream_iterator<T>(oss, delim.c_str()));
            // Now add the last element with no delimiter
            oss << vec.back();
        }
        return oss.str();
    }

//    template<class T> inline vector<double> vectorLog10 (vector<T> v) {
//        vector<double> log10v;
//        log10v.resize(v.size());
////        transform (v.begin(), v.end(), log10v.begin(), log10<double>());
//        return log10v;
//    };

    template<class T> inline double vectorSum (vector<T> v) {
        double sum = (double) accumulate(v.begin(), v.end(), 0.0);
        return sum;
    };
    
    template<class T> inline double vectorMean (vector<T> v) {
        double sum = (double) accumulate(v.begin(), v.end(), 0.0);
        return sum/(double) v.size();
    };
    
    template<class T> inline double vectorMean (vector<T> v, vector<T> w) {
        if (v.size() != w.size()) {
            cerr << "ERROR: vectorMean with v of size " << v.size() << " and w of size " << w.size() << endl;
            exit(0);
        }
        vector<double> vw(v.size());
        transform(v.begin(), v.end(), w.begin(), vw.begin(), multiplies<double>());
        
        double vw_sum = (double) accumulate(vw.begin(), vw.end(), 0.0);
        double w_sum = (double) accumulate(w.begin(), w.end(), 0.0);
        return vw_sum / w_sum;
    };

    template<class T> inline double vectorStd (vector<T> v) {
        double vector_mean = vectorMean(v);
        vector<double> mean_diff_squared(v.size());
        for (int i = 0; i < v.size(); i++) {
            mean_diff_squared[i] = pow(((double) v[i] - vector_mean), 2);
        }
        double SSE = accumulate(mean_diff_squared.begin(), mean_diff_squared.end(), 0.0);
        if (mean_diff_squared.size() > 1) {
            return sqrt(SSE/( mean_diff_squared.size() - 1));
        } else {
            return sqrt(SSE/( mean_diff_squared.size() ));
        }
    };

    /* Sqrt of
     * unbiased weighted sample variance, based on equation on wikipedia and described by GSL as
     * \Hat\sigma^2 = ((\sum w_i)/((\sum w_i)^2 - \sum (w_i^2))) \sum w_i (x_i - \Hat\mu)^2
     * Note that this expression reduces to an unweighted variance with the familiar 1/(N-1) factor
     * when there are N equal non-zero weights
     */
    template<class T> inline double vectorStd (vector<T> v, vector<T> w) {
        if (v.size() != w.size()) {
            cerr << "ERROR: vectorStd with v of size " << v.size() << " and w of size " << w.size() << endl;
            exit(0);
        }
        if (v.size() <= 1)
            return 0.0;
        
        double vector_mean = vectorMean(v, w);

        // normalize weights
        double sum_w = accumulate(w.begin(), w.end(), 0.0);
        vector<double> w_normd(v.size());
        transform(w.begin(), w.end(), w_normd.begin(), bind1st(multiplies<double>(), 1.0/sum_w));

        double sum_w_squared = pow(accumulate(w_normd.begin(), w_normd.end(), 0.0), 2);
        vector<double> squared_w(v.size());
        transform(w_normd.begin(), w_normd.end(), w_normd.begin(), squared_w.begin(), multiplies<double>());
        double sum_squared_w = (double) accumulate(squared_w.begin(), squared_w.end(), 0.0);

        vector<double> mean_diff_squared(v.size());
        for (int i = 0; i < v.size(); i++) {
            mean_diff_squared[i] = (double) w_normd[i] * pow(((double) v[i] - vector_mean), 2);
        }
        double SSE = accumulate(mean_diff_squared.begin(), mean_diff_squared.end(), 0.0);
        return sqrt( sum_w_squared / (sum_w_squared - sum_squared_w) * SSE );
    };

    double calcBinomialCdfUpper(int successes, double prob, int trials);

    double calcGaussianCdfUpper(double x, double sigma);

    double calcHygeoCdfUpper(unsigned int successes, unsigned int sub_population, unsigned int trials, unsigned int population);
    
    double calcTdistCdfLower(double t_stat, int df);
    
    //////////////////////////////////////////////////////////////////////////////
    //
    // process_mem_usage(double &, double &) - takes two doubles by reference,
    // attempts to read the system-dependent data for a process' virtual memory
    // size and resident set size, and return the results in MB.
    //
    // On failure, returns 0.0, 0.0
    // Modified from http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
    void procMemUsage(double& vm_usage, double& resident_set);


    typedef pair<int,double> int_double_pair;

    struct SecondCmp {
        bool operator()(const int_double_pair &lhs, const int_double_pair & rhs) {
            return lhs.second < rhs.second;
        }
    };

    typedef pair<unsigned int, unsigned int> Range;

    struct RangeCompare { //overlapping ranges are considered equivalent
        bool operator()(const Range& lhv, const Range& rhv) const {
            return lhv.second < rhv.first;
        }
    };

    inline bool inRange(const set<Range,RangeCompare>& ranges, unsigned int value) {
        return ranges.find(Range(value, value)) != ranges.end();
    }

    vector<size_t> benjaminiHochbergFDR (map<size_t,double> &pvals, double q_level);

    pair<double,double> numBinsFreedmanDiaconis (vector<int>);
}

#endif	/* _JMUTIL_HPP */

