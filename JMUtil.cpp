#include "JMUtil.hpp"

double JMUtil::calcBinomialCdfUpper(int successes, double prob, int trials) {
    if (successes > 0) {
        // using the cdf is about 10x faster than summing pdf's
        return gsl_cdf_binomial_Q(successes - 1, prob, trials);
    } else {
        // This is a little hack for the current bug in gsl
        return 1.0;
    }
}

//double JMUtil::calcBinomialCdfUpper(int successes, double prob, int trials) {
//    if (successes > 0) {
//        boost::math::binomial binom_distr(trials, prob);
//        return boost::math::cdf(boost::math::complement(binom_distr, successes - 1));
//    } else {
//        return 1.0;
//    }
//}

double JMUtil::calcHygeoCdfUpper(unsigned int successes, unsigned int sub_population, unsigned int trials, unsigned int population) {
    boost::math::hypergeometric hygeo(sub_population, trials, population);
    return boost::math::cdf(boost::math::complement(hygeo, successes));
}

double JMUtil::calcGaussianCdfUpper(double x, double sigma) {
    boost::math::normal gaussian(0, sigma);
    return boost::math::cdf(boost::math::complement(gaussian, x));
}

double JMUtil::calcTdistCdfLower(double t_stat, int df) {
    boost::math::students_t t_dist(df);
    return boost::math::cdf(t_dist, t_stat);
}

void JMUtil::procMemUsage(double& vm_usage, double& resident_set) {
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0 / 1024.0;
   resident_set = rss * page_size_kb / 1024.0;
   return;
}

vector<size_t> JMUtil::benjaminiHochbergFDR(map<size_t,double> &pvals, double q_level) {
    vector<int_double_pair> pval_pair_vec(pvals.begin(), pvals.end());
    sort(pval_pair_vec.begin(), pval_pair_vec.end(), SecondCmp());
    vector<size_t> fdr_passed_idx;
    for (size_t i = 0; i < pval_pair_vec.size(); i++) {
        double q = q_level * (i+1) / pval_pair_vec.size();
        if (pval_pair_vec[i].second <= q) {
            fdr_passed_idx.push_back(pval_pair_vec[i].first);
        }
    }
    sort(fdr_passed_idx.begin(), fdr_passed_idx.end());
    return fdr_passed_idx;
}

// Adapted from the MATLAB implementation /MATLAB/toolbox/stats/private/dfhistbins.m
pair<double,double> JMUtil::numBinsFreedmanDiaconis (const vector<int> sorted_vec) {
    double machine_epsilon = pow (2.0, -52);

    size_t n_obs = sorted_vec.size();
    double x_min = *sorted_vec.begin();
    double x_max = *sorted_vec.rbegin();
    double x_range = x_max - x_min;
    
    std::vector<double> double_vec(sorted_vec.begin(), sorted_vec.end());
    int quartile1 = roundTo( quantile(double_vec, 0.25), 0);
    int quartile2 = roundTo( quantile(double_vec, 0.75), 0);
    //cout << "quartile1: " << quartile1 << "; quartile2: " << quartile2 << endl;
    //exit(0);
    double iqr = quartile2 - quartile1;
    /* Guard against too small an IQR.  This may be because most observations are
     * censored, or because there are some extreme outliers.
     */
    if (iqr < x_range/10) {
        iqr = x_range/10;
    }

    double raw_bin_width = 2*iqr/pow((double) n_obs, 1.0/3.0);
//    cout << iqr << "\t" << n_obs << "\t" << x_range << "\t" << raw_bin_width << endl; ////
    double x_scale = max(fabs( x_max ), fabs( x_min ));
    // Make sure the bin width is not effectively zero
    raw_bin_width = max(raw_bin_width, machine_epsilon * x_scale);
//    cout << iqr << "\t" << n_obs << "\t" << x_range << "\t" << raw_bin_width << endl; ////

    double left_edge, right_edge, bin_width;
    double n_bins_actual = 1;
    if (x_range > max(sqrt(machine_epsilon)*x_scale, numeric_limits<double>::min())) {
        double pow_of_ten = pow(10.0, floor(log10(raw_bin_width))); // next lower power of 10
        double rel_size = raw_bin_width / pow_of_ten;
        if (rel_size < 1.5)
            bin_width = 1.0 * pow_of_ten;
        else if (rel_size < 2.5)
            bin_width = 2.0 * pow_of_ten;
        else if (rel_size < 4.0)
            bin_width = 3.0 * pow_of_ten;
        else if (rel_size < 7.5)
            bin_width = 5.0 * pow_of_ten;
        else
            bin_width = 10.0 * pow_of_ten;

        left_edge = min(bin_width*floor(x_min/bin_width), x_min);
        n_bins_actual = max(2.0, ceil((x_max - left_edge) / bin_width));
        right_edge = max(left_edge + n_bins_actual * bin_width, x_max);
    }
    //cout << "JMUtil::numBinsFreedmanDiaconis reporting " << n_bins_actual << " bins of width " << bin_width << endl; ////
    return pair<double,double>(n_bins_actual, bin_width);
}
