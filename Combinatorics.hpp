/* 
 * File:   Combinatorics.hpp
 * Author: jmintser
 *
 * Created on June 2, 2011, 6:06 PM
 */

#ifndef _COMBINATORICS_HPP
#define	_COMBINATORICS_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <numeric> // for accumulate
#include <functional> // for multiplies

using namespace std;

class Combinatorics {
public:


    /* Implements the combinations generator described in
     *  Mifsud: Algorithm 154: combination in lexicographical order. Commun. ACM 6(3): 103 (1963)
     *
     * Mifsud's algorithm is allegedly the fastest according to
     *  Akl, S.G. A Comparison of Combination Generation Methods. ACM TOMS: 7(1) (1981)
     *
     * If want Cn,r make sure vector I is of length r
     * This is modified from previous version so that a global bool is no longer required
     * Must be called with bool last = true and vector<int> I(k). Run until last == true again.  I ints start with 1
     */
    static inline void next_combination(int n, vector<int> &I, bool &last) {
        int r = I.size();
        if (last) {
            for (int j = 0; j < r; j++) {
                I[j] = j + 1;
                last = false; // having this inside the loop provides a cheap implicit check for r == 0
            }
            //last = false;
            return;
        }
        if (I[r - 1] < n) {
            I[r - 1]++;
            return;
        }
        for (int j = r; j >= 2; j--) {
            //cout << "j: " << j << endl; ////
            if (I[j - 2] < (n - r + j - 1)) {
                I[j - 2]++;
                for (int s = j; s <= r; s++) {
                    //cout << "s: " << s << endl; ////
                    I[s - 1] = I[j - 2] + s - (j - 1);
                }
                return;
            }
        }
        last = true;
        return;
    }

};

#endif	/* _COMBINATORICS_HPP */

