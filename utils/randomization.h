/*
   File:        randomization.h
   Author:      Artur Dubrawski
   Created:     November 1999
   Description: Collection of randomization test related tools

   Copyright 1999, Schenley Park Research Inc
*/

#ifndef RANDOMIZATION_H
#define RANDOMIZATION_H

#include "amdym.h"

#ifndef AWD_ZEROTOL 
#define AWD_ZEROTOL 1.0e-12
#endif

#define RTEST_NUMTYPES 3
#define RTEST_NONE -1
#define RTEST_TWO 0
#define RTEST_HI 1
#define RTEST_LO 2

//----------------------------------------------------------------------
// This is a general purpose significance test function.
// It tells how significantly 'val' differs from 'rmean',
// where 'rmean' and 'rsdev' are respectively mean and std. deviation of some
// reference distribution (presumably normal).
// 'type' determines the type of the test: put RTEST_TWO for the two-tailed,
// RTEST_HI for the one-tailed upper tail test (how significantly 'val' is larger 
// than 'rmean'), and RTEST_LO for the respective one-tailed lower tail test.
// Returned values:
// -1.0   - test failed (illegal 'type')
//  0.0   - no significance (below 95.0% level)
//  0.95  - sig level of at least 95.0%
//  0.99  - sig level of at least 99.0%
//  0.995 - sig level of at least 99.5%
//  0.998 - sig level of at least 99.8%
//
extern double rand_sig_test( int type, double val, double rmean, double rsdev );

char *mk_description_of_test(double test);

#endif /* RANDOMIZATION_H */

