/*
   File:        randomization.c
   Author:      Artur Dubrawski
   Created:     November 1999
   Description: Collection of randomization test related tools

   Copyright 1999, Schenley Park Research Inc
*/

#include "standard.h"
#include "randomization.h"
#include "amstr.h"

//---------------------------------------
static bool rand_legal_type( int type )
{
  if( type > 0 && type < RTEST_NUMTYPES )
    return(TRUE);
  else
    return(FALSE);
}

//-------------------------------------------
double rand_sig_test( int type, double val, double rmean, double rsdev )
{
  double result = -1.0;
  if( rand_legal_type(type) )
  {
    if( rsdev > AWD_ZEROTOL )
    {
      double z = fabs((rmean-val))/rsdev;
      result = 0.0;
      if( (type == RTEST_HI && rmean < val) || ((type==RTEST_LO && rmean > val) && (z > 1.645)))
      {
        if( z > 2.33 )
        {
          if( z > 2.58 )
          {
            if( z > 2.88 )
              result = 0.998; // z>2.88
            else
              result = 0.995;
          }
          else
            result = 0.99; // z>2.58
        }
        else
          result = 0.95; // z>2.33
      }
      
      if( type == RTEST_TWO && (z > 1.96) )
      {
        if( z > 2.58 )
        {
          if( z > 2.81 )
          {
            if( z > 3.08 )
              result = 0.998; // z>2.88
            else
              result = 0.995;
          }
          else
            result = 0.99; // z>2.58
        }
        else
          result = 0.95; // z>2.33
      }
    }
    else // sig is very small (an unusual case)
    {
      if( fabs(rmean-val) > AWD_ZEROTOL )
        result = 0.95;
    }
  }
  return(result);
}

char *mk_description_of_test(double test)
{
  char *s;
  if(test<.95) s = mk_copy_string("not significant (less than 95%)");
  else if(test<.99) s = mk_copy_string("SIGNIFICANT (95%-99%)");
  else if(test<.995) s = mk_copy_string("QUITE SIGNIFICANT (99%-99.5%)");
  else if(test<.998) s = mk_copy_string("VERY SIGNIFICANT (99.5%-99.8%)");
  else s = mk_copy_string("EXTREMELY SIGNIFICANT (99.8%-100%)");
  return s;
}
