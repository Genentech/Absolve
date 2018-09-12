
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	A collection of utility functions for timing/profiling and 
//	string manipulation
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include<utility>
#include<algorithm>
#include<string>

#include <boost/date_time/posix_time/posix_time.hpp>


#define ABS_MAX_AA_LEN  21845
#define ABS_MAX_NT_LEN  65536


typedef boost::posix_time::ptime TimePoint;

void ToUpper(char *str, int len);
void StringToUpper(std::string & data);

TimePoint GetCurTime();

int GetSecondsDiffToNow(TimePoint & sometimeago);

 
