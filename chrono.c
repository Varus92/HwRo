#include <time.h>
#include <stdlib.h>


#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

double myWallTime()
{
#ifdef __APPLE__
	static double timeConvert = 0.0;
	if ( timeConvert == 0.0 )
	{
		mach_timebase_info_data_t timeBase;
		mach_timebase_info(&timeBase);
		timeConvert = (double)timeBase.numer / (double)timeBase.denom / 1000000000.0;
	}
	return mach_absolute_time() * timeConvert;
#else
	
	return 0.0;
#endif // __APPLE__
}

double second()
{
#if 1
	double t = myWallTime();
	return(t);
#else
	return ((double)clock()/(double)CLK_TCK);
#endif
}
