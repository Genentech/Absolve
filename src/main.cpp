
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <string>


#include <boost/filesystem/operations.hpp>

#include <boost/filesystem/path.hpp>

#include <iostream>



#include "pipeline.h"

/////////////////////////////////////////////////////////

extern void TimeAligner();

namespace fs = boost::filesystem;

int main(int argc, char **argv, char **env){

	fs::path full_path( fs::initial_path<fs::path>() );

	full_path = fs::system_complete( fs::path( argv[0] ) );

	std::string exepath = full_path.parent_path().string();
	AbsolvePipeline pipeline;
	pipeline.ProcessOptions(argc, argv, env, exepath.c_str());
	pipeline.InitPipeline();
	pipeline.Run();
	return 0;
}
////////////////////////////////////////////////////////////////////////





 
