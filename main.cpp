//-------------------------------------------------------------------------------------------------------
// Copyright 2020 Yoshiya Usui
//
// This file is part of TetGen2Femtic.
//
// TetGen2Femtic is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TetGen2Femtic is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TetGen2Femtic. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------------------------------------
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "MeshData.h"

int main( int argc, char* argv[] ){

	if( argc < 2 ){
		std::cerr << "You must specify root name of tetgen file" << std::endl;
		exit(1);
	}

	const std::string fileNameRoot = argv[1];

	bool includeSediment(false);
	bool rotateModel(false);
	bool divAllElements(false);

	int numThread(1);

	for( int i = 2; i < argc; ++i ){
		if( strcmp(argv[i], "-sediment") == 0 ){
			includeSediment = true;
		}
		else if( strcmp(argv[i], "-rotate") == 0 ){
			rotateModel = true;
		}
		else if( strcmp(argv[i], "-thread") == 0 ){
			if( argc <= i + 1 ){
				std::cerr << "You must write total number of threads after -thread option." << std::endl;
				exit(1);
			}
			numThread = atoi(argv[++i]);
		}
		else if( strcmp(argv[i], "-div_all") == 0 ){
			divAllElements = true;
		}
	}

#ifdef _DEBUG_WRITE
	std::cout << "includeSeaConductivityAtExtendedRegion : " << includeSeaConductivityAtExtendedRegion << std::endl;
#endif

	MeshData::getInstance()->makeMeshDataForFemtic( fileNameRoot, includeSediment, rotateModel, numThread, divAllElements );

	return 0;

}
