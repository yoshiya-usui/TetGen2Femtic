//--------------------------------------------------------------------------
// MIT License
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//--------------------------------------------------------------------------
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
