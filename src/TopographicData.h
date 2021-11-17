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
#ifndef DBLDEF_TOPOGRAPHIC_DATA
#define DBLDEF_TOPOGRAPHIC_DATA

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

// Class of topographic data
class TopographicData{

public:
	TopographicData();

	~TopographicData();

	// Read topographical data
	void readTopographicData( const std::string& fileName );

	// Interpolate Z coordinate by inverse distance weighting
	bool interpolateZCoord( const double xCoord, const double yCoord, double& zCoord ) const;

	// Set maximum distance used for interpolating altitudes
	void setMaxDistanceForInterpolating( const double dist );

	// Set maximum number of points used for interpolating altitudes
	void setMaxNumOfPointsForInterpolating( const int num );

private:
	struct XYZ{
		double X;
		double Y;
		double Z;
	};

	// Name of the file containing the topographical data
	std::string m_fileName;

	// Total number of topographical data
	int m_numTopographicData;

	// Topographical data
	TopographicData::XYZ* m_topographicalData;

	// Maximum distance used for interpolating altitudes
	double m_maxDistanceForInterpolating;

	// Maximum number of points used for interpolating altitudes
	int m_maxNumOfPointsForInterpolating;

};

#endif
