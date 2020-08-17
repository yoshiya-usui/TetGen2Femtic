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
