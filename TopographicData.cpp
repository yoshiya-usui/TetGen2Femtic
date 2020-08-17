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
#include <algorithm>
#include <iomanip>
#include <math.h>
#include "TopographicData.h"

// Constructor
TopographicData::TopographicData():
	m_numTopographicData(0),
	m_topographicalData(NULL),
	m_maxDistanceForInterpolating(-1.0),
	m_maxNumOfPointsForInterpolating(-1)
{
};

// Destructor
TopographicData::~TopographicData(){

	if( m_topographicalData != NULL ){
		delete [] m_topographicalData;
		m_topographicalData = NULL;
	}

};

// Read topographical data
void TopographicData::readTopographicData( const std::string& fileName ){

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Read topographic data from " << fileName.c_str() << std::endl;

	ifs >> m_numTopographicData;
	std::cout << "Total number of topographic data : " << m_numTopographicData << std::endl;

	m_topographicalData = new TopographicData::XYZ[m_numTopographicData];

	for( int i = 0; i < m_numTopographicData; ++i ){
		ifs >> m_topographicalData[i].X;
		ifs >> m_topographicalData[i].Y;
		ifs >> m_topographicalData[i].Z;
		m_topographicalData[i].X *= 1000.0;
		m_topographicalData[i].Y *= 1000.0;
		//m_topographicalData[i].Z *= 1000.0;
	}

	ifs.close();

};

// Interpolate Z coordinate by inverse distance weighting
bool TopographicData::interpolateZCoord( const double xCoord, const double yCoord, double& zCoord ) const{

	std::vector< std::pair<double,int> > stackDistances;
	for( int i = 0; i < m_maxNumOfPointsForInterpolating; ++i ){
		stackDistances.push_back( std::make_pair( 1.0e12, -1 ) );// Initialize 
	}

	for( int iData = 0; iData < m_numTopographicData; ++iData ){

		const double distance = hypot( xCoord - m_topographicalData[iData].X, yCoord - m_topographicalData[iData].Y );

		if( distance > m_maxDistanceForInterpolating ){
			continue;
		}

		if( distance < stackDistances.back().first ){

			stackDistances.back().first = distance;
			stackDistances.back().second = iData;

			sort( stackDistances.begin(), stackDistances.end() );
		}

	}

	if( stackDistances.front().second < 0 ){
		//std::cerr << " Warning : No topology data were found near the point (X,Y) = ( " << xCoord << ", " << yCoord << " )" << std::endl;
		return false;
	}

	const double threthold = 1.0e-6;

	zCoord = 0.0;
	double weightSum(0.0);
	for( std::vector< std::pair<double,int> >::iterator itr = stackDistances.begin();
		itr != stackDistances.end(); ++itr ){

		if( itr->second < 0 ){
			continue;
		}

		const double weight = 1.0 / ( threthold + itr->first );
		zCoord += weight * m_topographicalData[itr->second].Z;
		weightSum += weight;

	}
	zCoord /= weightSum;

	return true;

}

// Set maximum distance used for interpolating altitudes
void TopographicData::setMaxDistanceForInterpolating( const double dist ){
	m_maxDistanceForInterpolating = dist;
}

// Set maximum number of points used for interpolating altitudes
void TopographicData::setMaxNumOfPointsForInterpolating( const int num ){
	m_maxNumOfPointsForInterpolating = num;
}
