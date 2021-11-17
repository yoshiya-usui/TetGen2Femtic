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
