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
#include "MeshData.h"

#include <assert.h>

const double MeshData::convkm2m = 1.0e3;

const double MeshData::PI = 3.14159265359;

const double MeshData::DEG2RAD = PI / 180.0;

// Return the instance of the class
MeshData* MeshData::getInstance()
{
   	static MeshData instance;// The only instance
  	return &instance;
}

// Make mesh data for femtic
void MeshData::makeMeshDataForFemtic( const std::string& rootName, const bool includeSediments, const bool rotateModel,
									 const int numThread, const bool divAllElements ){

	std::cout << "Total number of threas is " << numThread << std::endl;

	readNodeData( rootName, rotateModel);
	readElemData( rootName );
	readNeighElements( rootName );
	readFaceData( rootName );

#if 0
	readBoundaryAttribute();
#endif
	m_doesUseBoundaryMaker = false;

	readResisitivity();

	if( includeSediments ){
		std::cout << "Sedimentary layer is included." << std::endl;
		includeSedimentaryLayers();
	}
	if( rotateModel ){
		std::cout << "Model is rotated at an angle of 90 degrees around the z axis." << std::endl;
	}
	
	if( m_doesUseBoundaryMaker ){
		calcBoundaryFacesFromBoundaryMarker();
	}
	else{
		calcBoundaryFaces();
	}

#ifdef _ADDITIONAL_FIXED_AREA
	if( m_useAdditionalFixedRegion ){
		// Write faces above the additional area in which resistivity values are fixed
		writeFacesAboveAdditionalFixedRegionToVTK(rootName);

		// Fix resistivity in the additional area in which resistivity values are fixed
		fixResistivityInAdditionalFixedRegion();
	}
#endif

	writeMeshData();

	makeResistivityBlock(divAllElements);

	writeResisitivityData();

	writeMeshDataToVTK( rootName );
}

// Default constructer
MeshData::MeshData():
	m_doesUseBoundaryMaker(true),
	m_numNodes(0),
	m_nodeCoords(NULL),
	m_rotationAngle(0.0),
	m_elements(NULL),
	m_node2elem(NULL),
	m_faceTetGen(NULL),
	//m_faceFemtic(NULL),
	m_regionAttrAirLayer(0),
	m_elem2blk(NULL),
	m_numSphere(0),
	m_radius(NULL),
	m_maxBlockLength(NULL),
	m_oblatenessHorizontal(NULL),
	m_oblateness(NULL),
	m_numObsSite(0),
	m_obsSite(NULL),
#ifdef _ADDITIONAL_FIXED_AREA
	m_numThreads(1),
	m_useAdditionalFixedRegion(false)
#else
	m_numThreads(1)
#endif
{
	m_centerCoordSpheres.X = 0.0;
	m_centerCoordSpheres.Y = 0.0;
	m_centerCoordSpheres.Z = 0.0;
}

// Destructer
MeshData::~MeshData()
{

	if( m_nodeCoords != NULL ){
		delete [] m_nodeCoords;
		m_nodeCoords = NULL;
	}

	if( m_elements != NULL ){
		delete [] m_elements;
		m_elements = NULL;
	}

	if( m_node2elem != NULL ){
		delete [] m_node2elem;
		m_node2elem = NULL;
	}

	if( m_faceTetGen != NULL ){
		delete [] m_faceTetGen;
		m_faceTetGen = NULL;
	}

	if( m_elem2blk != NULL ){
		delete [] m_elem2blk;
		m_elem2blk = NULL;
	}
	
	if( m_radius != NULL ){
		delete [] m_radius;
		m_radius = NULL;
	}

	if( m_maxBlockLength != NULL ){
		delete [] m_maxBlockLength;
		m_maxBlockLength = NULL;
	}
	
	if( m_oblatenessHorizontal != NULL ){
		delete [] m_oblatenessHorizontal;
		m_oblatenessHorizontal = NULL;
	}

	if( m_oblateness != NULL ){
		delete [] m_oblateness;
		m_oblateness = NULL;
	}

	if( m_obsSite != NULL ){
		delete [] m_obsSite;
		m_obsSite = NULL;
	}

}

// Read node data
void MeshData::readNodeData( const std::string& rootName, const bool rotateModel ){

	std::string fileName = rootName;
	fileName += ".node";

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Read node data from " << fileName.c_str() << std::endl;
	
	ifs >> m_numNodes;
	std::cout << "Total number of nodes : " << m_numNodes << std::endl;
	
	m_nodeCoords = new MeshData::XYZ[m_numNodes];

	int ibuf(0);
	ifs >> ibuf;
	if( ibuf != 3 ){
		std::cerr << "Dimension must be 3 !!" << std::endl;
		exit(1);
	}

	int iAttribute(0);
	ifs >> iAttribute;
	if( iAttribute != 0 ){
		std::cerr << "Attributes of nodes are ignored in mesh of femtic." << std::endl;
	}

	int iBoundaryMarker(0);
	ifs >> iBoundaryMarker;
	if( iBoundaryMarker != 0 ){
		std::cerr << "Boundary marker of nodes are ignored in mesh of femtic." << std::endl;
	}
	
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){

		int nodeID(0);
		ifs >> nodeID;
		if( iNode + 1 != nodeID ){
			std::cerr << "Node ID must be sequence number from 1 !! nodeID = " << nodeID << std::endl;
			exit(1);
		}

		ifs >> m_nodeCoords[iNode].X >> m_nodeCoords[iNode].Y >> m_nodeCoords[iNode].Z;

		if( rotateModel ){
			const double xOrg = m_nodeCoords[iNode].X;
			const double yOrg = m_nodeCoords[iNode].Y;
			m_nodeCoords[iNode].X = yOrg;
			m_nodeCoords[iNode].Y = -xOrg;
		}

		m_nodeCoords[iNode].X *= MeshData::convkm2m;
		m_nodeCoords[iNode].Y *= MeshData::convkm2m;
		m_nodeCoords[iNode].Z *= MeshData::convkm2m;

		if( iAttribute != 0 ){
			ifs >> ibuf;
		}

		if( iBoundaryMarker != 0 ){
			ifs >> ibuf;
		}

	}

	ifs.close();

}

// Read element data
void MeshData::readElemData( const std::string& rootName ){

	std::string fileName = rootName;
	fileName += ".ele";

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Read element data from " << fileName.c_str() << std::endl;

	ifs >> m_numElements;
	std::cout << "Total number of elements : " << m_numElements << std::endl;
	
	m_elements = new MeshData::Element[m_numElements];

	int nodesPerTet(0);
	ifs >> nodesPerTet;

	if( nodesPerTet != 4 ){
		std::cerr << "Nodes number of a tetrahedron must be 4 !!" << std::endl;
		exit(1);
	}

	int iAttribute(0);
	ifs >> iAttribute;
	if( iAttribute != 1 ){
		std::cerr << "You must specify region attribute to elements." << std::endl;
		exit(1);
	}

	m_node2elem = new std::vector<int>[m_numNodes];

	for( int iElem = 0; iElem < m_numElements; ++iElem ){

		int elemID(0);
		ifs >> elemID;

		if( iElem + 1 != elemID ){
			std::cerr << "Element ID must be sequence number from 1!!" << std::endl;
			exit(1);
		}

		for( int i = 0; i < 4; ++i ){
			int nodeID(0);
			ifs >> nodeID;
			m_elements[iElem].nodeID[i] = nodeID;
			m_node2elem[nodeID-1].push_back( elemID );
		}

		ifs >> m_elements[iElem].attribute;

	}

	ifs.close();

}

// Read neighbor elements
void MeshData::readNeighElements( const std::string& rootName ){

	std::string fileName = rootName;
	fileName += ".neigh";

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Read neighboors of elements from " << fileName.c_str() << std::endl;

	int nElem(0);
	ifs >> nElem;
	if( nElem != m_numElements ){
		std::cerr << "Total element number written in .neigh file (" << nElem << " ) is different from the one written in .ele file ( " << m_numElements << " ) !!" << std::endl;
		exit(1);
	}

	if( m_elements == NULL ){
		std::cerr << "m_elements is NULL !!" << std::endl;
		exit(1);
	}

	int nodesPerTet(0);
	ifs >> nodesPerTet;

	if( nodesPerTet != 4 ){
		std::cerr << "Nodes number of a tetrahedron must be 4 !!" << std::endl;
		exit(1);
	}

	for( int iElem = 0; iElem < m_numElements; ++iElem ){

		int elemID(0);
		ifs >> elemID;
		if( iElem + 1 != elemID ){
			std::cerr << "Element ID must be sequence number from 1!!" << std::endl;
			exit(1);
		}

		for( int i = 0; i < 4; ++i ){
			ifs >> m_elements[iElem].neigh[i];
		}
		
	}

	ifs.close();

}

// Read face data
void MeshData::readFaceData( const std::string& rootName ){

	std::string fileName = rootName;
	fileName += ".face";

	std::ifstream ifs( fileName.c_str() );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Read face data from " << fileName.c_str() << std::endl;

	ifs >> m_numFaces;
	std::cout << "Total number of faces : " << m_numFaces << std::endl;
	
	m_faceTetGen = new MeshData::FaceTetGen[m_numFaces];

	int iBoundaryMarker(0);
	ifs >> iBoundaryMarker;
	if( iBoundaryMarker != 1 ){
		std::cerr << "Boundary marker of face must be one." << std::endl;
		exit(1);
	}

	for( int iFace = 0; iFace < m_numFaces; ++iFace ){

		int faceID(0);
		ifs >> faceID;

		for( int i = 0; i < 3; ++i ){
			ifs >> m_faceTetGen[iFace].nodeID[i];
		}
		ifs >> m_faceTetGen[iFace].attribute;
	}
	
	ifs.close();

}

// Read releationship between boundary of the model and attributes of faces
void MeshData::readBoundaryAttribute(){

	std::ifstream ifs( "bound_attr.dat" );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file bound_attr.dat !!" << std::endl;
		exit(1);
	}

	std::cout << "Read relationship between boundary and face from bound_attr.dat" << std::endl;

	int nBound(0);
	ifs >> nBound;

	if( nBound < 0 ){
		m_doesUseBoundaryMaker= false;
		std::cout << "Calculate boundary faces without boundary maker because the number of boundary is negative in bound_attr.dat" << std::endl;
		return;
	}

	m_doesUseBoundaryMaker= true;

	if( nBound != NUM_BOUNDARY ){
		std::cerr << "Number of boundary must be " << NUM_BOUNDARY << " !!" <<std::endl;
		exit(1);
	}

	for( int iBoun = 0; iBoun < NUM_BOUNDARY; ++iBoun ){

		int numAttr(0);
		ifs >> numAttr;

		for( int iAttr = 0; iAttr < numAttr; ++iAttr ){
			int attr(0);
			ifs >> attr;
			m_boundaryFaces[iBoun].attributes.push_back( attr );
		}

	}

	ifs.close();

}

// Read releationship between resistivity values and region attributes
void MeshData::readResisitivity(){

	std::ifstream ifs( "resistivity_attr.dat" );

	if( !ifs.is_open() ){
		std::cerr << "Cannot open file resistivity_attr.dat !!" << std::endl;
		exit(1);
	}

	std::cout << "Read relationship between resisitivity and region from resistivity_attr.dat." << std::endl;

	int numResistivity(0);
	ifs >> numResistivity;

	std::cout << "Number of different resisitivities : " << numResistivity << std::endl;

	if( numResistivity < 1 ){
		std::cerr << "Number of different resisitivities must be greater than or equal to one !!" << std::endl;
		exit(1);
	}

	m_resistivityData.reserve(numResistivity);

	double resistivityAir(0.0);
	for( int iRes = 0; iRes < numResistivity; ++iRes ){

		int attr(0);
		double resisitivity(0.0);
		ifs >> attr;
		m_regionAttr2ResistivityID.insert( std::make_pair( attr, iRes ) );

		ResistivityData resistivityDataBuf;
		resistivityDataBuf.attr = attr;

		int ibuf;
		ifs >> resistivityDataBuf.resistivity >> resistivityDataBuf.ndiv >> ibuf;
		
		if( iRes == 0 ){
			if( resistivityDataBuf.ndiv >= 0 ){
				std::cerr << "Number of division must be negative for the first resistivity data ( the air layer ) !!" << std::endl;
				exit(1);
			}
			m_regionAttrAirLayer = resistivityDataBuf.attr;
			resistivityAir = resistivityDataBuf.resistivity;
		}

		resistivityDataBuf.fix = ( ibuf > 0 ? true : false );

		if( resistivityDataBuf.ndiv < 0 ){
			std::cout << "Resistivity data " << iRes << " is the air or the sea because its division number is negative." << std::endl;
			m_regionAttrAirOrSea.insert( resistivityDataBuf.attr );
		}

		m_resistivityData.push_back(resistivityDataBuf);
	}

	std::cout << "First resistivity data must be the air layer." << std::endl;
	std::cout << "Therefore, the resisitivity of the air is set to be " << resistivityAir << " [Ohm-m]" << std::endl;

	//*******************************************************************************
	//***** Read data of spheres specifing maximum length of resistivity blocks *****
	//*******************************************************************************
	ifs >> m_centerCoordSpheres.X >> m_centerCoordSpheres.Y >> m_centerCoordSpheres.Z;
	m_centerCoordSpheres.X *= MeshData::convkm2m;
	m_centerCoordSpheres.Y *= MeshData::convkm2m;
	m_centerCoordSpheres.Z *= MeshData::convkm2m;
	std::cout << "Center coorinate of the spheres [m] : (X,Y,Z) = (" << m_centerCoordSpheres.X << ", " << m_centerCoordSpheres.Y <<  ", " << m_centerCoordSpheres.Z << ")" << std::endl;

	ifs >> m_rotationAngle;
	std::cout << "Rotation angle of horizontal ellipsoid [deg] : " << m_rotationAngle << std::endl;
	m_rotationAngle *= DEG2RAD;

	ifs >> m_numSphere;
	std::cout << "Total number of spheres : " << m_numSphere << std::endl;

	if( m_numSphere > 0 ){
		m_radius = new double[m_numSphere];
		m_maxBlockLength = new double[m_numSphere];
		m_oblatenessHorizontal = new double[m_numSphere];
		m_oblateness = new double[m_numSphere];
		std::cout << "<Radius[m]> <Edge Length[m]> <OblatenessHorizontal> <Oblateness>" << std::endl;
	}

	for( int i = 0; i < m_numSphere; ++i ){
		ifs >> m_radius[i] >> m_maxBlockLength[i] >> m_oblatenessHorizontal[i] >> m_oblateness[i];
		m_radius[i] *= MeshData::convkm2m;
		m_maxBlockLength[i] *= MeshData::convkm2m;
		std::cout << std::setw(15) << std::scientific << m_radius[i];
		std::cout << std::setw(15) << std::scientific << m_maxBlockLength[i];
		std::cout << std::setw(15) << std::scientific << m_oblatenessHorizontal[i];
		std::cout << std::setw(15) << std::scientific << m_oblateness[i] << std::endl;
	}
	for( int i = 1; i < m_numSphere; ++i ){
		if( m_radius[i] < m_radius[i-1] ){
			std::cerr << "Radius of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
		if( m_maxBlockLength[i] < m_maxBlockLength[i-1] ){
			std::cerr << "Edge length of the region " << i << " is smaller than that of the previous region." << std::endl;
			exit(1);
		}
		if( m_oblatenessHorizontal[i] < 0 || m_oblatenessHorizontal[i] > 1 ){
			std::cerr << "Oblateness of horizontal ellipsoid must be smaller than 1 and larger than 0." << std::endl;
			exit(1);
		}
		if( m_oblateness[i] < 0 || m_oblateness[i] > 1 ){
			std::cerr << "Oblateness must be smaller than 1 and larger than 0." << std::endl;
			exit(1);
		}
		if( m_radius[i]*(1.0-m_oblatenessHorizontal[i]) < m_radius[i-1]*(1.0-m_oblatenessHorizontal[i-1] ) ){
			std::cerr << "Length of shorter axis of horizontal ellipsoid " << i << " is less than that of the previous ellipsoid." << std::endl;
			exit(1);
		}
		if( m_radius[i]*(1.0-m_oblateness[i]) < m_radius[i-1]*(1.0-m_oblateness[i-1] ) ){
			std::cerr << "Depth of sphere " << i << " is shallower than that of the previous sphere in the earth." << std::endl;
			exit(1);
		}
	}

	//******************************************************************
	//***** Read data of maximum length near the observation sites *****
	//******************************************************************
	ifs >> m_numObsSite;

	std::cout << "Total number of observed site : " << m_numObsSite << std::endl;	

	m_obsSite = new ObservedSite[m_numObsSite];

	for( int iObs = 0; iObs < m_numObsSite; ++iObs ){
		ifs >> m_obsSite[iObs].X >> m_obsSite[iObs].Y >> m_obsSite[iObs].Z;
		m_obsSite[iObs].X *= MeshData::convkm2m;
		m_obsSite[iObs].Y *= MeshData::convkm2m;
		m_obsSite[iObs].Z *= MeshData::convkm2m;

		std::cout << "Locations of site [m] :" << std::endl;
		std::cout << std::setw(15) << std::scientific << m_obsSite[iObs].X
					<< std::setw(15) << std::scientific << m_obsSite[iObs].Y
					<< std::setw(15) << std::scientific << m_obsSite[iObs].Z << std::endl;

		ifs >> m_obsSite[iObs].numCircle;
		m_obsSite[iObs].radius.reserve(m_obsSite[iObs].numCircle);
		m_obsSite[iObs].length.reserve(m_obsSite[iObs].numCircle);

		std::cout << std::setw(15) << m_obsSite[iObs].numCircle << std::endl;

		for( int i = 0; i < m_obsSite[iObs].numCircle; ++i ){
			double radius(-1.0);
			double length(-1.0);
			ifs >> radius;
			ifs >> length;
			radius *= MeshData::convkm2m;
			length *= MeshData::convkm2m;
			if( i > 0 && radius < m_obsSite[iObs].radius.back() ){
				std::cerr << " Error : Radius must be specified in ascending order !! : " << radius << std::endl;
				exit(1);
			}
			if( i > 0 && length < m_obsSite[iObs].length.back() ){
				std::cerr << " Error : Edge length must increase with radius !! : " << length << std::endl;
				exit(1);
			}
			m_obsSite[iObs].radius.push_back( radius );
			m_obsSite[iObs].length.push_back( length );

			std::cout << std::setw(15) << std::scientific << m_obsSite[iObs].radius[i]
						<< std::setw(15) << std::scientific << m_obsSite[iObs].length[i] << std::endl;
		}
	}

#ifdef _ADDITIONAL_FIXED_AREA
	//***************************************************************************
	//***** Read data for the region where the resisitivity values are fixed *****
	//***************************************************************************
	int ibuf(0);
	ifs >> ibuf;
	if( ibuf != 0 ){
		m_useAdditionalFixedRegion = true;
		ifs >> m_paramForAdditionalFixedRegion.centerCoord.X;
		ifs >> m_paramForAdditionalFixedRegion.centerCoord.Y;
		ifs >> m_paramForAdditionalFixedRegion.centerCoord.Z;
		ifs >> m_paramForAdditionalFixedRegion.xLength;
		ifs >> m_paramForAdditionalFixedRegion.yLength;
		ifs >> m_paramForAdditionalFixedRegion.zLength;
		ifs >> m_paramForAdditionalFixedRegion.rotationAngle;
		ifs >> m_paramForAdditionalFixedRegion.minDepthOfEarthSurface;
		ifs >> m_paramForAdditionalFixedRegion.maxDepthOfEarthSurface;
		ifs >> m_paramForAdditionalFixedRegion.resistivity;

		m_paramForAdditionalFixedRegion.centerCoord.X *= MeshData::convkm2m;
		m_paramForAdditionalFixedRegion.centerCoord.Y *= MeshData::convkm2m;
		m_paramForAdditionalFixedRegion.centerCoord.Z *= MeshData::convkm2m;
		m_paramForAdditionalFixedRegion.xLength *= MeshData::convkm2m;
		m_paramForAdditionalFixedRegion.yLength *= MeshData::convkm2m;
		m_paramForAdditionalFixedRegion.zLength *= MeshData::convkm2m;
		m_paramForAdditionalFixedRegion.minDepthOfEarthSurface *= MeshData::convkm2m;
		m_paramForAdditionalFixedRegion.maxDepthOfEarthSurface *= MeshData::convkm2m;

		std::cout << "Center coorinate of additional fixed region [m] : (X,Y,Z) = ("
			<< m_paramForAdditionalFixedRegion.centerCoord.X << ", "
			<< m_paramForAdditionalFixedRegion.centerCoord.Y << ", "
			<< m_paramForAdditionalFixedRegion.centerCoord.Z << ")" << std::endl;
		std::cout << "Length along x direction [m] : " << m_paramForAdditionalFixedRegion.xLength << std::endl;
		std::cout << "Length along y direction [m] : " << m_paramForAdditionalFixedRegion.yLength << std::endl;
		std::cout << "Length along z direction [m] : " << m_paramForAdditionalFixedRegion.zLength << std::endl;
		std::cout << "Rotation angle [deg.] : " << m_paramForAdditionalFixedRegion.rotationAngle << std::endl;
		std::cout << "Minimum depth of the earth's surface [m] : " << m_paramForAdditionalFixedRegion.minDepthOfEarthSurface << std::endl;
		std::cout << "Maximum depth of the earth's surface [m] : " << m_paramForAdditionalFixedRegion.maxDepthOfEarthSurface << std::endl;		
		std::cout << "Resistivity [Ohm-m] : " << m_paramForAdditionalFixedRegion.resistivity << std::endl;		
		m_paramForAdditionalFixedRegion.rotationAngle *= DEG2RAD;

		int attrMax = 0;
		for( std::vector<ResistivityData>::const_iterator itr = m_resistivityData.begin(); itr != m_resistivityData.end(); ++itr ){
			attrMax = std::max( itr->attr, attrMax );
		}
		m_paramForAdditionalFixedRegion.attribute = attrMax + 1;
		std::cout << "Attribute of additional fixed region : " << m_paramForAdditionalFixedRegion.attribute << std::endl;

		ResistivityData resistivityDataBuf;
		resistivityDataBuf.attr = m_paramForAdditionalFixedRegion.attribute;
		resistivityDataBuf.fix = true;
		resistivityDataBuf.ndiv = -1;
		resistivityDataBuf.resistivity = m_paramForAdditionalFixedRegion.resistivity;
		m_resistivityData.push_back(resistivityDataBuf);
		m_regionAttr2ResistivityID.insert( std::make_pair( m_paramForAdditionalFixedRegion.attribute, numResistivity ) );
	}
#endif

	ifs.close();

}

// Include sedimentary layers
void MeshData::includeSedimentaryLayers(){

	const std::string fileName = "sediment_layer.dat";

	//// Center coordinate of the area including predifined sediments
	//XYZ centerCoordSedimentArea = { 0.0, 0.0, 0.0 };

	//// Radius of the area including predifined sediments
	//double radiusSedimentArea(0.0);

	//// Depth of sedimentary layer
	//double depthSedimentLayer(0.0);

	//// Attribute of sediment
	//int attrSedimentLayer(-1);

	//// Resistivity of sediment
	//double resistivitySedimentLayer(0.0);

	SedimentParameter sedimentParam;

	// Number of the region attribute pair for specifing the area including predifined sediments
	std::vector< std::pair<int,int> > attrPair;

	readParametersOfSedimentArea( fileName, sedimentParam, attrPair );

	std::vector<int> sedimentFaceSerials;
	findTrianglesAboveSediment( attrPair, sedimentParam.center, sedimentParam.radius, sedimentParam.oblateness, sedimentParam.rotation, sedimentFaceSerials );

	std::vector<int> elemSerials;
	findElementsInSediment( sedimentFaceSerials, sedimentParam.center, sedimentParam.radius, sedimentParam.depth, sedimentParam.oblateness, sedimentParam.rotation, elemSerials );

	changeAttributeOfSediment( elemSerials, sedimentParam.attr, sedimentParam.resistivity );
}

// Read parameter file
void MeshData::readParametersOfSedimentArea( const std::string& fileName, SedimentParameter& param, std::vector< std::pair<int,int> >& attrPair ) const{

	std::ifstream ifs( fileName.c_str() );

	if( ifs.fail() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Read parameters from " << fileName.c_str() << std::endl;

	ifs >> param.center.X;
	ifs >> param.center.Y;
	std::cout << "Center coordinate of the area including predifined sediments [km] : (X,Y) = ( " << param.center.X << ", " << param.center.Y << " )" << std::endl;
	param.center.X *= MeshData::convkm2m;
	param.center.Y *= MeshData::convkm2m;

	ifs >> param.radius;
	std::cout << "Radius of the area including predifined sediments [km] : " << param.radius << std::endl;
	param.radius *= MeshData::convkm2m;

	ifs >> param.depth;
	std::cout << "Depth of the sedimentary layer [km] : " << param.depth << std::endl;
	param.depth *= MeshData::convkm2m;

	ifs >> param.oblateness;
	std::cout << "Oblateness of the area including predifined sediments : " << param.oblateness << std::endl;

	ifs >> param.rotation;
	std::cout << "Rotation angle of the area including predifined sediments [deg.] : " << param.rotation << std::endl;
	param.rotation *= DEG2RAD;

	ifs >> param.attr;
	std::cout << "Region attribute of the sedimentary layer : " << param.attr << std::endl;

	ifs >> param.resistivity;
	std::cout << "Resistivity of the sedimentary layer [Ohm-m] : " << param.resistivity << std::endl;

	int numPair(0);
	ifs >> numPair;
	std::cout << "Number of the region attribute pair for specifing the area including predifined sediments : " << numPair << std::endl;

	for( int i = 0; i < numPair; ++i ){
		std::pair<int,int> pair;
		ifs >> pair.first >> pair.second;
		attrPair.push_back(pair);
	}

	std::cout << "Region attribute pair for specifing the area including predifined sediments : " << std::endl;
	for( std::vector< std::pair<int,int> >::iterator itr = attrPair.begin(); itr != attrPair.end(); ++itr ){
		std::cout << itr->first << " " << itr->second << std::endl;
	}

	ifs.close();

}

// Find triangles above the sediment
void MeshData::findTrianglesAboveSediment( const std::vector< std::pair<int,int> >& attrPair, const XY& centerCoord, 
	const double radius, const double oblateness, const double rotation, std::vector<int>& sedimentFaceSerials ) const{

	for( int iFace = 0; iFace < m_numFaces; ++iFace ){

		int elemID(0);
		int faceID(0);
		calcElementFace( iFace, elemID, faceID );

		const int elemIDFromZero = elemID - 1;
		const int elemIDNeighFromZero = m_elements[elemIDFromZero].neigh[faceID] - 1;

		const int regionAttr = m_elements[elemIDFromZero].attribute;
		const int regionAttrNeib = m_elements[elemIDNeighFromZero].attribute;

		if( elemIDNeighFromZero >= 0 && std::find( attrPair.begin(), attrPair.end(), std::make_pair(regionAttr,regionAttrNeib) ) != attrPair.end() ||
			elemIDNeighFromZero >= 0 && std::find( attrPair.begin(), attrPair.end(), std::make_pair(regionAttrNeib,regionAttr) ) != attrPair.end() ){// Sediment
			if( locateInEllipse( centerCoord, radius, oblateness, rotation, calcGravityCenterOfTetGenFace(iFace) ) ){
				sedimentFaceSerials.push_back(iFace);
			}
		}
		
	}

}

// Decide the specified point locates in the ellipse
bool MeshData::locateInEllipse( const XY& centerCoord, const double radius, const double oblateness, const double rotation, const XY& coord ) const{

	const double vecXOrg = coord.X - centerCoord.X;
	const double vecYOrg = coord.Y - centerCoord.Y;
	const double vecX = vecXOrg * cos( - rotation ) - vecYOrg * sin( - rotation );
	const double vecY = vecXOrg * sin( - rotation ) + vecYOrg * cos( - rotation );

	const double longAxisLength = radius;
	const double shortAxisLength = longAxisLength * ( 1.0 - oblateness );

	double val = pow( vecX / longAxisLength, 2 ) + pow( vecY / shortAxisLength, 2 );

	if( val <= 1.0 ){
		return true;
	}

	return false;

}

// Decide the specified point locates in the ellipse
bool MeshData::locateInEllipse( const XY& centerCoord, const double radius, const double oblateness, const double rotation, const XYZ& coord ) const{
	const MeshData::XY coord2D = {coord.X, coord.Y};
	return locateInEllipse( centerCoord, radius, oblateness, rotation, coord2D );
}

// Find elements corresponding to sediment
void MeshData::findElementsInSediment( const std::vector<int>& sedimentFaceSerials, const XY& centerCoord,
		const double radius, const double depth, const double oblateness, const double rotation, std::vector<int>& elemSerials ) const{

	for( int iElem = 0; iElem < m_numElements; ++iElem ){
		const MeshData::XYZ coord = calcGravityCenter(iElem);
		if( !locateInEllipse( centerCoord, radius, oblateness, rotation, coord ) ){// Out of the ellipse
			continue;
		}
		for( std::vector<int>::const_iterator itr = sedimentFaceSerials.begin(); itr != sedimentFaceSerials.end(); ++itr ){
			double coordZOnFace(-1.0);
			if( calcZCoordOnTetGenFace(coord, *itr, coordZOnFace) ){
				if( coord.Z - coordZOnFace > 0.0 && coord.Z - coordZOnFace < depth ){
					elemSerials.push_back(iElem);
				}
				continue;
			}
		}
	}

}

// Change attributes of the elements corresponding to sediment layer
void MeshData::changeAttributeOfSediment( const std::vector<int>& elemSerials, const int attrSediment, const double resistivitySediment ){

	for( std::vector<int>::const_iterator itr = elemSerials.begin(); itr != elemSerials.end(); ++itr ){
		m_elements[*itr].attribute = attrSediment;
	}

	ResistivityData resistivityDataBuf;
	resistivityDataBuf.attr = attrSediment;
	resistivityDataBuf.fix = 1;
	resistivityDataBuf.ndiv = 0;
	resistivityDataBuf.resistivity = resistivitySediment;

	m_resistivityData.push_back(resistivityDataBuf);
	m_regionAttr2ResistivityID.insert( std::make_pair( attrSediment, static_cast<int>(m_resistivityData.size()) - 1 ) );

}

// Calculate gravity center of the specified face of Tetgen
MeshData::XY MeshData::calcGravityCenterOfTetGenFace( const int iFace ) const{

	if( iFace < 0 || iFace >= m_numFaces ){
		std::cerr << "iFace ( " << iFace << " ) is wrong !!" << std::endl;
		exit(1);
	}

	MeshData::XY coord = { 0.0, 0.0 };
	for( int i= 0; i < 3; ++i ){
		const int iNode = m_faceTetGen[iFace].nodeID[i] - 1;
		coord.X += m_nodeCoords[iNode].X;
		coord.Y += m_nodeCoords[iNode].Y;
	}
	coord.X /= 3.0;
	coord.Y /= 3.0;

	return coord;
}

// Calculate Z coordinate of the point below the specified face of Tetgen
//double MeshData::calcZCoordOnTetGenFace( const XYZ& point, const int iFace ) const{
//
//	assert( iFace < 0 || iFace >= m_numFaces );
//
//	const XYZ nodeCoords[3] = {
//		m_nodeCoords[ m_faceTetGen[iFace].nodeID[0] - 1 ],
//		m_nodeCoords[ m_faceTetGen[iFace].nodeID[1] - 1 ],
//		m_nodeCoords[ m_faceTetGen[iFace].nodeID[2] - 1 ]
//	};
//
//	const double areaTotal = calcAreaOnXYPlane( nodeCoords[0], nodeCoords[1], nodeCoords[2] );
//
//	const double areaCoords[3] = {
//		calcAreaOnXYPlane( point, nodeCoords[1], nodeCoords[2] ) / areaTotal,
//		calcAreaOnXYPlane( nodeCoords[0], point, nodeCoords[2] ) / areaTotal,
//		calcAreaOnXYPlane( nodeCoords[0], nodeCoords[1], point ) / areaTotal
//	};
//
//	const double EPS = 1.0e-12;
//	assert( areaCoords[0] < -EPS || areaCoords[0] < -EPS || areaCoords[0] < -EPS );
//
//	double val(0.0);
//	for( int i = 0; i < 3; ++i ){
//		val += nodeCoords[i].Z * areaCoords[i]; 
//	}
//	return val;
//
//}
bool MeshData::calcZCoordOnTetGenFace( const XYZ& point, const int iFace, double& zCoord ) const{

	assert( iFace >= 0 || iFace < m_numFaces );

	const XYZ nodeCoords[3] = {
		m_nodeCoords[ m_faceTetGen[iFace].nodeID[0] - 1 ],
		m_nodeCoords[ m_faceTetGen[iFace].nodeID[1] - 1 ],
		m_nodeCoords[ m_faceTetGen[iFace].nodeID[2] - 1 ]
	};

	const double areaTotal = calcAreaOnXYPlane( nodeCoords[0], nodeCoords[1], nodeCoords[2] );

	const double areaCoords[3] = {
		calcAreaOnXYPlane( point, nodeCoords[1], nodeCoords[2] ) / areaTotal,
		calcAreaOnXYPlane( nodeCoords[0], point, nodeCoords[2] ) / areaTotal,
		calcAreaOnXYPlane( nodeCoords[0], nodeCoords[1], point ) / areaTotal
	};

	const double EPS = 1.0e-12;
	if( fabs( 1.0 - areaCoords[0] - areaCoords[1] - areaCoords[2] ) > EPS ){// Out ot triangle
		return false;
	}	

	zCoord = 0.0;
	for( int i = 0; i < 3; ++i ){
		zCoord += nodeCoords[i].Z * areaCoords[i]; 
	}	
	return true;

}

// Calculate area of triangle on the XY plnae
double MeshData::calcAreaOnXYPlane( const MeshData::XYZ& point1, const MeshData::XYZ& point2, const MeshData::XYZ& point3 ) const{
	return 0.5 * fabs( ( point2.X - point1.X ) * ( point3.Y - point1.Y ) - ( point2.Y - point1.Y ) * ( point3.X - point1.X ) ); 
}

//// Determine if the inputed point locate in the specified face of Tetgen
//bool MeshData::locateInTriangleOfXYPlane( const XYZ& point, const int iFace ) const{
//
//	assert( iFace < 0 || iFace >= m_numFaces );
//
//	const XYZ nodeCoords[3] = {
//		m_nodeCoords[ m_faceTetGen[iFace].nodeID[0] - 1 ],
//		m_nodeCoords[ m_faceTetGen[iFace].nodeID[1] - 1 ],
//		m_nodeCoords[ m_faceTetGen[iFace].nodeID[2] - 1 ]
//	};
//
//	if( locateLeftOfSegmentOnXYPlane( point, nodeCoords[0], nodeCoords[1] ) &&
//		locateLeftOfSegmentOnXYPlane( point, nodeCoords[1], nodeCoords[2] ) &&
//		locateLeftOfSegmentOnXYPlane( point, nodeCoords[2], nodeCoords[0] ) ){
//		return true;
//	}
//
//	return false;
//}
//
//// Determine if the inputed point locate at the left of the segment on the X-Y plane
//bool MeshData::locateLeftOfSegmentOnXYPlane( const XYZ& point, const XYZ& startPointOfSegment, const XYZ& endPointOfSegment ) const{
//
//	if( ( endPointOfSegment.Y - startPointOfSegment.Y )*( point.X - startPointOfSegment.X ) >= ( endPointOfSegment.X - startPointOfSegment.X )*( point.Y - startPointOfSegment.Y ) ){
//		return true;
//	}
//
//	return false;
//}

// Calculate intersection point of a specified line and an X-Y plane
void MeshData::calcIntersectionPoint( const XYZ& startCoord, const XYZ& endCoord, const double zCoordOfPlane, XYZ& intersectionPoint ) const{

	assert( ( zCoordOfPlane >= startCoord.Z && zCoordOfPlane <= endCoord.Z ) || ( zCoordOfPlane >= endCoord.Z && zCoordOfPlane <= startCoord.Z ) );

	const double factor = ( zCoordOfPlane - startCoord.Z ) / ( endCoord.Z - startCoord.Z );

	//if( factor > 0.0 ){
	//	intersectionPoint.X = startCoord.X + factor * ( endCoord.X - startCoord.X );
	//	intersectionPoint.Y = startCoord.Y + factor * ( endCoord.Y - startCoord.Y );
	//	intersectionPoint.Z = startCoord.Z + factor * ( endCoord.Z - startCoord.Z );
	//}
	//else{
	//	intersectionPoint.X = endCoord.X + factor * ( endCoord.X - startCoord.X );
	//	intersectionPoint.Y = endCoord.Y + factor * ( endCoord.Y - startCoord.Y );
	//	intersectionPoint.Z = endCoord.Z + factor * ( endCoord.Z - startCoord.Z );
	//}
	intersectionPoint.X = startCoord.X + factor * ( endCoord.X - startCoord.X );
	intersectionPoint.Y = startCoord.Y + factor * ( endCoord.Y - startCoord.Y );
	intersectionPoint.Z = startCoord.Z + factor * ( endCoord.Z - startCoord.Z );

}

// Calculate volume of an element
double MeshData::calculateVolumeOfElement( const int iElem ) const{

	return calculateVolumeOfTetrahedron(
		m_nodeCoords[m_elements[iElem].nodeID[0]-1],
		m_nodeCoords[m_elements[iElem].nodeID[1]-1],
		m_nodeCoords[m_elements[iElem].nodeID[2]-1],
		m_nodeCoords[m_elements[iElem].nodeID[3]-1]
	);

}

// Calculate volume of tetrahedron 
double MeshData::calculateVolumeOfTetrahedron( const XYZ& point1, const XYZ& point2, const XYZ& point3, const XYZ& point4 ) const{

	const double val = ( point2.X*point3.Y* point4.Z + point2.Y*point3.Z*point4.X + point2.Z*point3.X*point4.Y - point2.Z*point3.Y*point4.X - point2.X*point3.Z*point4.Y - point2.Y*point3.X*point4.Z )
		             - ( point1.X*point3.Y* point4.Z + point1.Y*point3.Z*point4.X + point1.Z*point3.X*point4.Y - point1.Z*point3.Y*point4.X - point1.X*point3.Z*point4.Y - point1.Y*point3.X*point4.Z )
		             + ( point1.X*point2.Y* point4.Z + point1.Y*point2.Z*point4.X + point1.Z*point2.X*point4.Y - point1.Z*point2.Y*point4.X - point1.X*point2.Z*point4.Y - point1.Y*point2.X*point4.Z )
		             - ( point1.X*point2.Y* point3.Z + point1.Y*point2.Z*point3.X + point1.Z*point2.X*point3.Y - point1.Z*point2.Y*point3.X - point1.X*point2.Z*point3.Y - point1.Y*point2.X*point3.Z );

	return fabs(val / 6.0);

}

//// Calculate boundary faces
//void MeshData::calcBoundaryFaces(){
//
//	for( int iFace = 0; iFace < m_numFaces; ++iFace ){
//
//		bool doesBoundary(false);
//		const int boundType = calcBoundTypeFromAttribute( m_faceFemtic[iFace].attribute, doesBoundary );
//
//		if( doesBoundary ){
//			m_boundaryFaces[boundType].elemIDs.push_back( m_faceFemtic[iFace].elemID );
//			m_boundaryFaces[boundType].faceIDs.push_back( m_faceFemtic[iFace].faceID );
//		}
//
//	}
//
//}


// Calculate boundary faces from boundary marker
void MeshData::calcBoundaryFacesFromBoundaryMarker(){

	for( int iFace = 0; iFace < m_numFaces; ++iFace ){

		bool doesBoundary(false);
		const int boundType = calcBoundTypeFromAttribute( m_faceTetGen[iFace].attribute, doesBoundary );

		if( doesBoundary ){
			int elemID(0);
			int faceID(0);
			calcElementFace( boundType,	iFace, elemID, faceID );
			m_boundaryFaces[boundType].elemIDs.push_back( elemID );
			m_boundaryFaces[boundType].faceIDs.push_back( faceID );
		}

	}

}

// Calculate boundary faces
void MeshData::calcBoundaryFaces(){

	for( int iFace = 0; iFace < m_numFaces; ++iFace ){

		int elemID(0);
		int faceID(0);
		calcElementFace( iFace, elemID, faceID );

		const int elemIDFromZero = elemID - 1;
		const int elemIDNeighFromZero = m_elements[elemIDFromZero].neigh[faceID] - 1;

		if( elemIDNeighFromZero >= 0 ){

			if( m_regionAttrAirOrSea.find( m_elements[elemIDFromZero].attribute ) == m_regionAttrAirOrSea.end() &&
				m_regionAttrAirOrSea.find( m_elements[elemIDNeighFromZero].attribute ) != m_regionAttrAirOrSea.end() ){// Earth surface
				m_boundaryFaces[EARTH_SURFACE].elemIDs.push_back( elemID );
				m_boundaryFaces[EARTH_SURFACE].faceIDs.push_back( faceID );
#ifdef _ADDITIONAL_FIXED_AREA
				addIfFaceIndexAboveAdditionalFixedRegion(iFace);
#endif
			}
			else if( m_regionAttrAirOrSea.find( m_elements[elemIDFromZero].attribute ) != m_regionAttrAirOrSea.end() &&
					m_regionAttrAirOrSea.find( m_elements[elemIDNeighFromZero].attribute ) == m_regionAttrAirOrSea.end() ){// Earth surface
				int faceIDNeigh = -1;
				for( int i = 0; i < 4; ++i ){
					if( m_elements[elemIDNeighFromZero].neigh[i] == elemID ){
						faceIDNeigh = i;
						break;
					}
				}
				if( faceIDNeigh < 0 ){
					std::cerr << "Face ID of neighboor element cannot be found." << std::endl;
					exit(1);
				}
				m_boundaryFaces[EARTH_SURFACE].elemIDs.push_back( elemIDNeighFromZero + 1 );
				m_boundaryFaces[EARTH_SURFACE].faceIDs.push_back( faceIDNeigh );
#ifdef _ADDITIONAL_FIXED_AREA
				addIfFaceIndexAboveAdditionalFixedRegion(iFace);
#endif
			}
			//else{
			//	std::cerr << "Face " << faceID << " of element " << elemID << " is not a boundary face." << std::endl;
			//	exit(1);
			//}

			continue;
		}

		const int boundaryType = calcBoundTypeFromCoordinate( iFace );
		m_boundaryFaces[boundaryType].elemIDs.push_back( elemID );
		m_boundaryFaces[boundaryType].faceIDs.push_back( faceID );
		
	}

	//for( int iElem = 0; iElem < m_numElements; ++iElem ){

	//	for( int i = 0; i < 4; ++i ){
	//		const int elemIDNeighFromZero = m_elements[iElem].neigh[i] - 1;
	//		if(  < 0 ;
	//	}
	//	
	//}

}

//// Change face data from TetGen to Femtic
//void MeshData::changeFaceData(){
//
//	m_faceFemtic = new FaceFemtic[m_numFaces];
//
//	for( int iFace = 0; iFace < m_numFaces; ++iFace ){
//
//		int elemID(0);
//		int faceID(0);
//		calcElementFace( iFace, elemID, faceID );
//
//		m_faceFemtic[iFace].elemID = elemID;
//		m_faceFemtic[iFace].faceID = faceID;
//		m_faceFemtic[iFace].attribute = m_faceTetGen[iFace].attribute;
//
//	}
//
//	if( m_faceTetGen != NULL ){
//		delete [] m_faceTetGen;
//		m_faceTetGen = NULL;
//	}
//
//}

//// Calculate element ID and face ID for Femtic
//void MeshData::calcElementFace( const int iFace, int& elemID, int& faceID ) const{
//
//	if( iFace < 0 || iFace >= m_numFaces ){
//		std::cerr << "iFace ( " << iFace << " ) is wrong !!" << std::endl;
//		exit(1);
//	}
//
//	const int nodeIDs[3] = {
//		m_faceTetGen[iFace].nodeID[0],
//		m_faceTetGen[iFace].nodeID[1],
//		m_faceTetGen[iFace].nodeID[2]
//	};
//
//#ifdef _DEBUG_WRITE
//	std::cout << "iFace : " << iFace << std::endl;
//	std::cout << "nodeIDs : " << nodeIDs[0] << " " <<  nodeIDs[1] << " " <<  nodeIDs[2] << std::endl;
//#endif
//
//	for( std::vector<int>::const_iterator itr = m_node2elem[ nodeIDs[0] - 1 ].begin();
//		 itr != m_node2elem[ nodeIDs[0] - 1 ].end();
//		 ++itr ){
//
//		elemID = *itr;
//
//#ifdef _DEBUG_WRITE
//		std::cout << "elemID : " << elemID << std::endl;
//#endif
//
//		std::vector<int> localNodeID;
//		int icount(0);
//		for( int iFaceNode = 0; iFaceNode < 3; ++iFaceNode ){
//
//#ifdef _DEBUG_WRITE
//			std::cout << "iFaceNode : " << iFaceNode << std::endl;
//			std::cout << "nodeIDs[iFaceNode] : " << nodeIDs[iFaceNode] << std::endl;
//#endif
//			bool found(false);
//			for( int iNode = 0; iNode < 4; ++iNode ){
//
//#ifdef _DEBUG_WRITE
//				//std::cout << "m_elements[elemID-1].attribute : " << m_elements[elemID-1].attribute << std::endl;
//				//std::cout << "m_elements[elemID-1].neigh[iNode] : " << m_elements[elemID-1].neigh[iNode] << std::endl;
//				std::cout << "iNode : " << iNode << std::endl;
//				std::cout << "m_elements[elemID-1].nodeID[iNode] : " << m_elements[elemID-1].nodeID[iNode] << std::endl;
//#endif
//				if( m_elements[elemID-1].nodeID[iNode] == nodeIDs[iFaceNode] ){
//					localNodeID.push_back( iNode );
//					found = true;
//					++icount;
//					break;
//				}
//
//			}
//			if( icount >= 3 ){
//				faceID = calcFaceIDOfFemtic( localNodeID );
//				return;
//			}
//			if( !found ){
//				break;
//			}
//		}
//
//	}
//
//	std::cerr << "Element face having node ";
//	std::cerr << nodeIDs[0] << " " << nodeIDs[1] << " " << nodeIDs[2] << " ";
//	std::cerr << "cannot be found !!" << std::endl; 
//	exit(1);
//
//}

// Calculate element ID and face ID for Femtic
void MeshData::calcElementFace( const int boundaryType, const int iFace, int& elemID, int& faceID ) const{

	if( iFace < 0 || iFace >= m_numFaces ){
		std::cerr << "iFace ( " << iFace << " ) is wrong !!" << std::endl;
		exit(1);
	}

	const int nodeIDs[3] = {
		m_faceTetGen[iFace].nodeID[0],
		m_faceTetGen[iFace].nodeID[1],
		m_faceTetGen[iFace].nodeID[2]
	};

//#ifdef _DEBUG_WRITE
//	std::cout << "iFace : " << iFace << std::endl;
//	std::cout << "nodeIDs : " << nodeIDs[0] << " " <<  nodeIDs[1] << " " <<  nodeIDs[2] << std::endl;
//#endif

	for( std::vector<int>::const_iterator itr = m_node2elem[ nodeIDs[0] - 1 ].begin();
		 itr != m_node2elem[ nodeIDs[0] - 1 ].end();
		 ++itr ){

		elemID = *itr;

//#ifdef _DEBUG_WRITE
//		std::cout << "elemID : " << elemID << std::endl;
//#endif

		if( boundaryType == MeshData::EARTH_SURFACE ){
//#ifdef _DEBUG_WRITE
//			std::cout << "m_elements[elemID-1].attribute : " << m_elements[elemID-1].attribute << std::endl;
//			std::cout << "m_regionAttrAirLayer : " << m_regionAttrAirLayer << std::endl;
//#endif
			//if( m_elements[elemID-1].attribute == m_regionAttrAirLayer ){
			if( m_regionAttrAirOrSea.find( m_elements[elemID-1].attribute ) != m_regionAttrAirOrSea.end() ){// The air or the sea
				continue;
			}
		}

		std::vector<int> localNodeID;
		int icount(0);
		for( int iFaceNode = 0; iFaceNode < 3; ++iFaceNode ){

//#ifdef _DEBUG_WRITE
//			std::cout << "iFaceNode : " << iFaceNode << std::endl;
//			std::cout << "nodeIDs[iFaceNode] : " << nodeIDs[iFaceNode] << std::endl;
//#endif
			bool found(false);
			for( int iNode = 0; iNode < 4; ++iNode ){

//#ifdef _DEBUG_WRITE
//				//std::cout << "m_elements[elemID-1].attribute : " << m_elements[elemID-1].attribute << std::endl;
//				//std::cout << "m_elements[elemID-1].neigh[iNode] : " << m_elements[elemID-1].neigh[iNode] << std::endl;
//				std::cout << "iNode : " << iNode << std::endl;
//				std::cout << "m_elements[elemID-1].nodeID[iNode] : " << m_elements[elemID-1].nodeID[iNode] << std::endl;
//#endif
				if( m_elements[elemID-1].nodeID[iNode] == nodeIDs[iFaceNode] ){
					localNodeID.push_back( iNode );
					found = true;
					++icount;
					break;
				}

			}
			if( icount >= 3 ){
				faceID = calcFaceIDOfFemtic( localNodeID );
				return;
			}
			if( !found ){
				break;
			}
		}

	}

	std::cerr << "Element face having node ";
	std::cerr << nodeIDs[0] << " " << nodeIDs[1] << " " << nodeIDs[2] << " ";
	std::cerr << "cannot be found !!" << std::endl; 
	exit(1);

}

// Calculate element ID and face ID for Femtic
void MeshData::calcElementFace( const int iFace, int& elemID, int& faceID ) const{

	if( iFace < 0 || iFace >= m_numFaces ){
		std::cerr << "iFace ( " << iFace << " ) is wrong !!" << std::endl;
		exit(1);
	}

	const int nodeIDs[3] = {
		m_faceTetGen[iFace].nodeID[0],
		m_faceTetGen[iFace].nodeID[1],
		m_faceTetGen[iFace].nodeID[2]
	};

//#ifdef _DEBUG_WRITE
//	std::cout << "iFace : " << iFace << std::endl;
//	std::cout << "nodeIDs : " << nodeIDs[0] << " " <<  nodeIDs[1] << " " <<  nodeIDs[2] << std::endl;
//#endif

	for( std::vector<int>::const_iterator itr = m_node2elem[ nodeIDs[0] - 1 ].begin();
		 itr != m_node2elem[ nodeIDs[0] - 1 ].end();
		 ++itr ){

		elemID = *itr;

//#ifdef _DEBUG_WRITE
//		std::cout << "elemID : " << elemID << std::endl;
//#endif

//		if( boundaryType == MeshData::EARTH_SURFACE ){
//#ifdef _DEBUG_WRITE
//			std::cout << "m_elements[elemID-1].attribute : " << m_elements[elemID-1].attribute << std::endl;
//			std::cout << "m_regionAttrAirLayer : " << m_regionAttrAirLayer << std::endl;
//#endif
//			if( m_elements[elemID-1].attribute == m_regionAttrAirLayer ){
//				continue;
//			}
//		}

		std::vector<int> localNodeID;
		int icount(0);
		for( int iFaceNode = 0; iFaceNode < 3; ++iFaceNode ){

//#ifdef _DEBUG_WRITE
//			std::cout << "iFaceNode : " << iFaceNode << std::endl;
//			std::cout << "nodeIDs[iFaceNode] : " << nodeIDs[iFaceNode] << std::endl;
//#endif
			bool found(false);
			for( int iNode = 0; iNode < 4; ++iNode ){

//#ifdef _DEBUG_WRITE
//				//std::cout << "m_elements[elemID-1].attribute : " << m_elements[elemID-1].attribute << std::endl;
//				//std::cout << "m_elements[elemID-1].neigh[iNode] : " << m_elements[elemID-1].neigh[iNode] << std::endl;
//				std::cout << "iNode : " << iNode << std::endl;
//				std::cout << "m_elements[elemID-1].nodeID[iNode] : " << m_elements[elemID-1].nodeID[iNode] << std::endl;
//#endif
				if( m_elements[elemID-1].nodeID[iNode] == nodeIDs[iFaceNode] ){
					localNodeID.push_back( iNode );
					found = true;
					++icount;
					break;
				}

			}
			if( icount >= 3 ){
				faceID = calcFaceIDOfFemtic( localNodeID );
				return;
			}
			if( !found ){
				break;
			}
		}

	}

	std::cerr << "Element face having node ";
	std::cerr << nodeIDs[0] << " " << nodeIDs[1] << " " << nodeIDs[2] << " ";
	std::cerr << "cannot be found !!" << std::endl; 
	exit(1);

}

// Calculate local face ID of Femtic
int MeshData::calcFaceIDOfFemtic( std::vector<int>& localNodeIDs ) const{

	if( static_cast<int>( localNodeIDs.size() ) != 3 ){
		std::cerr << "Number of local node IDs must be 3 !!" << std::endl; 
		exit(1);
	}

	//std::sort( localNodeIDs.begin(), localNodeIDs.end() ); 

	//if( localNodeIDs[0] == 1 &&
	//	localNodeIDs[1] == 2 &&
	//	localNodeIDs[2] == 3 ){
	//	return 0;
	//}
	//else
	//if( localNodeIDs[0] == 0 &&
	//	localNodeIDs[1] == 2 &&
	//	localNodeIDs[2] == 3 ){
	//	return 1;
	//}
	//else
	//if( localNodeIDs[0] == 0 &&
	//	localNodeIDs[1] == 1 &&
	//	localNodeIDs[2] == 3 ){
	//	return 2;
	//}
	//else
	//if( localNodeIDs[0] == 0 &&
	//	localNodeIDs[1] == 1 &&
	//	localNodeIDs[2] == 2 ){
	//	return 3;
	//}

	const int sum = localNodeIDs[0] + localNodeIDs[1] + localNodeIDs[2];

	switch(sum){
		case 6:
			return 0;
			break;
		case 5:
			return 1;
			break;
		case 4:
			return 2;
			break;
		case 3:
			return 3;
			break;
		default:
			std::cerr << "Locall node IDs (" << localNodeIDs[0] << " "  << localNodeIDs[1] << " "  << localNodeIDs[2] << " ) are wrong." << std::endl; 
			exit(1);
			return -1;
			break;
	}

	return -1;

}

// Calculate boundary type from attribute
int MeshData::calcBoundTypeFromAttribute( const int iAttr, bool& isBoundary ) const{

	for( int iBoun = 0; iBoun < NUM_BOUNDARY; ++iBoun ){

		for( std::vector<int>::const_iterator itr = m_boundaryFaces[iBoun].attributes.begin(); itr != m_boundaryFaces[iBoun].attributes.end(); ++itr ){
			if( iAttr == *itr ){
				isBoundary = true;
				return iBoun;
			}
		}

	}

	isBoundary = false;

#ifdef _DEBUG_WRITE
	std::cout << "Attribute ( " << iAttr << " ) is not assigned to a boundary !!" << std::endl; 
#endif

	return 0;

}

// Calculate boundary type from coordinate
int MeshData::calcBoundTypeFromCoordinate( const int iFace ) const{

	const double EPS = 1.0e-9;

	if( iFace < 0 || iFace >= m_numFaces ){
		std::cerr << "iFace ( " << iFace << " ) is wrong !!" << std::endl;
		exit(1);
	}

	const int nodeIDFromZero[3] = {
		m_faceTetGen[iFace].nodeID[0] - 1,
		m_faceTetGen[iFace].nodeID[1] - 1,
		m_faceTetGen[iFace].nodeID[2] - 1
	};

	if( fabs( m_nodeCoords[nodeIDFromZero[0]].X - m_nodeCoords[nodeIDFromZero[1]].X ) < EPS &&
		fabs( m_nodeCoords[nodeIDFromZero[0]].X - m_nodeCoords[nodeIDFromZero[2]].X ) < EPS ){
		// On Y-Z plane

		if( m_nodeCoords[nodeIDFromZero[0]].X > 0 ){
			return YZPlus;
		}
		else{
			return YZMinus;
		}
	}
	else
	if( fabs( m_nodeCoords[nodeIDFromZero[0]].Y - m_nodeCoords[nodeIDFromZero[1]].Y ) < EPS &&
		fabs( m_nodeCoords[nodeIDFromZero[0]].Y - m_nodeCoords[nodeIDFromZero[2]].Y ) < EPS ){
		// On Z-X plane

		if( m_nodeCoords[nodeIDFromZero[0]].Y > 0 ){
			return ZXPlus;
		}
		else{
			return ZXMinus;
		}
	}
	else
	if( fabs( m_nodeCoords[nodeIDFromZero[0]].Z - m_nodeCoords[nodeIDFromZero[1]].Z ) < EPS &&
		fabs( m_nodeCoords[nodeIDFromZero[0]].Z - m_nodeCoords[nodeIDFromZero[2]].Z ) < EPS ){
		// On X-Y plane

		if( m_nodeCoords[nodeIDFromZero[0]].Z > 0 ){
			return XYPlus;
		}
		else{
			return XYMinus;
		}
	}

	std::cerr << "Face " << iFace << " is not on the boundary plane !!" << std::endl;
	exit(1);

}

// Write mesh data
void MeshData::writeMeshData() const{

	std::ofstream ofs( "mesh.dat" );

	if( !ofs.is_open() ){
		std::cerr << "Cannot open file mesh.dat !!" << std::endl;
		exit(1);
	}

	std::cout << "Output mesh data to mesh.dat ." << std::endl;

	ofs << "TETRA" << std::endl;

	ofs << std::setw(10) << m_numNodes << std::endl;

	ofs.precision(9);
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		ofs << std::setw(10) << iNode
			<< std::setw(20) << std::scientific << m_nodeCoords[iNode].X
			<< std::setw(20) << std::scientific << m_nodeCoords[iNode].Y
			<< std::setw(20) << std::scientific << m_nodeCoords[iNode].Z
			<< std::endl;
	}

	ofs << std::setw(10) << m_numElements << std::endl;
	for( int iElem = 0; iElem < m_numElements; ++iElem ){
		ofs << std::setw(10) << iElem
			<< std::setw(10) << ( m_elements[iElem].neigh[0] > 0 ? m_elements[iElem].neigh[0] - 1 : -1 )
			<< std::setw(10) << ( m_elements[iElem].neigh[1] > 0 ? m_elements[iElem].neigh[1] - 1 : -1 ) 
			<< std::setw(10) << ( m_elements[iElem].neigh[2] > 0 ? m_elements[iElem].neigh[2] - 1 : -1 )
			<< std::setw(10) << ( m_elements[iElem].neigh[3] > 0 ? m_elements[iElem].neigh[3] - 1 : -1 )
			<< std::setw(10) << m_elements[iElem].nodeID[0] - 1
			<< std::setw(10) << m_elements[iElem].nodeID[1] - 1
			<< std::setw(10) << m_elements[iElem].nodeID[2] - 1
			<< std::setw(10) << m_elements[iElem].nodeID[3] - 1
			<< std::endl;
	}

	for( int iBoun = 0; iBoun < NUM_BOUNDARY; ++iBoun ){
		const int numFace = static_cast<int>( m_boundaryFaces[iBoun].elemIDs.size() );
		ofs << std::setw(10) << numFace << std::endl;
		for( int iFace = 0; iFace < numFace; ++iFace ){
			ofs << std::setw(10) << m_boundaryFaces[iBoun].elemIDs[iFace] - 1
				<< std::setw(10) << m_boundaryFaces[iBoun].faceIDs[iFace]
				<< std::endl;
		}
	}

	ofs.close();

}

// Write mesh data to vtk file
void MeshData::writeMeshDataToVTK( const std::string& rootName ) const{

	std::string fileName = rootName;
	fileName += ".femtic.vtk";

	std::ofstream vtkFile( fileName.c_str() );

	if( !vtkFile.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Output mesh data to " << fileName.c_str() << std::endl;

	vtkFile << "# vtk DataFile Version 2.0" << std::endl;
	vtkFile << "MeshData" << std::endl;
	vtkFile << "ASCII" << std::endl;
	vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtkFile << "POINTS " << m_numNodes << " double" << std::endl;

	vtkFile.precision(9);
	for( int iNode = 0; iNode < m_numNodes; ++iNode ){
		vtkFile << std::setw(20) << std::scientific << m_nodeCoords[iNode].X
			    << std::setw(20) << std::scientific << m_nodeCoords[iNode].Y
				<< std::setw(20) << std::scientific << m_nodeCoords[iNode].Z
				<< std::endl;
	}

	vtkFile << "CELLS " << m_numElements << " " << m_numElements * 5 << std::endl;
	for( int iElem = 0; iElem < m_numElements; ++iElem ){
		vtkFile << std::setw(5)  << 4
				<< std::setw(10) << m_elements[iElem].nodeID[0] - 1
				<< std::setw(10) << m_elements[iElem].nodeID[1] - 1
				<< std::setw(10) << m_elements[iElem].nodeID[2] - 1
				<< std::setw(10) << m_elements[iElem].nodeID[3] - 1 << std::endl;
	}

	vtkFile << "CELL_TYPES " << m_numElements << std::endl;
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		vtkFile << std::setw(5) << 10 << std::endl;
	}

	vtkFile << "CELL_DATA " << m_numElements << std::endl;
	vtkFile << "SCALARS BlockSerial int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		vtkFile << std::setw(10) << m_elem2blk[iElem] << std::endl;
	}

	vtkFile << "SCALARS ElemPerBlock int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		const int iBlk = m_elem2blk[iElem];
		vtkFile << std::setw(10) << static_cast<int>( m_blk2elem[iBlk].size() ) << std::endl;
	}

	vtkFile << "SCALARS Resistivity[Ohm-m] double" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	vtkFile.precision(9);
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		vtkFile << std::setw(20) << m_blk2resistivity[ m_elem2blk[iElem] ].first << std::endl;
	}

	vtkFile << "SCALARS ElemSerial int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		vtkFile << iElem << std::endl;
	}

	vtkFile << "SCALARS Fixed int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		vtkFile << std::setw(20) << static_cast<int>( m_blk2resistivity[ m_elem2blk[iElem] ].second ) << std::endl;
	}

#ifdef _DEBUG_WRITE
	int* boundFlag = new int[m_numElements];
	for( int iBoun = 0; iBoun < NUM_BOUNDARY; ++iBoun ){
		for( int i = 0; i < m_numElements; ++i ){
			boundFlag[i] = 0;
		}
		const int numFace = static_cast<int>( m_boundaryFaces[iBoun].elemIDs.size() );
		for( int iFace = 0; iFace < numFace; ++iFace ){
			boundFlag[ m_boundaryFaces[iBoun].elemIDs[iFace] - 1 ] = 1;
		}
		vtkFile << "SCALARS ";
		switch( iBoun ){
			case MeshData::YZMinus:
				vtkFile << "YZMinus";
				break;
			case MeshData::YZPlus:
				vtkFile << "YZPlus";
				break;
			case MeshData::ZXMinus:
				vtkFile << "ZXMinus";
				break;
			case MeshData::ZXPlus:
				vtkFile << "ZXPlus";
				break;
			case MeshData::XYMinus:
				vtkFile << "XYMinus";
				break;
			case MeshData::XYPlus:
				vtkFile << "XYPlus";
				break;
			case MeshData::EARTH_SURFACE:
				vtkFile << "EARTH_SURFACE";
				break;
			default:
				std::cerr << "Wrong boundary type !!" << std::endl;
				exit(1);
				break;
		}
		vtkFile << " int" <<  std::endl;
		vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
			vtkFile << boundFlag[iElem] << std::endl;
		}
	}

	delete [] boundFlag;
#endif

	vtkFile << "POINT_DATA " << m_numNodes << std::endl;
	vtkFile << "SCALARS NodeSerial int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( int iNode = 0 ; iNode < m_numNodes; ++iNode ){
		vtkFile << iNode << std::endl;
	}

	vtkFile.close();

}


//// Write resistivity data
//void MeshData::writeResisitivityData() const{
//
//	std::ofstream ofs( "resistivity_model_iter0.dat" );
//
//	if( !ofs.is_open() ){
//		std::cerr << "Cannot open file resistivity_model_iter0.dat !!" << std::endl;
//		exit(1);
//	}
//
//	std::cout << "Output mesh data to resistivity_model_iter0.dat ." << std::endl;
//
//	//const int numResisitivityBlk = static_cast<int>( m_regionAttr2Resistivity.size() );
//	ofs << std::setw(10) << m_numElements
//		<< std::setw(10) << m_numResistivity
//		<< std::endl;
//
//	for( int iElem = 0; iElem < m_numElements; ++iElem ){
//
//		ofs << std::setw(10) << iElem
//			<< std::setw(10) << (m_regionAttr2Resistivity.find( m_elements[iElem].attribute )->second).first
//			<< std::endl;
//		
//	}
//
//	ofs.precision(9);
//	for( int iBlk = 0; iBlk < numResisitivityBlk; ++iBlk ){
//
//		double val(0.0);
//		bool found(false);
//		for( std::map< int, std::pair<int,double> >::const_iterator itr = m_regionAttr2Resistivity.begin();
//			 itr != m_regionAttr2Resistivity.end();
//			 ++itr ){
//
//			 if( iBlk == (itr->second).first ){
//				 val = (itr->second).second;
//				 found = true;
//				 break;
//			 }
//
//		}
//		if( !found ){
//			std::cerr << "Resistivity block ID " << iBlk << " cannot be found !!" << std::endl;
//			exit(1);
//		}
//		
//		ofs << std::setw(10) << iBlk
//			<< std::setw(20) << std::scientific << val
//			<< std::setw(10) << 0
//			<< std::endl;
//		
//	}
//
//	ofs.close();
//
//}

// Write resistivity data
void MeshData::writeResisitivityData() const{

	std::ofstream ofs( "resistivity_block_iter0.dat" );

	if( !ofs.is_open() ){
		std::cerr << "Cannot open file resistivity_block_iter0.dat !!" << std::endl;
		exit(1);
	}

	std::cout << "Output mesh data to resistivity_block_iter0.dat ." << std::endl;

	//const int numResisitivityBlk = static_cast<int>( m_regionAttr2Resistivity.size() );
	ofs << std::setw(10) << m_numElements
		<< std::setw(10) << static_cast<int>( m_blk2resistivity.size() )
		<< std::endl;

	for( int iElem = 0; iElem < m_numElements; ++iElem ){

		ofs << std::setw(10) << iElem
			<< std::setw(10) << m_elem2blk[iElem]
			<< std::endl;
		
	}

	ofs.precision(6);
	int icount(0);
	for( std::vector< std::pair<double,bool> >::const_iterator itr = m_blk2resistivity.begin(); itr != m_blk2resistivity.end(); ++itr, ++icount ){
		ofs << std::setw(10) << icount
			<< std::setw(15) << std::scientific << itr->first
#ifdef _OLD
#else
			<< std::setw(15) << 1.0e-20 << std::setw(15) << 1.0e+20 << std::setw(15) << 1.0
#endif
			<< std::setw(10) << ( itr->second ? 1 : 0 )
			<< std::endl;
	}

	ofs.close();

}

// Make resistivity block by partitioning resistivity block by RCB
void MeshData::makeResistivityBlock( const bool divAllElements ){

	const int numResistivity = static_cast<int>( m_resistivityData.size() );
	ElementVecs* elementVecsArray = new ElementVecs[numResistivity];

	for( int iRes = 0; iRes < numResistivity; ++iRes ){
		std::vector<int> dum;
		elementVecsArray[iRes].push_back( std::make_pair( dum, false ) );
	}

	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		std::map< int, int >::iterator itr = m_regionAttr2ResistivityID.find( m_elements[iElem].attribute );
		if (itr == m_regionAttr2ResistivityID.end()){
			std::cerr << "Error : Region attiribute << " << m_elements[iElem].attribute << " of element " << iElem << " cannot be found !!" << std::endl;
			exit(1);
		}
		const int iRes = itr->second;
		elementVecsArray[iRes][0].first.push_back( iElem );
	}

	int numResistivityBlocks = 1;// For the air layer
	for( int iRes = 1; iRes < numResistivity; ++iRes ){// Skip iRes = 0 ( the air layer )

		if( m_resistivityData[iRes].ndiv < 0 ){
			++numResistivityBlocks;
			continue;
		}

		if(divAllElements){
			divideElementVecsAll( elementVecsArray[iRes] );
			numResistivityBlocks += static_cast<int>( elementVecsArray[iRes].size() );
		}else{
			divideElementVecs( elementVecsArray[iRes],  m_resistivityData[iRes].ndiv );
			if( m_numSphere > 0  ){
				divideBlockUntilRequiredLength( elementVecsArray[iRes] );
			}
			numResistivityBlocks += static_cast<int>( elementVecsArray[iRes].size() );
		}

	}

	m_elem2blk = new int[m_numElements];
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		m_elem2blk[iElem] = -1;
	}

	std::vector<int> tmpVec;
	m_blk2elem.resize(numResistivityBlocks, tmpVec);
	//m_blkID2resistivity.resize(m_numResistivityBlocks, std::make_pair(0.0,false));
	m_blk2resistivity.reserve(numResistivityBlocks);

	//std::cout << "<blockID> <#elements> <resistivity> <fix>" << std::endl; 

	int icount(0);
	for( int iRes = 0; iRes < numResistivity; ++iRes ){
		for( ElementVecs::iterator itrVecs = elementVecsArray[iRes].begin(); itrVecs != elementVecsArray[iRes].end(); ++itrVecs ){
			for( std::vector<int>::iterator itr = (itrVecs->first).begin(); itr != (itrVecs->first).end(); ++itr ){
				m_elem2blk[*itr] = icount;
				m_blk2elem[icount].push_back(*itr);
			}
			m_blk2resistivity.push_back( std::make_pair(m_resistivityData[iRes].resistivity,m_resistivityData[iRes].fix) );
			//std::cout << std::setw(10) << static_cast<int>( m_blk2resistivity.size() )
			//		  << std::setw(10) << static_cast<int>( m_blk2elem[icount].size() )
			//		  << std::setw(20) << std::scientific << m_blk2resistivity.back().first
			//		  << std::setw(5) << std::scientific << m_blk2resistivity.back().second
			//		  << std::endl; 
			++icount;
		}

	}

	delete [] elementVecsArray;

}

// Divide elements into sub-elements
void MeshData::divideElementVecs( ElementVecs& elementVecs, const int ndiv ) const{

	for( int idiv = 0; idiv < ndiv; ++idiv ){

		const int num = static_cast<int>( elementVecs.size() );
		bool everyElementHasDifferentResistivity(true);
		for( int i = 0; i < num; ++i ){
			if( static_cast<int>( elementVecs[i].first.size() ) < 2 ){
				continue;
			}
			elementVecs.push_back( std::make_pair( divide( elementVecs[i].first ), false ) );
			everyElementHasDifferentResistivity = false;
		}

		if(everyElementHasDifferentResistivity){
			return;
		}

	}
	
}

// Assign each element to different parameter cell
void MeshData::divideElementVecsAll( ElementVecs& elementVecs ) const{

	ElementVecs elementVecsNew;
	const int num = static_cast<int>( elementVecs.size() );
	for( int i = 0; i < num; ++i ){
		for( std::vector<int>::iterator itr = elementVecs[i].first.begin(); itr != elementVecs[i].first.end(); ++itr ){
			std::vector<int> elemVec;
			elemVec.push_back(*itr);
			elementVecsNew.push_back( std::make_pair( elemVec, false ) );
		}
	}
	elementVecs.swap(elementVecsNew);

}

// Divide resistivity blocks to the required length
void MeshData::divideBlockUntilRequiredLength( ElementVecs& elementVecs ) const{

	int icount = 0;
	for( ; icount < m_numElements; ++icount ){

		bool satisfyAll = true;

		int numElemVec = static_cast<int>( elementVecs.size() );

		for( int ivec = 0; ivec < numElemVec; ++ivec ){

			if( elementVecs[ivec].second ){
				// Length of resistivity blocks is shorter than the required length

				//std::cout << "Length of resistivity blocks is shorter than the required length" << std::endl;

				continue;
			}
			
			if( static_cast<int>( elementVecs[ivec].first.size() ) < 2 ){
				elementVecs[ivec].second = true;

				//std::cout << "elementVecs[ivec].first.size() : " << elementVecs[ivec].first.size() << std::endl;

				continue;
			}

			std::vector< std::pair<double,int> > stack[3];

			for( std::vector<int>::iterator itr = (elementVecs[ivec].first).begin(); itr != (elementVecs[ivec].first).end(); ++itr ){
				stack[0].push_back( std::make_pair( calcGravityCenterX(*itr), *itr ) );  
				stack[1].push_back( std::make_pair( calcGravityCenterY(*itr), *itr ) );  
				stack[2].push_back( std::make_pair( calcGravityCenterZ(*itr), *itr ) );  
			}

			std::vector< std::pair<double,int> > stackDis;
			for( int i = 0; i < 3; ++i ){
				sort( stack[i].begin(), stack[i].end() );
				stackDis.push_back( std::make_pair( fabs( stack[i].back().first - stack[i].front().first ), i ) );
			}

			sort( stackDis.begin(), stackDis.end() );

			const int idirMax = stackDis.back().second;
			const double disMax = stackDis.back().first;

			if( disMax < calcShorterBlockLength( stack[idirMax].front().second, stack[idirMax].back().second ) ){
				elementVecs[ivec].second = true;

				//std::cout << "disMax required : " << disMax << " " << calcShorterBlockLength( stack[idirMax].front().second, stack[idirMax].back().second ) << std::endl;

				continue;
			}

			satisfyAll = false;

			const int iEnd = static_cast<int>( stack[idirMax].size() );
			const int iMid = iEnd / 2 + iEnd % 2;

			std::vector<int> newElements1;
			for( int i = 0; i < iMid; ++i ){
				newElements1.push_back( stack[idirMax][i].second );
			}
			(elementVecs[ivec].first).swap( newElements1 );

			std::vector<int> newElements2;
			for( int i = iMid; i < iEnd; ++i ){
				newElements2.push_back( stack[idirMax][i].second );
			}
			
			elementVecs.push_back( std::make_pair( newElements2, false ) );
			
		}

		if( satisfyAll ){
			break;
		}

	}

	if( icount >= m_numElements ){
		std::cerr << "Loop count reach maximum value : " << icount << std::endl;
		exit(1);
	}

}

// Divide specified elements for the longest direction
std::vector<int> MeshData::divide( std::vector<int>& elements ) const{

	if( static_cast<int>( elements.size() ) < 2 ){
		std::cerr << "Error : Number of specified elements is less than 2" << std::endl;
		exit(1);
	}

	std::vector< std::pair<double,int> > stack[3];

	for( std::vector<int>::iterator itr = elements.begin(); itr != elements.end(); ++itr ){
		stack[0].push_back( std::make_pair( calcGravityCenterX(*itr), *itr ) );  
		stack[1].push_back( std::make_pair( calcGravityCenterY(*itr), *itr ) );  
		stack[2].push_back( std::make_pair( calcGravityCenterZ(*itr), *itr ) );  
	}

	std::vector< std::pair<double,int> > stackDis;
	for( int i = 0; i < 3; ++i ){
		sort( stack[i].begin(), stack[i].end() );
		stackDis.push_back( std::make_pair( fabs( stack[i].back().first - stack[i].front().first ), i ) );
	}

	sort( stackDis.begin(), stackDis.end() );

	const int idirMax = stackDis.back().second;

	const int iEnd = static_cast<int>( stack[idirMax].size() );
	const int iMid = iEnd / 2 + iEnd % 2;

	std::vector<int> newElements1;
	for( int i = 0; i < iMid; ++i ){
		newElements1.push_back( stack[idirMax][i].second );
	}
	elements.swap( newElements1 );

	std::vector<int> newElements2;
	for( int i = iMid; i < iEnd; ++i ){
		newElements2.push_back( stack[idirMax][i].second );
	}
	return newElements2;

}

// Calculate shorter block length of elements
double MeshData::calcShorterBlockLength( const int elem1, const int elem2 ) const{

	const XYZ coord1 = {
		calcGravityCenterX(elem1),
		calcGravityCenterY(elem1),
		calcGravityCenterZ(elem1)
	};

	const XYZ coord2 = {
		calcGravityCenterX(elem2),
		calcGravityCenterY(elem2),
		calcGravityCenterZ(elem2)
	};

	return std::min( calcBlockLength( coord1 ), calcBlockLength( coord2 ) );

}

// Calculate maximum block length in the transition region
double MeshData::calcBlockLength( const XYZ& coord ) const{

	double length(1.0e20);
	if( locateInSphere( coord, 0 ) ){
		length = m_maxBlockLength[0];
	}
	else if( !locateInSphere( coord, m_numSphere-1 ) ){
		length = m_maxBlockLength[m_numSphere-1];
	}
	else{
		for( int iSphere = 1; iSphere < m_numSphere; ++iSphere ){
			if( locateInSphere( coord, iSphere ) ){
				length =calcBlockLengthTransitionRegion( coord, iSphere );
				break;
			}
		}
	}

	return std::min( calcBlockLengthNearSite(coord), length );
}

// Decide whether the spacified point locate within the specified hemi-sphere
bool MeshData::locateInSphere( const XYZ& coord, const int iSphere ) const{

	if( iSphere < 0 || iSphere >= m_numSphere ){
		std::cerr << "Wrong shepre ID:  " << iSphere << std::endl;
		exit(1);
	}

	//const double depth = m_radius[iSphere] * ( 1.0 - m_oblateness[iSphere] );

	//double val = pow( ( coord.X - m_centerCoordSpheres.X ) / m_radius[iSphere], 2 )
	//			+ pow( ( coord.Y - m_centerCoordSpheres.Y ) / m_radius[iSphere], 2 )
	//			+ pow( ( coord.Z - m_centerCoordSpheres.Z ) / depth, 2 );
	
	const double vecXOrg = coord.X - m_centerCoordSpheres.X;
	const double vecYOrg = coord.Y - m_centerCoordSpheres.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoordSpheres.Z;

	const double longAxisLength = m_radius[iSphere];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iSphere] );
	const double depth = longAxisLength * ( 1.0 - m_oblateness[iSphere] );

	double val = pow( vecX / longAxisLength, 2 )
			   + pow( vecY / shortAxisLength, 2 )
			   + pow( vecZ / depth, 2 );

	if( val <= 1.0 ){
		return true;
	}

	return false;

}

// Calculate maximum block length in the transition region
double MeshData::calcBlockLengthTransitionRegion( const XYZ& coord, const int iSphere ) const{
	
	if( iSphere < 1 || iSphere >= m_numSphere ){
		std::cerr << "Wrong shepre ID:  " << iSphere << std::endl;
		exit(1);
	}

	//const double radiusInner = m_radius[iSphere-1];
	//const double radiusOuter = m_radius[iSphere];
	//const double edgeLengthInner = m_maxBlockLength[iSphere-1];
	//const double edgeLengthOuter = m_maxBlockLength[iSphere];
	//const double oblatenessInner = m_oblateness[iSphere-1];
	//const double oblatenessOuter = m_oblateness[iSphere];

	//const double radiusHorizontal = hypot( coord.X - m_centerCoordSpheres.X, coord.Y - m_centerCoordSpheres.Y );
	//const double radius = hypot( radiusHorizontal, coord.Z - m_centerCoordSpheres.Z );
	//const double slope = fabs( coord.Z - m_centerCoordSpheres.Z ) / radiusHorizontal;

	//const double radius0 = radiusInner * ( 1.0 - oblatenessInner ) * sqrt( ( slope*slope + 1.0 ) / ( slope*slope + (1.0-oblatenessInner)*(1.0-oblatenessInner) ) );
	//const double radius1 = radiusOuter * ( 1.0 - oblatenessOuter ) * sqrt( ( slope*slope + 1.0 ) / ( slope*slope + (1.0-oblatenessOuter)*(1.0-oblatenessOuter) ) );

	//if( radius1 < radius0 ){
	//	std::cerr << "Specified coordinate ( " << coord.X << ", " << coord.Y << ", " << coord.Z << " ) does't locate in the trasition region" << std::endl;
	//	exit(1);
	//}

	//if( radius < radius0 ){
	//	std::cerr << "Specified coordinate ( " << coord.X << ", " << coord.Y << ", " << coord.Z << " ) does't locate in the trasition region" << std::endl;
	//	std::cerr << "iSphere : " << iSphere << std::endl;
	//	std::cerr << "radiusHorizontal : " << radiusHorizontal << std::endl;
	//	std::cerr << "radius : " << radius << std::endl;
	//	std::cerr << "radius0 : " << radius0 << std::endl;
	//	std::cerr << "radius1 : " << radius1 << std::endl;
	//	exit(1);
	//}

	//return ( radius - radius0 ) / ( radius1 - radius0 ) * ( edgeLengthOuter - edgeLengthInner ) + edgeLengthInner;

	const double vecXOrg = coord.X - m_centerCoordSpheres.X;
	const double vecYOrg = coord.Y - m_centerCoordSpheres.Y;
	// Coordinate transform
	const double vecX = vecXOrg * cos( - m_rotationAngle ) - vecYOrg * sin( - m_rotationAngle );
	const double vecY = vecXOrg * sin( - m_rotationAngle ) + vecYOrg * cos( - m_rotationAngle );
	const double vecZ = coord.Z - m_centerCoordSpheres.Z;

	const double angleHorizontal = atan2( vecY, vecX );
	const double lengthHorizontal = hypot( vecY, vecX );
	const double angleVertical = atan2( vecZ, lengthHorizontal );
	const double length = hypot( lengthHorizontal, vecZ );

	const double length0 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iSphere-1 );
	const double length1 = calculateLengthOnEllipsoid( angleHorizontal, angleVertical, iSphere );

	if( length < length0 ){
		std::cerr << "Specified coordinate ( " << coord.X << ", " << coord.Y << ", " << coord.Z << " ) does't locate in the trasition region" << std::endl;
		std::cerr << "iSphere : " << iSphere << std::endl;

		std::cerr << "vecX : " << vecX << std::endl;
		std::cerr << "vecY : " << vecY << std::endl;
		std::cerr << "vecZ : " << vecZ << std::endl;

		std::cerr << "lengthHorizontal      : " << lengthHorizontal << std::endl;
		std::cerr << "angleH : " << angleHorizontal << std::endl;
		std::cerr << "angleV : " << angleVertical << std::endl;

		std::cerr << "length  : " << length << std::endl;
		std::cerr << "length0 : " << length0 << std::endl;
		std::cerr << "length1 : " << length1 << std::endl;
		exit(1);
	}

	return ( length - length0 ) / ( length1 - length0 ) * ( m_maxBlockLength[iSphere] - m_maxBlockLength[iSphere-1] ) + m_maxBlockLength[iSphere-1];


}

// Calculate length on ellipsoid
double MeshData::calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iSphere ) const{

	if( iSphere < 0 || iSphere >= m_numSphere ){
		std::cerr << "Wrong shepre ID:  " << iSphere << std::endl;
		exit(1);
	}

	if( angleH < -PI || angleH > PI ){
		std::cerr << "Horizontal angle is improper : " << angleH << std::endl;
		exit(1);
	}

	if( angleV < -PI || angleV > PI ){
		std::cerr << "Vertical angle is improper : " << angleV << std::endl;
		exit(1);
	}

	const double longAxisLength = m_radius[iSphere];
	const double shortAxisLength = longAxisLength * ( 1.0 - m_oblatenessHorizontal[iSphere] );
	const double verticalLength = longAxisLength * ( 1.0 - m_oblateness[iSphere] );
	
	const double eps = 1.0e-9;
	double lengthH(-1.0);
	if( fabs(angleH - PI*0.5) < eps || fabs(angleH + PI*0.5) < eps ){
		lengthH = shortAxisLength;
	}
	else{
		const double constValH = longAxisLength * shortAxisLength / hypot( shortAxisLength, longAxisLength * tan(angleH) );
		lengthH = hypot( constValH, constValH * tan(angleH) );	
	}

	if( fabs(angleV - PI*0.5) < eps || fabs(angleV + PI*0.5) < eps ){
		return verticalLength;
	}
	else{
		const double constValV = lengthH * verticalLength / hypot( verticalLength, lengthH * tan(angleV) );
		return hypot( constValV, constValV * tan(angleV) );
	}

}

// Calculate maximum block length near site
double MeshData::calcBlockLengthNearSite( const XYZ& coord ) const{

	double length(1.0e20);

	for( int iObs = 0; iObs < m_numObsSite; ++iObs ){
		const double radius = sqrt( pow( coord.X - m_obsSite[iObs].X, 2 ) + pow( coord.Y - m_obsSite[iObs].Y, 2 ) + pow( coord.Z - m_obsSite[iObs].Z, 2 ) );
		
		if( radius <= m_obsSite[iObs].radius[0] ){
			length = std::min(length, m_obsSite[iObs].length[0] );
			continue;
		}

		for( int i = m_obsSite[iObs].numCircle - 1; i >= 1; --i ){
			if( radius <= m_obsSite[iObs].radius[i] ){
				const double lengthTemp = m_obsSite[iObs].length[i-1] + ( radius - m_obsSite[iObs].radius[i-1] ) / (  m_obsSite[iObs].radius[i] -  m_obsSite[iObs].radius[i-1] ) * ( m_obsSite[iObs].length[i] - m_obsSite[iObs].length[i-1] );
				length = std::min(length, lengthTemp );
			}else{
				break;
			}
		}

	}

	return length;

}

// Calculate gravity center of element
MeshData::XYZ MeshData::calcGravityCenter( const int iElem ) const{

	if( iElem < 0 || iElem >= m_numElements ){
		std::cerr << "iElem " << iElem << " is out of range !!" << std::endl;
		exit(1);
	}

	//MeshData::XYZ coord = { 0.0, 0.0, 0.0 };

	//for( int i = 0; i < 4; ++i ){
	//	const int iNode = m_elements[iElem].nodeID[i] - 1;
	//	coord.X += m_nodeCoords[ iNode ].X;
	//	coord.Y += m_nodeCoords[ iNode ].Y;
	//	coord.Z += m_nodeCoords[ iNode ].Z;
	//}

	//coord.X *= 0.25;
	//coord.Y *= 0.25;
	//coord.Z *= 0.25;

	MeshData::XYZ coord = {
		calcGravityCenterX(iElem),	
		calcGravityCenterY(iElem),
		calcGravityCenterZ(iElem)
	};

	return coord;

}

// Calculate X coordinate of gravity center of element
double MeshData::calcGravityCenterX( const int iElem ) const{

	if( iElem < 0 || iElem >= m_numElements ){
		std::cerr << "iElem " << iElem << " is out of range !!" << std::endl;
		exit(1);
	}

	double val(0.0);
	for( int i = 0; i < 4; ++i ){
		val += m_nodeCoords[ m_elements[iElem].nodeID[i] - 1 ].X;
	}

	return val * 0.25;

}

// Calculate Y coordinate of gravity center of element
double MeshData::calcGravityCenterY( const int iElem ) const{

	if( iElem < 0 || iElem >= m_numElements ){
		std::cerr << "iElem " << iElem << " is out of range !!" << std::endl;
		exit(1);
	}

	double val(0.0);
	for( int i = 0; i < 4; ++i ){
		val += m_nodeCoords[ m_elements[iElem].nodeID[i] - 1 ].Y;
	}

	return val * 0.25;

}

// Calculate Z coordinate of gravity center of element
double MeshData::calcGravityCenterZ( const int iElem ) const{

	if( iElem < 0 || iElem >= m_numElements ){
		std::cerr << "iElem " << iElem << " is out of range !!" << std::endl;
		exit(1);
	}

	double val(0.0);
	for( int i = 0; i < 4; ++i ){
		val += m_nodeCoords[ m_elements[iElem].nodeID[i] - 1 ].Z;
	}

	return val * 0.25;

}

#ifdef _ADDITIONAL_FIXED_AREA
// Add if the inputted faces of the earth's surface below which resistivity values are fixed
void MeshData::addIfFaceIndexAboveAdditionalFixedRegion( const int iFace ){

	if( iFace < 0 || iFace >= m_numFaces ){
		std::cerr << "iFace ( " << iFace << " ) is wrong !!" << std::endl;
		exit(1);
	}

	const int nodeIndex[3] = {
		m_faceTetGen[iFace].nodeID[0] - 1,
		m_faceTetGen[iFace].nodeID[1] - 1,
		m_faceTetGen[iFace].nodeID[2] - 1
	};

	const XYZ coordCenter = {
		( m_nodeCoords[nodeIndex[0]].X + m_nodeCoords[nodeIndex[1]].X + m_nodeCoords[nodeIndex[2]].X ) / 3.0,
		( m_nodeCoords[nodeIndex[0]].Y + m_nodeCoords[nodeIndex[1]].Y + m_nodeCoords[nodeIndex[2]].Y ) / 3.0,
		( m_nodeCoords[nodeIndex[0]].Z + m_nodeCoords[nodeIndex[1]].Z + m_nodeCoords[nodeIndex[2]].Z ) / 3.0
	};
	
	if( checkIfCoordInAdditionalFixedRegion( coordCenter ) ){
		m_faceIndexAboveAdditionalFixedRegion.push_back(iFace);
	}

}

// Check if the inputted coordinate locates in the additional areal below which resistivity values are fixed
bool MeshData::checkIfCoordInAdditionalFixedRegion( const XYZ& coord ) const{

	const XYZ coordShifted = { 
		coord.X - m_paramForAdditionalFixedRegion.centerCoord.X,
		coord.Y - m_paramForAdditionalFixedRegion.centerCoord.Y,
		coord.Z
	}; 

	const XYZ coordRotated = {
		coordShifted.X * cos( - m_paramForAdditionalFixedRegion.rotationAngle ) - coordShifted.Y * sin( - m_paramForAdditionalFixedRegion.rotationAngle ),
		coordShifted.X * sin( - m_paramForAdditionalFixedRegion.rotationAngle ) + coordShifted.Y * cos( - m_paramForAdditionalFixedRegion.rotationAngle ),
		coordShifted.Z
	};
	
	if( fabs(coordRotated.X) <= m_paramForAdditionalFixedRegion.xLength &&
		fabs(coordRotated.Y) <= m_paramForAdditionalFixedRegion.yLength &&
		coordRotated.Z >= m_paramForAdditionalFixedRegion.minDepthOfEarthSurface &&
		coordRotated.Z <= m_paramForAdditionalFixedRegion.maxDepthOfEarthSurface ){
		return true;
	}

	return false;

}

// Check if the inputted coordinate locates under the faces below which resistivity values are fixed
bool MeshData::checkIfCoordUnderTheFacesConstrutingAdditionalFixedRegion( const XY& coord ) const{

	const double EPS = 1.0e-6;
	for( std::vector<int>::const_iterator itrFace = m_faceIndexAboveAdditionalFixedRegion.begin(); itrFace != m_faceIndexAboveAdditionalFixedRegion.end(); ++itrFace ){
		const int faceIndex = *itrFace;
		const int nodeIndex[3] = {
			m_faceTetGen[faceIndex].nodeID[0] - 1,
			m_faceTetGen[faceIndex].nodeID[1] - 1,
			m_faceTetGen[faceIndex].nodeID[2] - 1
		};

		const XY coordTri[3] = {
			{ m_nodeCoords[nodeIndex[0]].X, m_nodeCoords[nodeIndex[0]].Y },
			{ m_nodeCoords[nodeIndex[1]].X, m_nodeCoords[nodeIndex[1]].Y },
			{ m_nodeCoords[nodeIndex[2]].X, m_nodeCoords[nodeIndex[2]].Y }
		};

		//const double det0 = calcOuterProduct( coordTri[0], coordTri[1], coord  );
		//const double det1 = calcOuterProduct( coordTri[1], coordTri[2], coord  );
		//const double det2 = calcOuterProduct( coordTri[2], coordTri[0], coord  );

		//if( det0 > -EPS && det1 > -EPS && det2 > -EPS ){
		//	// Locate in triangle
		//	return true;
		//}

		const double det  = fabs( calcOuterProduct( coordTri[0], coordTri[1], coordTri[2] ) );
		const double det0 = fabs( calcOuterProduct( coordTri[0], coordTri[1], coord ) );
		const double det1 = fabs( calcOuterProduct( coordTri[1], coordTri[2], coord ) );
		const double det2 = fabs( calcOuterProduct( coordTri[2], coordTri[0], coord ) );

		if( fabs( det0 + det1 + det2 - det ) / det < EPS ){
			// Locate in triangle
			return true;
		}
	}

	return false;

}

// Fix resistivity in the additional area in which resistivity values are fixed
void MeshData::fixResistivityInAdditionalFixedRegion(){

	const int attrAdditionalFixedRegion = m_paramForAdditionalFixedRegion.attribute;

	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){		
		if( m_regionAttrAirOrSea.find( m_elements[iElem].attribute ) != m_regionAttrAirOrSea.end() ){
			// Exclude sea or air
			continue;
		}
		const XY coordXY = { calcGravityCenterX(iElem), calcGravityCenterY(iElem) };
		const double coordZ = calcGravityCenterZ(iElem);
		if( checkIfCoordUnderTheFacesConstrutingAdditionalFixedRegion( coordXY ) ){
			// Locate in triangle
			if( fabs( coordZ - m_paramForAdditionalFixedRegion.centerCoord.Z ) < m_paramForAdditionalFixedRegion.zLength ){
				// Locate in fixed region
				m_elements[iElem].attribute = attrAdditionalFixedRegion;
			}
		}
	}

}

// Calculate gravity centers of resistivity blocks
MeshData::XYZ MeshData::calcGravityCenterOfResistivityBlock( const int blockIndex ) const{

	XYZ coordGravityCenter = { 0.0, 0.0, 0.0 };
	double sumVolume(0.0);

	for( std::vector<int>::const_iterator itr = m_blk2elem[blockIndex].begin(); itr != m_blk2elem[blockIndex].end(); ++itr ){
		const int elemIndex = *itr;
		const double volume = calculateVolumeOfElement(elemIndex);
		MeshData::XYZ center = {
			calcGravityCenterX(elemIndex),	
			calcGravityCenterY(elemIndex),
			calcGravityCenterZ(elemIndex)
		};
		coordGravityCenter.X += center.X * volume;
		coordGravityCenter.Y += center.Y * volume;
		coordGravityCenter.Z += center.Z * volume;
		sumVolume += volume;
	}

	const double div = 1.0 / sumVolume;
	coordGravityCenter.X *= div;
	coordGravityCenter.Y *= div;
	coordGravityCenter.Z *= div;

	return coordGravityCenter;

}


// Calculate outer product
double MeshData::calcOuterProduct( const MeshData::XY& startCoord, const MeshData::XY& endCoord, const MeshData::XY& coord ) const{

	return ( endCoord.X - startCoord.X )*( coord.Y - startCoord.Y ) - ( endCoord.Y - startCoord.Y )*( coord.X - startCoord.X );

}

// Write faces above the additional area in which resistivity values are fixed
void MeshData::writeFacesAboveAdditionalFixedRegionToVTK( const std::string& rootName ) const{

	std::map<int,int> nodeIDGlobal2Local;
	std::map<int,int> nodeIDLocal2Global;
	for( std::vector<int>::const_iterator itrFace = m_faceIndexAboveAdditionalFixedRegion.begin(); itrFace != m_faceIndexAboveAdditionalFixedRegion.end(); ++itrFace ){
		const int faceIndex = *itrFace;
		const int nodeID[3] = {
			m_faceTetGen[faceIndex].nodeID[0],
			m_faceTetGen[faceIndex].nodeID[1],
			m_faceTetGen[faceIndex].nodeID[2]
		};
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDGlobal = m_faceTetGen[faceIndex].nodeID[iNode];
			if( nodeIDGlobal2Local.find(nodeIDGlobal) == nodeIDGlobal2Local.end() ){
				// Not found
				const int nodeIDLocal = static_cast<int>( nodeIDGlobal2Local.size() ) + 1;
				nodeIDGlobal2Local.insert( std::make_pair( nodeIDGlobal, nodeIDLocal ) );
				nodeIDLocal2Global.insert( std::make_pair( nodeIDLocal, nodeIDGlobal ) );
			}
		}
	}

	std::string fileName = rootName;
	fileName += ".faceAboveFixedRegion.vtk";

	std::ofstream vtkFile( fileName.c_str() );

	if( !vtkFile.is_open() ){
		std::cerr << "Cannot open file " << fileName.c_str() << std::endl;
		exit(1);
	}

	std::cout << "Output faces above the additional area to " << fileName.c_str() << std::endl;

	vtkFile << "# vtk DataFile Version 2.0" << std::endl;
	vtkFile << "FaceData" << std::endl;
	vtkFile << "ASCII" << std::endl;
	vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtkFile.precision(9);

	const int numNode = static_cast<int>( nodeIDGlobal2Local.size() );
	vtkFile << "POINTS " << numNode << " double" << std::endl;
	for( std::map<int,int>::const_iterator itrNode = nodeIDLocal2Global.begin(); itrNode != nodeIDLocal2Global.end(); ++itrNode ){
		const int nodeIDLocal = itrNode->first;
		std::map<int,int>::const_iterator itr = nodeIDLocal2Global.find(nodeIDLocal);
		if( itr == nodeIDLocal2Global.end() ){
			// Not found
			std::cerr << "Node ID " << nodeIDLocal << " is NOT found in nodeIDLocal2Global !!" << std::endl;
			exit(1);
		}
		const int nodeIndexGlobal = itr->second - 1;
		vtkFile << std::setw(20) << std::scientific << m_nodeCoords[nodeIndexGlobal].X
			    << std::setw(20) << std::scientific << m_nodeCoords[nodeIndexGlobal].Y
				<< std::setw(20) << std::scientific << m_nodeCoords[nodeIndexGlobal].Z
				<< std::endl;
	}

	const int numTriangles = static_cast<int>( m_faceIndexAboveAdditionalFixedRegion.size() );
	vtkFile << "CELLS " << numTriangles << " " << numTriangles * 4 << std::endl;
	for( std::vector<int>::const_iterator itrFace = m_faceIndexAboveAdditionalFixedRegion.begin(); itrFace != m_faceIndexAboveAdditionalFixedRegion.end(); ++itrFace ){
		vtkFile << std::setw(5)  << 3;
		const int faceIndex = *itrFace;
		for( int iNode = 0; iNode < 3; ++iNode ){
			const int nodeIDGlobal = m_faceTetGen[faceIndex].nodeID[iNode];
			std::map<int,int>::const_iterator itr = nodeIDGlobal2Local.find(nodeIDGlobal);
			if( itr == nodeIDGlobal2Local.end() ){
				// Not found
				std::cerr << "Node ID " << nodeIDGlobal << " is NOT found in nodeIDGlobal2Local !!" << std::endl;
				exit(1);
			}
			vtkFile << std::setw(10) << itr->second - 1;
		}
		vtkFile << std::endl;
	}

	vtkFile << "CELL_TYPES " << m_numElements << std::endl;
	for( int iElem = 0 ; iElem < m_numElements; ++iElem ){
		vtkFile << std::setw(5) << 5 << std::endl;
	}

	vtkFile << "CELL_DATA " << m_numElements << std::endl;
	vtkFile << "SCALARS FaceIndex int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( std::vector<int>::const_iterator itrFace = m_faceIndexAboveAdditionalFixedRegion.begin(); itrFace != m_faceIndexAboveAdditionalFixedRegion.end(); ++itrFace ){
		vtkFile << *itrFace << std::endl;
	}

	vtkFile << "POINT_DATA " << m_numNodes << std::endl;
	vtkFile << "SCALARS NodeID int" <<  std::endl;
	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
	for( std::map<int,int>::const_iterator itr = nodeIDLocal2Global.begin(); itr != nodeIDLocal2Global.end(); ++itr ){
		vtkFile << itr->second << std::endl;
	}

	vtkFile.close();

}

#endif
