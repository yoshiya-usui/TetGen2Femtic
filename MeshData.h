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
#ifndef DBLDEF_TETGEN_MESH_DATA
#define DBLDEF_TETGEN_MESH_DATA

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include "TopographicData.h"
// Class of mesh data of tetgen
class MeshData{

public:
	// Return the the instance of the class
    static MeshData* getInstance();

	// Make mesh data for femtic
	void makeMeshDataForFemtic( const std::string& rootName, const bool includeSediments,
		const bool rotation, const int numThread, const bool divAllElements );

private:

	const static double convkm2m;

	const static double PI;

	const static double DEG2RAD;

	struct XY{
		double X;
		double Y;
	};

	struct XYZ{
		double X;
		double Y;
		double Z;
	};

	struct Element{
		int nodeID[4];
		int attribute;
		int neigh[4];
	};

	struct FaceTetGen{
		int nodeID[3];
		int attribute;
	};
	
	//struct FaceFemtic{
	//	int elemID;
	//	int faceID;
	//	int attribute;
	//};

	struct BoundaryFace{
		std::vector<int> elemIDs;
		std::vector<int> faceIDs;
		std::vector<int> attributes;
	};

	struct ResistivityData{
		double resistivity;
		int ndiv;
		int attr;
		bool fix;
	};

	enum BoundaryType{
		YZMinus = 0,
		YZPlus,
		ZXMinus,
		ZXPlus,
		XYMinus,
		XYPlus,
		EARTH_SURFACE,
		NUM_BOUNDARY
	};

	struct ObservedSite{
		double X;
		double Y;
		double Z;
		int numCircle;
		std::vector<double> radius;
		std::vector<double> length;
	};

	enum TypeID{
		EARTH = 0,
		AIR = 1
	};

	struct XYMinMax{
		double xMin;
		double xMax;
		double yMin;
		double yMax;
	};

	typedef std::vector< std::pair< std::vector<int>, bool> > ElementVecs;

	struct SedimentParameter{
		XY center;// Coordinate of the center coordinate of the area including predifined sediments
		double radius;// Radius of the area including predifined sediments
		double oblateness;// Oblateness of the area including predifined sediments
		double rotation;// Rotation angle of the area including predifined sediments
		double depth;// Depth of sedimentary layer
		int attr;// Attribute of sediment
		double resistivity;// Resistivity of sediment
	};

#ifdef _ADDITIONAL_FIXED_AREA
	struct ParametersForAdditionalFixedRegion{
		double minDepthOfEarthSurface;
		double maxDepthOfEarthSurface;
		XYZ centerCoord;
		double xLength;
		double yLength;
		double zLength;
		double rotationAngle;
		double resistivity;
		int attribute;
	};
#endif

	bool m_doesUseBoundaryMaker;

	// Total number of nodes
	int m_numNodes;

	// Coordinates of nodes
	XYZ* m_nodeCoords;

	// Total number of elements
	int m_numElements;

	// Element data
	Element* m_elements;

	// Array convert node IDs to element IDs
	std::vector<int>* m_node2elem;

	// Total number of faces
	int m_numFaces;

	// Face data of TetGen
	FaceTetGen* m_faceTetGen;

	//// Face data of Femtic
	//FaceFemtic* m_faceFemtic;

	// Region attribute of the air layer
	int m_regionAttrAirLayer;

	// Region attribute of the air or the sea
	std::set<int> m_regionAttrAirOrSea;

	//// Array convert boundary ID to attributes of face data
	//std::vector<int> m_bound2attribute[NUM_BOUNDARY];

	// Boundary faces of Femtic
	BoundaryFace m_boundaryFaces[NUM_BOUNDARY];
	
	//// Map convert region attribute to resistivity data
	//std::map< int, std::pair<int,double> > m_regionAttr2Resistivity;
	// Map convert region attribute to resistivity ID
	std::map< int, int > m_regionAttr2ResistivityID;

	// Resistivity data
	std::vector<ResistivityData> m_resistivityData;

	// Array convert element ID to resistivity block ID
	int* m_elem2blk;

	// Array convert resistivity block serial to resistivity value
	std::vector< std::pair<double,bool> > m_blk2resistivity;

	// Array convert resistivity block serial to element serials
	std::vector< std::vector<int> > m_blk2elem;

	// Center coordinate of spheres
	XYZ m_centerCoordSpheres;

	// Rotation angle
	double m_rotationAngle;

	// Total number of spheres
	int m_numSphere;

	// Radiuses of horizontal circules
	double* m_radius;

	// Maximum length of resisitivity block within each spheres
	double* m_maxBlockLength;

	// Oblateness of hemi-spheres on horizontal plane
	double* m_oblatenessHorizontal;

	// Oblateness of hemi-spheres
	double* m_oblateness;

	// Total number of observation sites
	int m_numObsSite;

	// Pointer to the object of the observation sites
	ObservedSite* m_obsSite;

	// Get total number of threads
	int m_numThreads;

#ifdef _ADDITIONAL_FIXED_AREA
	// Index of the faces of the earth's surface below which resistivity values are fixed
	std::vector<int> m_faceIndexAboveAdditionalFixedRegion;

	// Flag indicating whether additonal fixed region was took into account or not
	bool m_useAdditionalFixedRegion;

	// Parameters for additional fixed region
	ParametersForAdditionalFixedRegion m_paramForAdditionalFixedRegion;
#endif

	// Constructer
	explicit MeshData();

	// Destructer
	~MeshData();

	// Copy constructer
	explicit MeshData(const MeshData& rhs);

	// Assignment operator
	MeshData& operator=(const MeshData& rhs);

	// Read node data
	void readNodeData( const std::string& rootName, const bool rotateModel );

	// Read element data
	void readElemData( const std::string& rootName );

	// Read neighbor elements
	void readNeighElements( const std::string& rootName );

	// Read face data
	void readFaceData( const std::string& rootName );

	// Read releationship between boundary of the model and attributes of faces
	void readBoundaryAttribute();

	// Read releationship between resistivity values and region attributes
	void readResisitivity();

	// Include sedimentary layers
	void includeSedimentaryLayers();

	// Read parameter file
	void readParametersOfSedimentArea( const std::string& fileName, SedimentParameter& param, std::vector< std::pair<int,int> >& attrPair ) const;

	// Find triangles above the sediment
	void findTrianglesAboveSediment( const std::vector< std::pair<int,int> >& attrPair, const XY& centerCoord, 
		const double radius, const double oblateness, const double rotation, std::vector<int>& sedimentFaceSerials ) const;

	// Decide the specified point locates in the ellipse
	bool locateInEllipse( const XY& centerCoord, const double radius, const double oblateness, const double rotation, const XY& coord ) const;

	// Decide the specified point locates in the ellipse
	bool locateInEllipse( const XY& centerCoord, const double radius, const double oblateness, const double rotation, const XYZ& coord ) const;

	// Find elements corresponding to sediment
	void findElementsInSediment( const std::vector<int>& sedimentFaceSerials, const XY& centerCoord,
		const double radius, const double depth, const double oblateness, const double rotation, std::vector<int>& elemSerials ) const;

	// Change attributes of the elements corresponding to sediment layer
	void changeAttributeOfSediment( const std::vector<int>& elemSerials, const int attrSediment, const double resistivitySediment );

	// Calculate gravity center of the specified face of Tetgen
	MeshData::XY calcGravityCenterOfTetGenFace( const int iFace ) const;

	// Calculate Z coordinate of the point below the specified face of Tetgen
	//double calcZCoordOnTetGenFace( const XYZ& point, const int iFace ) const;
	bool calcZCoordOnTetGenFace( const XYZ& point, const int iFace, double& zCoord ) const;

	// Calculate area of triangle on the XY plnae
	double calcAreaOnXYPlane( const XYZ& point1, const MeshData::XYZ& pointt2, const MeshData::XYZ& point3 ) const;

	//// Determine if the inputed point locate in the specified face of Tetgen
	//bool locateInTriangleOfXYPlane( const XYZ& point, const int iFace ) const;

	//// Determine if the inputed point locate at the left of the segment on the X-Y plane
	//bool locateLeftOfSegmentOnXYPlane( const XYZ& point, const XYZ& startPointOfSegment, const XYZ& endPointOfSegment ) const;

	// Calculate intersection point of a specified line and an X-Y plane
	void calcIntersectionPoint( const XYZ& startCoord, const XYZ& endCoord, const double zCoordOfPlane, XYZ& intersectionPoint ) const;

	// Calculate volume of an element
	double calculateVolumeOfElement( const int iElem ) const;

	// Calculate volume of tetrahedron 
	double calculateVolumeOfTetrahedron( const XYZ& coord0, const XYZ& coord1, const XYZ& coord2, const XYZ& coord3 ) const;

	// Calculate boundary faces from boundary marker
	//void calcBoundaryFaces();
	void calcBoundaryFacesFromBoundaryMarker();

	// Calculate boundary faces
	void calcBoundaryFaces();

	//// Change face data from TetGen to Femtic
	//void changeFaceData();

	//// Calculate element ID and face ID for Femtic
	//void calcElementFace( const int iFace, int& elemID, int& faceID ) const;

	// Calculate element ID and face ID for Femtic
	void calcElementFace( const int boundaryType, const int iFace, int& elemID, int& faceID ) const;

	// Calculate element ID and face ID for Femtic
	void calcElementFace( const int iFace, int& elemID, int& faceID ) const;

	// Calculate local face ID of Femtic
	int calcFaceIDOfFemtic( std::vector<int>& localNodeIDs ) const;

	// Calculate boundary type from attribute
	int calcBoundTypeFromAttribute( const int iAttr, bool& isBoundary ) const;

	// Calculate boundary type from coordinate
	int calcBoundTypeFromCoordinate( const int iFace ) const;

	// Write mesh data
	void writeMeshData() const;

	// Write mesh data to vtk file
	void writeMeshDataToVTK( const std::string& rootName ) const;

	// Write resistivity data
	void writeResisitivityData() const;

	// Make resistivity block by partitioning resistivity block by RCB
	void makeResistivityBlock( const bool divAllElements );

	// Divide elements into sub-elements
	void divideElementVecs( ElementVecs& elementVecs, const int ndiv ) const;

	// Assign each element to different parameter cell
	void divideElementVecsAll( ElementVecs& elementVecs ) const;

	// Divide resistivity blocks to the required length
	void divideBlockUntilRequiredLength( ElementVecs& elementVecs ) const;

	// Divide specified elements
	std::vector<int> divide( std::vector<int>& elements ) const;

	// Calculate shorter block length of elements
	double calcShorterBlockLength( const int elem1, const int elem2 ) const;

	// Calculate maximum block length in the transition region
	double calcBlockLength( const XYZ& coord ) const;

	// Decide whether the spacified point locate within the specified hemi-sphere
	bool locateInSphere( const XYZ& coord, const int iSphere ) const;

	// Calculate maximum block length in the transition region
	double calcBlockLengthTransitionRegion( const XYZ& coord, const int iSphere ) const;

	// Calculate length on ellipsoid
	double calculateLengthOnEllipsoid( const double angleH, const double angleV, const int iSphere ) const;

	// Calculate maximum block length near site
	double calcBlockLengthNearSite( const XYZ& coord ) const;

	// Calculate gravity center of element
	XYZ calcGravityCenter( const int iElem ) const;

	// Calculate X coordinate of gravity center of element
	double calcGravityCenterX( const int iElem ) const;

	// Calculate Y coordinate of gravity center of element
	double calcGravityCenterY( const int iElem ) const;

	// Calculate Z coordinate of gravity center of element
	double calcGravityCenterZ( const int iElem ) const;

#ifdef _ADDITIONAL_FIXED_AREA
	// Add if the inputted faces of the earth's surface below which resistivity values are fixed
	void addIfFaceIndexAboveAdditionalFixedRegion( const int iFace );

	// Check if the inputted coordinate locates above the additional area below which resistivity values are fixed
	bool checkIfCoordInAdditionalFixedRegion( const XYZ& coord ) const;

	// Check if the inputted coordinate locates under the faces below which resistivity values are fixed
	bool checkIfCoordUnderTheFacesConstrutingAdditionalFixedRegion( const XY& coord ) const;

	// Fix resistivity in the additional area in which resistivity values are fixed
	void fixResistivityInAdditionalFixedRegion();

	// Calculate gravity centers of resistivity blocks
	XYZ calcGravityCenterOfResistivityBlock( const int blockIndex ) const;

	// Calculate outer product
	double calcOuterProduct( const MeshData::XY& startCoord, const MeshData::XY& endCoord, const MeshData::XY& coord ) const;

	// Write faces above the additional area in which resistivity values are fixed
	void writeFacesAboveAdditionalFixedRegionToVTK( const std::string& rootName ) const;
#endif

};

#endif
