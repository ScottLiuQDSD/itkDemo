// testITKImageReader.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ImageBase\include\itkImageFileReader.h"
#include "Common\include\itkImage.h"
#include "Common\include\itkRGBPixel.h"
#include "Common\include\itkPointSet.h"
#include "Common\include\itkCovariantVector.h"
#include "Mesh\include\itkMesh.h"
#include "Common\include\itkLineCell.h"
#include "Common\include\itkDefaultStaticMeshTraits.h"
#include "Common\include\itkTetrahedronCell.h"
#include "Common\include\itkTriangleCell.h"
#include "Mesh\include\itkAutomaticTopologyMeshSource.h"
#include "Common\include\itkCellInterfaceVisitor.h"
#include "Filtering\Path\include\itkPolyLineParametricPath.h"

/*
///Region 4.3.9 4.3.10
typedef float PixelType;
typedef itk::Mesh<PixelType, 3> MeshType;
typedef MeshType::PointType PointType;
typedef MeshType::CellType CellType;
typedef itk::VertexCell<CellType> VertexType;
typedef itk::LineCell<CellType> LineType;
typedef itk::TriangleCell<CellType> TriangleType;
typedef itk::TetrahedronCell<CellType> TetrahedronType;

/// for 4.3.9
class CustomTriangleVisitor
{
public:
	typedef itk::TriangleCell<CellType> TriangleType;
	void Visit(unsigned long cellId, TriangleType *t)
	{
		std::cout << "Cell = " << cellId << "is a TriangleType.";
		std::cout << t->GetNumberOfPoints() << std::endl;
	}
	CustomTriangleVisitor() {}
	virtual ~CustomTriangleVisitor() {}
};
class CustomVertexVisitor
{
public:
	void Visit(unsigned long cellId, VertexType *t)
	{
		std::cout << "Cell = " << cellId << "is a Vertex.";
		std::cout << " associated with point id = ";
		std::cout << t->GetPointId() << std::endl;
	}
	CustomVertexVisitor() {}
	virtual ~CustomVertexVisitor() {}
};
class CustomLineVisitor
{
public:
	CustomLineVisitor():m_Mesh(nullptr) {}
	virtual ~CustomLineVisitor() {}
	void SetMesh(MeshType* mesh) { m_Mesh = mesh; }
	void Visit(unsigned long cellId, LineType *t)
	{
		std::cout << "Cell = " << cellId << "is a Line.";
		LineType::PointIdIterator pit = t->PointIdsBegin();
		MeshType::PointType p0;
		MeshType::PointType p1;
		m_Mesh->GetPoint(*pit++, &p0);
		m_Mesh->GetPoint(*pit++, &p1);
		const double length = p0.EuclideanDistanceTo(p1);
		std::cout << " length = " << length << std::endl;
	}
private:
	MeshType::Pointer m_Mesh;
};
/// for 4.3.10
class CustomTriangleVisitor
{
public:
	void Visit(unsigned long cellId, TriangleType *t)
	{
		std::cout << "Cell = " << cellId << "is a TriangleType.";
		LineType::PointIdIterator pit = t->PointIdsBegin();
		LineType::PointIdIterator end = t->PointIdsEnd();
		while (end != pit ) {
			std::cout << "point id = " << *pit << std::endl;
			++pit;
		}
	}
	CustomTriangleVisitor() {}
	virtual ~CustomTriangleVisitor() {}
};
class CustomTetrahedronVisitor
{
public:
	void Visit(unsigned long cellId, TetrahedronType *t)
	{
		std::cout << "Cell = " << cellId << "is a Tetrahedron. ";
		std::cout << " number of faces = ";
		std::cout << t->GetNumberOfFaces() << std::endl;
	}
	CustomTetrahedronVisitor() {}
	virtual ~CustomTetrahedronVisitor() {}
};

typedef itk::CellInterfaceVisitorImplementation<
	PixelType,
	MeshType::CellTraits,
	VertexType,
	CustomVertexVisitor
>VertexVisitorInterfaceType;
typedef itk::CellInterfaceVisitorImplementation<
	PixelType,
	MeshType::CellTraits,
	LineType,
	CustomLineVisitor
>LineVisitorInterfaceType;
typedef itk::CellInterfaceVisitorImplementation<
	PixelType,
	MeshType::CellTraits,
	TriangleType,
	CustomTriangleVisitor
	>TriangleVisitorInterfaceType;
typedef itk::CellInterfaceVisitorImplementation<
	PixelType,
	MeshType::CellTraits,
	TetrahedronType,
	CustomTetrahedronVisitor
>TetrahedronVisitorInterfaceType;
*/
void CompositeExample8();

int _tmain(int argc, _TCHAR* argv[])
{
	CompositeExample8();
	/// region 0
	/*
	typedef unsigned char PixelType;
	const unsigned int Dimension = 2;
	typedef itk::Image<PixelType, Dimension> ImageType;

	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	const char* fileName = "F:/11.bmp";
	reader->SetFileName(fileName);
	reader->Update();
	ImageType::Pointer image = reader->GetOutput();

	const ImageType::IndexType pixelIndex = { { 27, 39 } };
	ImageType::PixelType pixelValue = image->GetPixel(pixelIndex);
	image->SetPixel(pixelIndex, pixelValue + 3);

	ImageType::DirectionType direction;
	direction.SetIdentity();
	image->SetDirection(direction);

	typedef itk::Point<double, ImageType::ImageDimension> PointType;
	PointType point;
	point[0] = 1.45;
	point[1] = 7.21;
	point[2] = 9.28;

	ImageType::IndexType pixelIndex0;
	const bool isInside = image->TransformPhysicalPointToIndex(point, pixelIndex0);
	if (true == isInside) {
		ImageType::PixelType pixelVal = image->GetPixel(pixelIndex0);
		pixelVal += 5;
		image->SetPixel(pixelIndex0, pixelVal);
	}
	typedef itk::Vector<float, 3> PixelType;
	*/

	///Region 1
	/*
	typedef itk::RGBPixel<float> PixelType;
	typedef itk::PointSet<PixelType, 3> PointSetType;

	PointSetType::Pointer pointSet = PointSetType::New();

	PointSetType::PixelType pixel;
	PointSetType::PointType point;

	unsigned int pointId = 0;
	const double  radius = 3.0f;
	for (unsigned int i = 0; i < 360; i ++) {
		const double angle = i * itk::Math::pi / 180.0;
		point[0] = radius * std::sin(angle);
		point[0] = radius * std::cos(angle);
		point[2] = 1.0;

		pixel.SetRed(point[0] * 2.0);
		pixel.SetGreen(point[1] * 2.0);
		pixel.SetBlue(point[2] * 2.0);
		pointSet->SetPoint(pointId, point);
		pointSet->SetPointData(pointId, pixel);
		pointId++;

	}

	typedef PointSetType::PointsContainer::ConstIterator PointIterator;
	PointIterator ptIter = pointSet->GetPoints()->Begin();
	PointIterator ptEnd = pointSet->GetPoints()->End();

	while (ptEnd != ptIter) {
		point = ptIter.Value();
		std::cout << point << std::endl;
		++ptIter;
	}

	typedef PointSetType::PointDataContainer::ConstIterator PointDataIterator;
	PointDataIterator pxIter = pointSet->GetPointData()->Begin();
	PointDataIterator pxEnd = pointSet->GetPointData()->End();
	while (pxEnd != pxIter) {
		pixel = pxIter.Value();
		std::cout << pixel << std::endl;
		++pxIter;
	}
	*/

	///Region 2
	/*
	const unsigned int Dimension = 3;
	typedef itk::Vector<float, Dimension> PixelType;
	typedef itk::PointSet<PixelType, Dimension> PointSetType;

	PointSetType::Pointer pointSet = PointSetType::New();

	PointSetType::PixelType tangent;
	PointSetType::PointType point;

	unsigned int pointId = 0;
	const double  radius = 300.0f;
	for (unsigned int i = 0; i < 360; i ++) {
	const double angle = i * itk::Math::pi / 180.0;
	point[0] = radius * std::sin(angle);
	point[0] = radius * std::cos(angle);
	point[2] = 1.0;

	tangent[0] = std::cos(angle);
	tangent[1] = -std::sin(angle);
	tangent[2] = 0.0;
	pointSet->SetPoint(pointId, point);
	pointSet->SetPointData(pointId, tangent);
	pointId++;

	}

	typedef PointSetType::PointDataContainer::ConstIterator PointDataIterator;
	PointDataIterator pxIter = pointSet->GetPointData()->Begin();
	PointDataIterator pxEnd = pointSet->GetPointData()->End();
	typedef PointSetType::PointsContainer::ConstIterator PointIterator;
	PointIterator ptIter = pointSet->GetPoints()->Begin();
	PointIterator ptEnd = pointSet->GetPoints()->End();
	while (pxEnd != pxIter && ptEnd != ptIter) {
		ptIter->Value() = itk::Point<float>(ptIter.Value() + pxIter.Value());
		++pxIter;
		++ptIter;
	}
	*/

	///Region 3 normal vector as pixel
	/*
	const unsigned int Dimension = 3;
	typedef itk::CovariantVector<float, Dimension> PixelType;
	typedef itk::PointSet<PixelType, Dimension> PointSetType;

	PointSetType::Pointer pointSet = PointSetType::New();

	PointSetType::PixelType gradient;
	PointSetType::PointType point;

	unsigned int pointId = 0;
	const double  radius = 300.0f;
	for (unsigned int i = 0; i < 360; i ++) {
		const double angle = i * itk::Math::pi / 180.0;
		point[0] = radius * std::sin(angle);
		point[0] = radius * std::cos(angle);
		point[2] = 1.0;

		gradient[0] = std::sin(angle);
		gradient[1] = std::cos(angle);
		gradient[2] = 0.0;
		pointSet->SetPoint(pointId, point);
		pointSet->SetPointData(pointId, gradient);
		pointId++;

	}

	typedef PointSetType::PointDataContainer::ConstIterator PointDataIterator;
	PointDataIterator pxIter = pointSet->GetPointData()->Begin();
	PointDataIterator pxEnd = pointSet->GetPointData()->End();
	typedef PointSetType::PointsContainer::ConstIterator PointIterator;
	PointIterator ptIter = pointSet->GetPoints()->Begin();
	PointIterator ptEnd = pointSet->GetPoints()->End();
	int idx = 0;
	while (pxEnd != pxIter && ptEnd != ptIter) {
		point = ptIter.Value();
		gradient = pxIter.Value();
		for (unsigned int i = 0; i < Dimension; i ++) {
			point[i] += gradient[i];
		}
		pointSet->GetPoints()->SetElement(idx, point);
		++pxIter;
		++ptIter;
		++idx;
	}
	*/

	/// Region 4 Mesh
	/*
	const unsigned int Dimension = 3;
	typedef float PixelType;
	typedef itk::Mesh<PixelType, Dimension> MeshType;

	MeshType::Pointer mesh = MeshType::New();

	MeshType::PointType point0;
	MeshType::PointType point1;
	MeshType::PointType point2;
	MeshType::PointType point3;
	point0[0] = -1.0, point0[1] = -1.0, point0[2] = 0.0;
	point1[0] = 1.0, point1[1] = -1.0, point1[2] = 0.0;
	point2[0] = 1.0, point2[1] = 1.0, point2[2] = 0.0;
	point3[0] = -1.0, point3[1] = 1.0, point3[2] = 0.0;
	mesh->SetPoint(0, point0);
	mesh->SetPoint(1, point1);
	mesh->SetPoint(2, point2);
	mesh->SetPoint(3, point3);

	std::cout << "Points = " << mesh->GetNumberOfPoints() << std::endl;


	typedef MeshType::PointsContainer::Iterator PointsIterator;
	PointsIterator ptIter = mesh->GetPoints()->Begin();
	PointsIterator ptEnd = mesh->GetPoints()->End();
	while (ptEnd != ptIter) {
		MeshType::PointType p = ptIter.Value();
		std::cout << p << std::endl;
		++ptIter;
	}
	*/

	/// Region 5 Mesh insert unit
	/*
	const unsigned int Dimension = 3;
	typedef float PixelType;
	typedef itk::Mesh<PixelType, Dimension> MeshType;
	typedef MeshType::CellType CellType;
	typedef itk::LineCell<CellType> LineType;
	typedef CellType::CellAutoPointer CellAutoPointer;


	MeshType::Pointer mesh = MeshType::New();

	MeshType::PointType point0;
	MeshType::PointType point1;
	MeshType::PointType point2;
	MeshType::PointType point3;
	point0[0] = -1.0, point0[1] = -1.0, point0[2] = 0.0;
	point1[0] = 1.0, point1[1] = -1.0, point1[2] = 0.0;
	point2[0] = 1.0, point2[1] = 1.0, point2[2] = 0.0;
	point3[0] = -1.0, point3[1] = 1.0, point3[2] = 0.0;
	mesh->SetPoint(0, point0);
	mesh->SetPoint(1, point1);
	mesh->SetPoint(2, point2);
	mesh->SetPoint(3, point3);

	CellAutoPointer line0;
	CellAutoPointer line1;
	CellAutoPointer line2;
	CellAutoPointer line3;

	line0.TakeNoOwnership(new LineType);
	line1.TakeNoOwnership(new LineType);
	line2.TakeNoOwnership(new LineType);
	line3.TakeNoOwnership(new LineType);

	line0->SetPointId(0, 0);
	line0->SetPointId(1, 1);
	line1->SetPointId(0, 1);
	line1->SetPointId(1, 2);
	line2->SetPointId(0, 2);
	line2->SetPointId(1, 3);
	line3->SetPointId(0, 3);
	line3->SetPointId(1, 0);

	mesh->SetCell(0, line0);
	mesh->SetCell(1, line1);
	mesh->SetCell(2, line2);
	mesh->SetCell(3, line3);

	std::cout << "Cells = " << mesh->GetNumberOfCells() << std::endl;


	typedef MeshType::CellsContainer::Iterator PointsIterator;
	PointsIterator clIter = mesh->GetCells()->Begin();
	PointsIterator clEnd = mesh->GetCells()->End();
	while (clEnd != clIter) {
		MeshType::CellType* cellPtr = clIter.Value();
		LineType *line = dynamic_cast<LineType *>(cellPtr);
		if (ITK_NULLPTR == line) {
			continue;
		}
		std::cout << line->GetNumberOfPoints() << std::endl;
		++clIter;
	}
	*/
	
	/// Region 4.3.3 Mesh manage data in unit
	/*
	const unsigned int Dimension = 2;
	typedef float PixelType;
	typedef itk::Mesh<PixelType, Dimension> MeshType;
	typedef MeshType::CellType CellType;
	typedef itk::LineCell<CellType> LineType;
	typedef CellType::CellAutoPointer CellAutoPointer;

	MeshType::Pointer mesh = MeshType::New();
	typedef MeshType::PointType PointType;

	PointType point;
	const unsigned int numberOfPoints = 10;
	for (unsigned int id = 0; id < numberOfPoints; id ++) {
		point[0] = static_cast<PointType::ValueType>(id);
		point[1] = std::log(static_cast<double>(id) +itk::Math::eps);
		mesh->SetPoint(id, point);
	}
	// not a suggest way
	CellType::CellAutoPointer line;
	const unsigned int numberOfCells = numberOfPoints - 1;
	for (unsigned int cellId = 0; cellId < numberOfCells; cellId ++) {
		line.TakeNoOwnership(new LineType);
		line->SetPointId(0, cellId);
		line->SetPointId(1, cellId + 1);
		mesh->SetCell(cellId, line);
		mesh->SetCellData(cellId, static_cast<PixelType>(cellId*cellId)); // not effective
	}
	for (unsigned int cellId = 0; cellId < numberOfCells; cellId++) {
		PixelType value = static_cast<PixelType>(0.0);
		mesh->GetCellData(cellId, &value);
		std::cout << "Cell " << cellId << " = " << value << std::endl;
	}

	// the way suggested
	typedef MeshType::CellDataContainer::Iterator PointsIterator;
	PointsIterator clIter = mesh->GetCellData()->Begin();
	PointsIterator clEnd = mesh->GetCellData()->End();
	while (clEnd != clIter) {
		PixelType value = clIter.Value();
		std::cout << line->GetNumberOfPoints() << std::endl;
		++clIter;
	}
	*/

	/// Region 4.3.4 custom mesh
	/*
	const unsigned int PointDimension = 3;
	const unsigned int MaxTopologicalDimension = 2;
	typedef itk::Vector<double, 4> PixelType;
	typedef itk::Matrix<double, 4, 3> CellDataType;
	typedef double CoordinateType;
	typedef double InterpolationWeightType;
	typedef itk::DefaultStaticMeshTraits<
		PixelType, PointDimension, MaxTopologicalDimension,
		CoordinateType, InterpolationWeightType, CellDataType> MeshTraits;
	typedef itk::Mesh<PixelType, PointDimension, MeshTraits> MeshType;
	typedef MeshType::CellType CellType;
	typedef itk::LineCell<CellType> LineType;
	typedef CellType::CellAutoPointer CellAutoPointer;

	MeshType::Pointer mesh = MeshType::New();
	typedef MeshType::PointType PointType;
	PointType point;
	const unsigned int numberOfPoints = 10;
	for (unsigned int id = 0; id < numberOfPoints; id ++) {
		point[0] = 1.565;
		point[1] = 3.647;
		point[2] = 4.129;
		mesh->SetPoint(id, point);
	}
	// not a suggest way
	CellType::CellAutoPointer line;
	const unsigned int numberOfCells = numberOfPoints - 1;
	for (unsigned int cellId = 0; cellId < numberOfCells; cellId ++) {
		line.TakeNoOwnership(new LineType);
		line->SetPointId(0, cellId);
		line->SetPointId(1, cellId + 1);
		mesh->SetCell(cellId, line);
		CellDataType value;
		mesh->SetCellData(cellId, value);
	}
	for (unsigned int cellId = 0; cellId < numberOfCells; cellId++) {
		CellDataType value;
		mesh->GetCellData(cellId, &value);
		std::cout << "Cell " << cellId << " = " << value << std::endl;
	}

	// the way suggested
	typedef MeshType::CellDataContainer::Iterator PointsIterator;
	PointsIterator clIter = mesh->GetCellData()->Begin();
	PointsIterator clEnd = mesh->GetCellData()->End();
	while (clEnd != clIter) {
		CellDataType value = clIter.Value();
		std::cout << line->GetNumberOfPoints() << std::endl;
		++clIter;
	}
	*/

	/// Region 4.3.5 Topological & K-Complex wave
	/*
	const unsigned int Dimension = 3;
	typedef float PixelType;
	typedef itk::Mesh<PixelType, Dimension> MeshType;
	typedef MeshType::CellType CellType;
	typedef itk::VertexCell<CellType> VertexType;
	typedef itk::LineCell<CellType> LineType;
	typedef itk::TriangleCell<CellType> TriangleCell;
	typedef itk::TetrahedronCell<CellType> TetrahedronType;

	MeshType::Pointer mesh = MeshType::New();
	MeshType::PointType point0;
	MeshType::PointType point1;
	MeshType::PointType point2;
	MeshType::PointType point3;
	point0[0] = -1.0, point0[1] = -1.0, point0[2] = 0.0;
	point1[0] = 1.0, point1[1] = -1.0, point1[2] = 0.0;
	point2[0] = 1.0, point2[1] = 1.0, point2[2] = 0.0;
	point3[0] = -1.0, point3[1] = 1.0, point3[2] = 0.0;
	mesh->SetPoint(0, point0);
	mesh->SetPoint(1, point1);
	mesh->SetPoint(2, point2);
	mesh->SetPoint(3, point3);

	CellType::CellAutoPointer cellPointer;
	cellPointer.TakeNoOwnership(new TetrahedronType);
	cellPointer->SetPointId(0, 0);
	cellPointer->SetPointId(1, 1);
	cellPointer->SetPointId(2, 2);
	cellPointer->SetPointId(3, 3);
	mesh->SetCell(0, cellPointer);

	cellPointer.TakeNoOwnership(new TriangleCell);
	cellPointer->SetPointId(0, 0);
	cellPointer->SetPointId(1, 1);
	cellPointer->SetPointId(2, 2);
	mesh->SetCell(1, cellPointer);
	cellPointer.TakeNoOwnership(new TriangleCell);
	cellPointer->SetPointId(0, 0);
	cellPointer->SetPointId(1, 2);
	cellPointer->SetPointId(2, 3);
	mesh->SetCell(2, cellPointer);
	cellPointer.TakeNoOwnership(new TriangleCell);
	cellPointer->SetPointId(0, 0);
	cellPointer->SetPointId(1, 3);
	cellPointer->SetPointId(2, 1);
	mesh->SetCell(3, cellPointer);
	cellPointer.TakeNoOwnership(new TriangleCell);
	cellPointer->SetPointId(0, 3);
	cellPointer->SetPointId(1, 2);
	cellPointer->SetPointId(2, 1);
	mesh->SetCell(4, cellPointer);

	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 0);
	cellPointer->SetPointId(1, 1);
	mesh->SetCell(5, cellPointer);
	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 1);
	cellPointer->SetPointId(1, 2);
	mesh->SetCell(6, cellPointer);
	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 2);
	cellPointer->SetPointId(1, 0);
	mesh->SetCell(7, cellPointer);
	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 1);
	cellPointer->SetPointId(1, 3);
	mesh->SetCell(8, cellPointer);
	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 3);
	cellPointer->SetPointId(1, 2);
	mesh->SetCell(9, cellPointer);
	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 3);
	cellPointer->SetPointId(1, 0);
	mesh->SetCell(10, cellPointer);

	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 0);
	mesh->SetCell(11, cellPointer);
	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 1);
	mesh->SetCell(12, cellPointer);
	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 2);
	mesh->SetCell(13, cellPointer);
	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 3);
	mesh->SetCell(14, cellPointer);

	typedef MeshType::PointsContainer::ConstIterator PointIterator;
	PointIterator ptIter = mesh->GetPoints()->Begin();
	PointIterator ptEnd = mesh->GetPoints()->End();
	while (ptEnd != ptIter) {
		std::cout << ptIter.Value() << std::endl;
		++ptIter;
	}
	
	typedef MeshType::CellsContainer::ConstIterator CellIterator;
	CellIterator cellIter = mesh->GetCells()->Begin();
	CellIterator cellEnd = mesh->GetCells()->End();
	while (cellEnd != cellIter) {
		CellType *cell = cellIter.Value();
		std::cout << cell->GetNumberOfPoints() << std::endl;
		typedef CellType::PointIdIterator PointIdIterator;
		PointIdIterator ptIter = cell->PointIdsBegin();
		PointIdIterator ptEnd = cell->PointIdsEnd();
		while (ptEnd != ptIter) {
			std::cout << *ptIter << std::endl;
			++ptIter;
		}
		++cellIter;
	}

	/// this should be changed if we want identify all cells
	/// here is only a demo, do not do thing like below as cellId keep as zero.
	MeshType::CellIdentifier cellId = 0;
	int dimension = 0;
	MeshType::CellFeatureIdentifier featureId = 0;
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 11);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 12);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 13);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 14);

	cellId = 0;
	dimension = 1;
	featureId = 0;
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 5);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 6);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 7);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 8);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 9);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 10);

	cellId = 0;
	dimension = 2;
	featureId = 0;
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 1);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 2);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 3);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 4);

	cellId = 0;
	MeshType::CellFeatureCount n0;
	MeshType::CellFeatureCount n1;
	MeshType::CellFeatureCount n2;
	n0 = mesh->GetNumberOfCellBoundaryFeatures(0, cellId);
	n1 = mesh->GetNumberOfCellBoundaryFeatures(1, cellId);
	n2 = mesh->GetNumberOfCellBoundaryFeatures(2, cellId);

	dimension = 0;
	for (unsigned int b0 = 0; b0 < n0; b0 ++) {
		MeshType::CellIdentifier id;
		bool found = mesh->GetBoundaryAssignment(dimension, cellId, b0, &id);
		if (true == found) {
			std::cout << id << std::endl;		}
	}

	cellId = 2;
	dimension = 1;
	featureId = 0;
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 7);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 9);
	mesh->SetBoundaryAssignment(dimension, cellId, featureId++, 10);
	*/


	/// Region 4.3.6 polyline
	/*
	const unsigned int Dimension = 2;
	typedef float PixelType;
	typedef itk::Mesh<PixelType, Dimension> MeshType;
	typedef MeshType::CellType CellType;
	typedef itk::VertexCell<CellType> VertexType;
	typedef itk::LineCell<CellType> LineType;
	typedef CellType::CellAutoPointer CellAutoPointer;

	MeshType::Pointer mesh = MeshType::New();

	MeshType::PointType point0;
	MeshType::PointType point1;
	MeshType::PointType point2;
	MeshType::PointType point3;
	point0[0] = -1.0, point0[1] = -1.0;
	point1[0] = 1.0, point1[1] = -1.0;
	point2[0] = 1.0, point2[1] = 1.0;
	point3[0] = -1.0, point3[1] = 1.0;
	mesh->SetPoint(0, point0);
	mesh->SetPoint(1, point1);
	mesh->SetPoint(2, point2);
	mesh->SetPoint(3, point3);

	CellAutoPointer cellPointer;

	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 0);
	cellPointer->SetPointId(1, 1);
	mesh->SetCell(0, cellPointer);

	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 1);
	cellPointer->SetPointId(1, 2);
	mesh->SetCell(1, cellPointer);

	cellPointer.TakeNoOwnership(new LineType);
	cellPointer->SetPointId(0, 2);
	cellPointer->SetPointId(1, 0);
	mesh->SetCell(2, cellPointer);

	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 0);
	mesh->SetCell(3, cellPointer);

	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 1);
	mesh->SetCell(4, cellPointer);

	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 2);
	mesh->SetCell(5, cellPointer);

	cellPointer.TakeNoOwnership(new VertexType);
	cellPointer->SetPointId(0, 3);
	mesh->SetCell(6, cellPointer);

	std::cout << "Cells = " << mesh->GetNumberOfCells() << std::endl;
	
	typedef MeshType::PointsContainer::Iterator PointsIterator;
	PointsIterator ptIter = mesh->GetPoints()->Begin();
	PointsIterator ptEnd = mesh->GetPoints()->End();
	while (ptEnd != ptIter) {
		std::cout << ptIter.Value() << std::endl;
		++ptIter;
	}

	typedef MeshType::CellsContainer::Iterator CellIterator;
	CellIterator clIter = mesh->GetCells()->Begin();
	CellIterator clEnd = mesh->GetCells()->End();
	while (clEnd != clIter) {
		CellType *cell = clIter.Value();
		std::cout << cell->GetNumberOfPoints() << std::endl;
		++clIter;
	}
	*/
	
	/// Region 4.3.7 Automatic Mesh
	/*
	typedef float PixelType;
	typedef itk::Mesh<PixelType, 3> MeshType;
	typedef MeshType::PointType PointType;
	typedef itk::AutomaticTopologyMeshSource<MeshType> MeshSourceType;
	typedef MeshSourceType::IdentifierArrayType IdentifierArrayType;
	MeshSourceType::Pointer meshSource;
	meshSource = MeshSourceType::New();

	meshSource->AddTetrahedron(
		meshSource->AddPoint(-1, -1, -1),
		meshSource->AddPoint(1, 1, -1),
		meshSource->AddPoint(1, -1, 1),
		meshSource->AddPoint(-1, 1, 1)
		);
	*/

	/// Region 4.3.8 mesh cells iteration
	/*
	typedef float PixelType;
	typedef itk::Mesh<PixelType, 3> MeshType;
	typedef MeshType::PointType PointType;
	typedef itk::AutomaticTopologyMeshSource<MeshType> MeshSourceType;
	typedef MeshSourceType::IdentifierArrayType IdentifierArrayType;
	MeshSourceType::Pointer meshSource;
	meshSource = MeshSourceType::New();

	meshSource->AddTetrahedron(
		meshSource->AddPoint(-1, -1, -1),
		meshSource->AddPoint(1, 1, -1),
		meshSource->AddPoint(1, -1, 1),
		meshSource->AddPoint(-1, 1, 1)
		);
	typedef MeshType::CellsContainer::ConstIterator CellIterator;
	typedef MeshType::CellType CellType;
	typedef itk::VertexCell<CellType> VertexType;
	typedef itk::LineCell<CellType> LineType;
	typedef itk::TriangleCell<CellType> TriangleType;
	MeshType *mesh = meshSource->GetOutput();
	CellIterator cellIter = mesh->GetCells()->Begin();
	CellIterator cellEnd = mesh->GetCells()->End();
	while (cellEnd != cellIter) {
		CellType *cell = cellIter->Value();
		std::cout << cell->GetNumberOfPoints() << std::endl;
// 		if (CellType::LINE_CELL == cell->GetType()) {
// 			LineType *line = static_cast<LineType *>(cell);
// 			std::cout << "dimension" << line->GetDimension();
// 			std::cout << "# points = " << line->GetNumberOfPoints();
// 			std::cout << std::endl;
// 		}
		switch (cell->GetType()) {
			case CellType::VERTEX_CELL:
			{
				std::cout << "VertexCell: " << std::endl;
				VertexType *vertex = static_cast<VertexType *>(cell);
				std::cout << "dimension" << vertex->GetDimension();
				std::cout << "# points = " << vertex->GetNumberOfPoints();
				std::cout << std::endl;
				break;
			}
			case CellType::LINE_CELL:
			{
				std::cout << "LineCell: " << std::endl;
				LineType *line = static_cast<LineType *>(cell);
				std::cout << "dimension" << line->GetDimension();
				std::cout << "# points = " << line->GetNumberOfPoints();
				std::cout << std::endl;
				break;
			}
			case CellType::TRIANGLE_CELL:
			{
				std::cout << "LineCell: " << std::endl;
				TriangleType *trianlge = static_cast<TriangleType *>(cell);
				std::cout << "dimension" << trianlge->GetDimension();
				std::cout << "# points = " << trianlge->GetNumberOfPoints();
				std::cout << std::endl;
				break;
			}
			default:
			{
				std::cout << "Cell with more than 3 points : " << std::endl;
				std::cout << "dimension" << cell->GetDimension();
				std::cout << "# points = " << cell->GetNumberOfPoints();
				std::cout << std::endl;
				break;
			}
		}
		++cellIter;
	}
	*/

	/// Region 4.3.9 access cells mesh cell visitor
	/*
	typedef float PixelType;
	typedef itk::Mesh<PixelType, 3> MeshType;
	typedef MeshType::PointType PointType;
	typedef itk::AutomaticTopologyMeshSource<MeshType> MeshSourceType;
	typedef MeshSourceType::IdentifierArrayType IdentifierArrayType;
	MeshSourceType::Pointer meshSource;
	meshSource = MeshSourceType::New();

	meshSource->AddTetrahedron(
		meshSource->AddPoint(-1, -1, -1),
		meshSource->AddPoint(1, 1, -1),
		meshSource->AddPoint(1, -1, 1),
		meshSource->AddPoint(-1, 1, 1)
		);

	MeshType *mesh = meshSource->GetOutput();
	TriangleVisitorInterfaceType::Pointer triangleVisitor = TriangleVisitorInterfaceType::New();
	typedef CellType::MultiVisitor CellMultiVisitorType;
	CellMultiVisitorType::Pointer multiVisitor = CellMultiVisitorType::New();
	multiVisitor->AddVisitor(triangleVisitor);
	mesh->Accept(multiVisitor);
	*/

	/// Region 4.3.10 visitor more
	/*
	typedef float PixelType;
	typedef itk::Mesh<PixelType, 3> MeshType;
	typedef MeshType::PointType PointType;
	typedef itk::AutomaticTopologyMeshSource<MeshType> MeshSourceType;
	typedef MeshSourceType::IdentifierArrayType IdentifierArrayType;
	MeshSourceType::Pointer meshSource;
	meshSource = MeshSourceType::New();

	meshSource->AddTetrahedron(
		meshSource->AddPoint(-1, -1, -1),
		meshSource->AddPoint(1, 1, -1),
		meshSource->AddPoint(1, -1, 1),
		meshSource->AddPoint(-1, 1, 1)
		);

	MeshType *mesh = meshSource->GetOutput();
	VertexVisitorInterfaceType::Pointer vertexVisitor = VertexVisitorInterfaceType::New();
	LineVisitorInterfaceType::Pointer lineVisitor = LineVisitorInterfaceType::New();
	TriangleVisitorInterfaceType::Pointer triangleVisitor = TriangleVisitorInterfaceType::New();
	TriangleVisitorInterfaceType::Pointer tetrahedromVisitor = TriangleVisitorInterfaceType::New();
	lineVisitor->SetMesh(mesh);

	typedef CellType::MultiVisitor CellMultiVisitorType;
	CellMultiVisitorType::Pointer multiVisitor = CellMultiVisitorType::New();
	multiVisitor->AddVisitor(vertexVisitor);
	multiVisitor->AddVisitor(lineVisitor);
	multiVisitor->AddVisitor(triangleVisitor);
	multiVisitor->AddVisitor(tetrahedromVisitor);
	mesh->Accept(multiVisitor);
	*/
















return 0;
}

