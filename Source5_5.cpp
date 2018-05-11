#include "SpatialObjects\include\itkArrowSpatialObject.h"
#include "SpatialObjects\include\itkBlobSpatialObject.h"
#include "SpatialObjects\include\itkSpatialObjectPoint.h"
#include "SpatialObjects\include\itkCylinderSpatialObject.h"
#include "SpatialObjects\include\itkEllipseSpatialObject.h"
#include "SpatialObjects\include\itkSpatialObject.h"
#include "SpatialObjects\include\itkGaussianSpatialObject.h"
#include "SpatialObjects\include\itkGroupSpatialObject.h"
#include "SpatialObjects\include\itkImageSpatialObject.h"
#include "Common\include\itkImageRegionIterator.h"
#include "SpatialObjects\include\itkImageMaskSpatialObject.h"
#include "SpatialObjects\include\itkLandmarkSpatialObject.h"
#include "SpatialObjects\include\itkLineSpatialObject.h"
#include "SpatialObjects\include\itkSpatialObjectToImageFilter.h"
#include "SpatialObjects\include\itkMeshSpatialObject.h"
#include "SpatialObjects\include\itkSpatialObjectReader.h"
#include "SpatialObjects\include\itkSpatialObjectWriter.h"
#include "Common\include\itkDefaultDynamicMeshTraits.h"

void ArrowSpatialObject551()
{
	typedef itk::ArrowSpatialObject<3> ArrowType;
	ArrowType::Pointer myArrow = ArrowType::New();
	myArrow->SetLength(2);

	ArrowType::VectorType direction;
	direction.Fill(0);
	direction[1] = 1.0;
	myArrow->SetDirection(direction);
}

void BlobSpatialObject552()
{
	typedef itk::BlobSpatialObject<3> BlobType;
	typedef BlobType::Pointer BlobPointer;
	typedef itk::SpatialObjectPoint<3> BlobPointType;

	BlobType::PointListType list;
	for (unsigned int i = 0; i < 4; i ++) {
		BlobPointType p;
		p.SetPosition(i, i + 1, i + 2);
		p.SetRed(1);
		p.SetGreen(0);
		p.SetBlue(0);
		p.SetAlpha(1.0);
		list.push_back(p);
	}
	BlobPointer blob = BlobType::New();
	blob->GetProperty()->SetName("My Blob");
	blob->SetId(1);
	blob->SetPoints(list);

	BlobType::PointListType ptList = blob->GetPoints();
	std::cout << "the blob contains " << ptList.size();
	std::cout << " points " << std::endl;

	BlobType::PointListType::const_iterator iter = blob->GetPoints().begin();
	while (blob->GetPoints().end() != iter) {
		std::cout << "position = " << (*iter).GetPosition() << std::endl;
		std::cout << "color = " << (*iter).GetColor() << std::endl;
		++iter;
	}

}

void CylinderSpatialObject553()
{
	typedef itk::CylinderSpatialObject CylinderType;
	CylinderType::Pointer myCylinder = CylinderType::New();
	double radius = 4.0;
	myCylinder->SetRadius(radius);
	double height = 12.0;
	myCylinder->SetHeight(height);

	itk::Point<double, 3> insidePt;
	insidePt[0] = 1;
	insidePt[1] = 2;
	insidePt[2] = 0;
	std::cout << "Is my point " << insidePt << " inside the Cylinder? :"
		<< myCylinder->IsInside(insidePt) << std::endl;
	myCylinder->Print(std::cout);
};

void EllipseSpatialObject554()
{
	typedef itk::EllipseSpatialObject<3> EllipseType;
	EllipseType::Pointer myEllipse = EllipseType::New();
	EllipseType::ArrayType radius;
	for (unsigned int i = 0; i < 3; i ++) {
		radius[i] = i;
	}
	myEllipse->SetRadius(radius);
	EllipseType::ArrayType curRadius = myEllipse->GetRadius();
	std::cout << "Current radius is " << curRadius << std::endl;

	itk::Point<double, 3> insidePt;
	insidePt.Fill(0);
	if (myEllipse->IsInside(insidePt)) {
		std::cout << "the point " << insidePt
			<< " is really inside the ellipse" << std::endl;
	}
	itk::Point<double, 3> outsidePt;
	outsidePt.Fill(3.0);
	if (myEllipse->IsInside(outsidePt)) {
		std::cout << "the point " << outsidePt
			<< " is really outside the ellipse" << std::endl;
	}

	if (myEllipse->IsEvaluableAt(insidePt)) {
		std::cout << "the point " << insidePt
			<< " is evaluable at the point" << std::endl;
		double value;
		myEllipse->ValueAt(insidePt, value);
	}
	myEllipse->ComputeBoundingBox();
	EllipseType::BoundingBoxType *boundingBox = myEllipse->GetBoundingBox();
	std::cout << "the bounding box " << boundingBox->GetBounds() << std::endl;

}

void GaussianSpatialObject555()
{
	typedef itk::GaussianSpatialObject<3> GaussianType;
	GaussianType::Pointer myGaussian = GaussianType::New();

	myGaussian->SetMaximum(2);
	myGaussian->SetRadius(4.0);
	itk::Point<double, 3> pt;
	pt[0] = 1.0;
	pt[1] = 2.0;
	pt[2] = 1.0;

	double value;
	myGaussian->ValueAt(pt, value);
	std::cout << "ValueAt ( " << pt << " ) = " << value << std::endl;
}

void GroupSpatialObject556()
{
	typedef itk::GroupSpatialObject<3> GroupType;
	GroupType::Pointer myGroup = GroupType::New();

	typedef itk::EllipseSpatialObject<3> EllipseType;
	EllipseType::Pointer myEllipse = EllipseType::New();
	myEllipse->SetRadius(2.0);
	myGroup->AddSpatialObject(myEllipse);
	
	GroupType::VectorType offset;
	offset.Fill(0);
	myGroup->GetObjectToParentTransform()->SetOffset(offset);
	myGroup->ComputeObjectToWorldTransform();

	GroupType::PointType point;
	point.Fill(10);
	std::cout << "Is my point " << point << " inside? " 
		<< myGroup->IsInside(point, 2) << std::endl;

	myGroup->RemoveSpatialObject(myEllipse);

}

void ImageSpatialObject557()
{
	typedef itk::Image<short, 2 > Image;
	Image::Pointer image = Image::New();
	Image::SizeType size = { { 10, 10 } };
	Image::RegionType region;
	region.SetSize(size);
	image->SetRegions(region);
	image->Allocate();

	typedef itk::ImageRegionIterator<Image> Iterator;
	Iterator it(image, region);
	short pixelVal = 0;
	for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++pixelVal) {
		it.Set(pixelVal);
	}
	typedef itk::ImageSpatialObject<2, short> ImageSpatialObject;
	ImageSpatialObject::Pointer imageS0 = ImageSpatialObject::New();
	imageS0->SetImage(image);

	typedef itk::Point<double, 2> Point;
	Point insidePt;
	insidePt.Fill(9);
	if (imageS0->IsInside(insidePt)) {
		std::cout << "The point " << insidePt << " is inside the image." << std::endl;
	}

	double retValue;
	imageS0->ValueAt(insidePt, retValue);
	std::cout << "ValueAt ( " << insidePt << " ) = "
		<< retValue << std::endl;

	ImageSpatialObject::OutputVectorType retDerivative;
	imageS0->DerivativeAt(insidePt, 1, retDerivative);
	std::cout << "1st derivative at " << insidePt
		<< " = " << retDerivative << std::endl;

}

void ImageMaskSpatialObject558()
{
	typedef itk::ImageMaskSpatialObject<3> ImageMaskSpatialObject;
	typedef ImageMaskSpatialObject::PixelType PixelType;
	typedef ImageMaskSpatialObject::ImageType ImageType;
	typedef itk::ImageRegionIterator<ImageType> Iterator;
	ImageType::Pointer image = ImageType::New();
	ImageType::SizeType size = { { 50, 50, 50 } };
	ImageType::IndexType index = { { 0, 0, 0 } };
	ImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(index);
	image->SetRegions(region);
	image->Allocate(true);
	ImageType::RegionType insideRegion;
	ImageType::SizeType insideSize = { { 30, 30, 30 } };
	ImageType::IndexType insideIndex = { { 10, 10, 10 } };
	insideRegion.SetSize(size);
	insideRegion.SetIndex(index);
	Iterator iter(image, insideRegion);
	iter.GoToBegin();
	while (!iter.IsAtEnd()) {
		iter.Set(itk::NumericTraits<PixelType>::max());
		++iter;
	}
	ImageMaskSpatialObject::Pointer maskS0 = ImageMaskSpatialObject::New();
	maskS0->SetImage(image);
	ImageMaskSpatialObject::PointType inside;
	inside.Fill(20);
	std::cout << "Is my point " << inside << " inside my mask ? " 
		<< maskS0->IsInside(inside) << std::endl;
	ImageMaskSpatialObject::PointType outside;
	outside.Fill(50);
	std::cout << "Is my point " << outside << " outside my mask ? "
		<< maskS0->IsInside(outside) << std::endl;

}

void LandmarkSpatialObject559()
{
	typedef itk::LandmarkSpatialObject<3> LandmarkType;
	typedef LandmarkType::Pointer LandmarkPointer;
	typedef itk::SpatialObjectPoint<3> LandmarkPointType;
	LandmarkPointer landmark = LandmarkType::New();
	landmark->GetProperty()->SetName("landmark1");
	landmark->SetId(1);
	LandmarkType::PointListType list;
	for (unsigned int i = 0; i < 5; i++) {
		LandmarkPointType p;
		p.SetPosition(i, i + 1, i + 2);
		p.SetColor(1, 0, 0);
		list.push_back(p);
	}
	landmark->SetPoints(list);
	size_t nPoints = landmark->GetPoints().size();
	std::cout << "Number of points is " << nPoints << std::endl;
	LandmarkType::PointListType::const_iterator it
		= landmark->GetPoints().begin();
	while (landmark->GetPoints().end () != it) {
		std::cout << "position :" << (*it).GetPosition() << std::endl;
		std::cout << "color :" << (*it).GetColor() << std::endl;
		++it;
	}

}

void LineSpatialObject5510()
{
	typedef itk::LineSpatialObject<3> LineType;
	typedef LineType::Pointer LinePointer;
	typedef itk::LineSpatialObjectPoint<3> LinePointType;
	typedef itk::CovariantVector<double, 3> VectorType;
	LinePointer line = LineType::New();

	LineType::PointListType list;
	for (unsigned int i = 0; i < 3; i ++) {
		LinePointType p;
		p.SetPosition(i, i + 2, i + 4);
		p.SetColor(1, 0, 0, 1);
		VectorType normal1;
		VectorType normal2;
		for (unsigned int j = 0; j < 3; j ++) {
			normal1[j] = j;
			normal2[j] = j * 2;
		}
		p.SetNormal(normal1, 0);
		p.SetNormal(normal2, 1);
		list.push_back(p);
	}
	line->GetProperty()->SetName("Line1");
	line->SetId(1);
	line->SetPoints(list);

	LineType::PointListType pointlist = line->GetPoints();
	std::cout << "line pt number : " << pointlist.size() << std::endl;
	LineType::PointListType::const_iterator it = line->GetPoints().begin();
	while (line->GetPoints().end() != it) {
		std::cout << "pos =" << (*it).GetPosition() << std::endl;
		std::cout << "clr =" << (*it).GetColor() << std::endl;
		std::cout << "1st normal =" << (*it).GetNormal(0) << std::endl;
		std::cout << "2nd normal =" << (*it).GetNormal(1) << std::endl;
		++it;
	}

}

void MeshSpatialObject5511()
{
	typedef itk::DefaultDynamicMeshTraits<float, 3, 3> MeshTrait;
	typedef itk::Mesh<float, 3, MeshTrait> MeshType;
	typedef MeshType::CellTraits CellTraits;
	typedef itk::CellInterface<float, CellTraits> CellInterfaceType;
	typedef itk::TetrahedronCell<CellInterfaceType> TetraCellType;
	typedef MeshType::PointType PointType;
	typedef MeshType::CellType CellType;
	typedef CellType::CellAutoPointer CellAutoPointer;
	MeshType::Pointer myMesh = MeshType::New();

	MeshType::CoordRepType testPointCoords[4][3] =
		{ { 0, 0, 0 }, { 9, 0, 0 }, { 9, 9, 0 }, { 0, 0, 9 } };
	MeshType::PointIdentifier tetraPoints[4] = { 0, 1, 2, 4 };
	for (int i = 0; i < 4; ++i) {
		myMesh->SetPoint(i, PointType(testPointCoords[i]));
	}
	myMesh->SetCellsAllocationMethod(MeshType::CellsAllocatedDynamicallyCellByCell);
	CellAutoPointer testCell1;
	testCell1.TakeOwnership(new TetraCellType);
	testCell1->SetPointIds(tetraPoints);
	myMesh->SetCell(0, testCell1);
	typedef itk::MeshSpatialObject<MeshType> MeshSpatialObjectType;
	MeshSpatialObjectType::Pointer myMeshSpatialObject = MeshSpatialObjectType::New();
	myMeshSpatialObject->SetMesh(myMesh);

	myMeshSpatialObject->GetMesh();
	std::cout << "Mesh bounds : " <<
		myMeshSpatialObject->GetBoundingBox()->GetBounds() << std::endl;
	MeshSpatialObjectType::PointType myPhysicalPoint;
	myPhysicalPoint.Fill(1);
	std::cout << "Is my physical point inside? : " <<
		myMeshSpatialObject->IsInside(myPhysicalPoint) << std::endl;
	typedef itk::SpatialObjectWriter< 3, float, MeshTrait > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(myMeshSpatialObject);
	writer->SetFileName("myMesh.meta");
	writer->Update();
	typedef itk::SpatialObjectReader< 3, float, MeshTrait > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("myMesh.meta");
	reader->Update();
	typedef itk::Image< unsigned char, 3 > ImageType;
	typedef itk::GroupSpatialObject< 3 >   GroupType;
	typedef itk::SpatialObjectToImageFilter< GroupType, ImageType >
		SpatialObjectToImageFilterType;
	SpatialObjectToImageFilterType::Pointer imageFilter =
		SpatialObjectToImageFilterType::New();
	imageFilter->SetInput(reader->GetGroup());
	imageFilter->Update();
	ImageType::Pointer myBinaryMeshImage = imageFilter->GetOutput();

}

void SurfaceSpatialObject5512()
{

}

