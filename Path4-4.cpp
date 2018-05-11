#include "Filtering\Path\include\itkPolyLineParametricPath.h"
#include "Common\include\itkImage.h"
#include "ImageBase\include\itkImageFileReader.h"

void PolylineParametricPath1()
{
	const unsigned int Dimension = 2;
	typedef itk::Image< unsigned char, Dimension > ImageType;
	typedef itk::PolyLineParametricPath<Dimension> PathType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("F:/1.jpg");

	ImageType::ConstPointer image = reader->GetOutput();
	PathType::Pointer path = PathType::New();
	path->Initialize();
	typedef PathType::ContinuousIndexType ContinuousIndexType;
	ContinuousIndexType cindex;
	typedef ImageType::PointType ImagePointType;
	ImagePointType origin = image->GetOrigin();
	ImageType::SpacingType spacing = image->GetSpacing();
	ImageType::SizeType size = image->GetBufferedRegion().GetSize();
	ImagePointType point;
	point[0] = origin[0] + spacing[0] * size[0];
	point[1] = origin[1] + spacing[1] * size[1];
	image->TransformPhysicalPointToContinuousIndex(origin, cindex);
	path->AddVertex(cindex);
	image->TransformPhysicalPointToContinuousIndex(point, cindex);
	path->AddVertex(cindex);

}