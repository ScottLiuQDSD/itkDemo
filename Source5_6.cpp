#include "SpatialObjects\include\itkSceneSpatialObject.h"
#include "SpatialObjects\include\itkEllipseSpatialObject.h"
#include "Common\include\itkImage.h"
#include "TestKernel\include\itkRandomImageSource.h"
#include "SpatialObjects\include\itkSpatialObjectToImageStatisticsCalculator.h"

void SceneSpatialObject56()
{
	typedef itk::SceneSpatialObject<3> SceneSpatialObjectType;
	SceneSpatialObjectType::Pointer scene = SceneSpatialObjectType::New();
	typedef itk::EllipseSpatialObject<3> EllipseType;
	EllipseType::Pointer ellipse1 = EllipseType::New();
	ellipse1->SetRadius(1);
	ellipse1->SetId(1);
	EllipseType::Pointer ellipse2 = EllipseType::New();
	ellipse2->SetRadius(2);
	ellipse2->SetId(2);
	scene->AddSpatialObject(ellipse1);
	scene->AddSpatialObject(ellipse2);
	std::cout << "Number of objects in the SceneSpatialObject = ";
	std::cout << scene->GetNumberOfObjects() << std::endl;
	std::cout << "Object in the SceneSpatialObject with an ID == 2: "
		<< std::endl;
	scene->GetObjectById(2)->Print(std::cout);
	scene->RemoveSpatialObject(ellipse1);

	SceneSpatialObjectType::ObjectListType * myObjectList = scene->GetObjects();
	std::cout << "Number of objects in the SceneSpatialObject = ";
	std::cout << myObjectList->size() << std::endl;
	scene->FixHierarchy();
	scene->Clear();

}

void SpatialObjectToImageStatistcsCalculator58()
{
	typedef itk::Image< unsigned char, 2 >      ImageType;
	typedef itk::RandomImageSource< ImageType > RandomImageSourceType;
	RandomImageSourceType::Pointer randomImageSource
		= RandomImageSourceType::New();
	ImageType::SizeValueType size[2];
	size[0] = 10;
	size[1] = 10;
	randomImageSource->SetSize(size);
	randomImageSource->Update();
	ImageType::Pointer image = randomImageSource->GetOutput();
	typedef itk::EllipseSpatialObject<2> EllipseType;
	EllipseType::Pointer ellipse = EllipseType::New();
	ellipse->SetRadius(2);
	EllipseType::VectorType offset;
	offset.Fill(5);
	ellipse->GetIndexToObjectTransform()->SetOffset(offset);
	ellipse->ComputeObjectToParentTransform();
	typedef itk::SpatialObjectToImageStatisticsCalculator<
		ImageType, EllipseType > CalculatorType;
	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetImage(image);
	calculator->SetSpatialObject(ellipse);
	calculator->Update();
	std::cout << "Sample mean = " << calculator->GetMean() << std::endl;
	std::cout << "Sample covariance = " << calculator->GetCovarianceMatrix();

}
