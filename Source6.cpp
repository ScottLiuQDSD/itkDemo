#include "Common\include\itkImage.h"
#include "Common\include\itkImageRegionConstIterator.h"
#include "Common\include\itkImageRegionIterator.h"
#include "SpatialObjects\include\itkGroupSpatialObject.h"
#include "SpatialObjects\include\itkSpatialObjectToImageFilter.h"
#include "ImageBase\include\itkImageFileReader.h"
#include "ImageBase\include\itkImageFileWriter.h"
#include "Common\include\itkImageRegionIteratorWithIndex.h"
#include "Common\include\itkImageLinearIteratorWithIndex.h"
#include "Common\include\itkImageSliceConstIteratorWithIndex.h"
#include "Common\include\itkImageSliceIteratorWithIndex.h"
#include "Common\include\itkImageRandomConstIteratorWithIndex.h"
#include "Common\include\itkNeighborhoodIterator.h"
#include "Common\include\itkConstNeighborhoodIterator.h"
#include "Filtering\ImageIntensity\include\itkRescaleIntensityImageFilter.h"
#include "Common\include\itkSobelOperator.h"
#include "Common\include\itkNeighborhoodInnerProduct.h"
#include "Common\include\itkNeighborhoodAlgorithm.h"
#include "Common\include\itkGaussianOperator.h"
#include "Filtering\FastMarching\include\itkFastMarchingImageFilter.h"
#include "Filtering\ImageIntensity\include\itkAddImageFilter.h"
#include "TestKernel\include\itkRandomImageSource.h"
#include "Common\include\itkConstShapedNeighborhoodIterator.h"

void CreateIterator6_2_1()
{

	typedef itk::Image<float, 3> ImageType;
	typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	typedef itk::GroupSpatialObject< 3 >   GroupType;
	typedef itk::SpatialObjectToImageFilter< GroupType, ImageType >
		SpatialObjectToImageFilterType;
	SpatialObjectToImageFilterType::Pointer imageFilter =
		SpatialObjectToImageFilterType::New();
	imageFilter->Update();

	ImageType::Pointer image = imageFilter->GetOutput();
	ConstIteratorType constIter(image, image->GetRequestedRegion());
	IteratorType iter(image, image->GetRequestedRegion());

}
void MoveIterator622()
{
	// GotoBegin(), GotoEnd(), ++, --
	// +=, -=, SetPosition()
	//bool isAtEnd();
	//bool isAtBegin();
	//IndexType GetIndex();

}
void AccessData623()
{
	//PixelType Get();
	//void Set(PixelType);
	// aviod it.Set(it.Get() +1);
	// PixelType& Value();
	// it.Value()++;
}
void LoopIterator624()
{
// 	ConstIteratorType in(inputImage, inputImage->GetRequestedRegion());
// 	IteratorType out(outputImage, inputImage->GetRequestedRegion());
// 	for (in.GotoBegin(), out.GoToBegin(); !in.isAtEnd(); ++in, ++out) {
// 		out.Set(in.Get() * in.Get());
// 	}
// 	in.GotoEnd();
// 	out.GoToEnd();
// 	while (!in.IsAtBegin()) {
// 		--in, --out;
// 		out.Set(in.Get() * in.Get());
// 	}
}
void ImageRegionIterator631()
{
	const unsigned int Dimension = 2;
	typedef unsigned char PixelType;
	typedef itk::Image<PixelType, Dimension> ImageType;
	typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;
	
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriteType;

	ImageType::RegionType inputRegion;
	ImageType::RegionType::IndexType inputStart;
	ImageType::RegionType::SizeType size;
	inputStart[0] = 20;
	inputStart[1] = 70;
	size[0] = 210;
	size[1] = 140;
	inputRegion.SetSize(size);
	inputRegion.SetIndex(inputStart);
	
	ImageType::RegionType outputRegion;
	ImageType::RegionType::IndexType outputStart;
	outputStart[0] = 0;
	outputStart[1] = 0;
	outputRegion.SetSize(size);
	outputRegion.SetIndex(outputStart);

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\FatMRISlice.png");
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject& err)
	{
		std::cerr << "Exception caught!" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

	//check the reagion is contained within the input image or not.
	if ( ! reader->GetOutput()->GetRequestedRegion().IsInside(inputRegion)) {
		std::cerr << "Error!" << std::endl;
		std::cerr << "The region " << inputRegion << " is not contained within the input image region..." << std::endl;
	}

	ImageType::Pointer outputImage = ImageType::New();
	outputImage->SetRegions(outputRegion);
	const ImageType::SpacingType& spacing = reader->GetOutput()->GetSpacing();
	const ImageType::PointType& inputOrigin = reader->GetOutput()->GetOrigin();
	double outputOrigin[Dimension];
	for (unsigned int i = 0; i < Dimension; i++) {
		outputOrigin[i] = inputOrigin[i] + spacing[i] * inputStart[i];
	}
	outputImage->SetSpacing(spacing);
	outputImage->SetOrigin(outputOrigin);
	outputImage->Allocate();

	ConstIteratorType inputIter(reader->GetOutput(), inputRegion);
	IteratorType outputIter(outputImage, outputRegion);
	inputIter.GoToBegin();
	outputIter.GoToBegin();
	while (!inputIter.IsAtEnd()) {
		outputIter.Set(inputIter.Get());
		++inputIter;
		++outputIter;
	}

	WriteType::Pointer writer = WriteType::New();
	writer->SetFileName("G:/itkOutput/ImageData/outFatMRISlice.png");
	writer->SetInput(outputImage);

	try
	{
		writer->Update();
	}catch (itk::ExceptionObject &err) {
		std::cerr << "Exception caught!" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void ImageRegionIteratorWithIndex632()
{
	const unsigned int Dimension = 2;
	typedef itk::RGBPixel<unsigned char> RGBPixelType;
	typedef itk::Image<RGBPixelType, Dimension> ImageType;
	typedef itk::ImageRegionIteratorWithIndex<ImageType> IterWithIndexType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriteType;
	ImageType::ConstPointer inputImage;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\VisibleWomanEyeSlice.png");
	try {
		reader->Update();
		inputImage = reader->GetOutput();
	} catch (itk::ExceptionObject& err) {
		std::cerr << "Exception caught!" << std::endl;
		std::cerr << err << std::endl;
		return;
	}


	ImageType::Pointer outputImage = ImageType::New();
	outputImage->SetRegions(inputImage->GetRequestedRegion());
	outputImage->CopyInformation(inputImage);
	outputImage->Allocate();

	IterWithIndexType outputIt(outputImage, outputImage->GetRequestedRegion());
	ImageType::IndexType requestedIndex =
		outputImage->GetRequestedRegion().GetIndex();
	ImageType::SizeType requestedSize =
		outputImage->GetRequestedRegion().GetSize();
	for (outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt) {
		ImageType::IndexType idx = outputIt.GetIndex();
		idx[0] = requestedIndex[0] + requestedSize[0] - 1 - idx[0];
		outputIt.Set(inputImage->GetPixel(idx));
	}

	WriteType::Pointer writer = WriteType::New();
	writer->SetFileName("G:/itkOutput/ImageData/out632.png");
	writer->SetInput(outputImage);
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

}
void ImageLinearIteratorWithIndex633_1()
{
	const unsigned int Dimension = 2;
	typedef itk::RGBPixel<unsigned char> RGBPixelType;
	typedef itk::Image<RGBPixelType, Dimension> ImageType;
	typedef itk::ImageLinearIteratorWithIndex<ImageType> IteratorType;
	typedef itk::ImageLinearConstIteratorWithIndex<ImageType> ConstIteratorType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<ImageType> WriteType;
	ImageType::ConstPointer inputImage;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\VisibleWomanEyeSlice.png");
	try {
		reader->Update();
		inputImage = reader->GetOutput();
	} catch (itk::ExceptionObject& err) {
		std::cerr << "Exception caught!" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

	ImageType::Pointer outputImage = ImageType::New();
	outputImage->SetRegions(inputImage->GetRequestedRegion());
	outputImage->CopyInformation(inputImage);
	outputImage->Allocate();

	ConstIteratorType inputIt(inputImage, inputImage->GetRequestedRegion());
	IteratorType outputIt(outputImage, inputImage->GetRequestedRegion());
	inputIt.SetDirection(0);
	outputIt.SetDirection(0);

	for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); outputIt.NextLine(),inputIt.NextLine()) {
		inputIt.GoToBeginOfLine();
		outputIt.GoToEndOfLine();
		while (!inputIt.IsAtEndOfLine()) {
			--outputIt;
			outputIt.Set(inputIt.Get());
			++inputIt;
		}
	}
	WriteType::Pointer writer = WriteType::New();
	writer->SetFileName("G:/itkOutput/ImageData/out633.png");
	writer->SetInput(outputImage);
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void ImageLinearIteratorWithIndex633_2()
{
	typedef unsigned char PixelType;
	typedef itk::Image<PixelType, 3> Image3DType;
	typedef itk::Image<PixelType, 4> Image4DType;
	typedef itk::ImageFileReader<Image4DType> Reader4DType;
	typedef itk::ImageFileWriter<Image3DType> Writer3DType;

	Reader4DType::Pointer reader4D = Reader4DType::New();
	reader4D->SetFileName("");
	try {
		reader4D->Update();
	} catch (itk::ExceptionObject& err) {
		std::cerr << "Exception caught!" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	Image4DType::ConstPointer image4D = reader4D->GetOutput();

	Image3DType::Pointer image3D = Image3DType::New();
	typedef Image3DType::IndexType Index3DType;
	typedef Image3DType::SizeType Size3DType;
	typedef Image3DType::RegionType Region3DType;
	typedef Image3DType::SpacingType Spacing3DType;
	typedef Image3DType::PointType Origin3DType;
	typedef Image4DType::IndexType Index4DType;
	typedef Image4DType::SizeType Size4DType;
	typedef Image4DType::RegionType Region4DType;
	typedef Image4DType::SpacingType  Spacing4DType;
	typedef Image4DType::PointType Origin4DType;

	Index3DType index3D;
	Size3DType size3D;
	Region3DType region3D;
	Spacing3DType spacing3D;
	Origin3DType origin3D;
	Region4DType region4D = image4D->GetBufferedRegion();
	Index4DType index4D = region4D.GetIndex();
	Size4DType size4D = region4D.GetSize();
	Spacing4DType spacing4D = image4D->GetSpacing();
	Origin4DType origin4D = image4D->GetOrigin();

	for (unsigned int i = 0; i < 3; i ++) {
		size3D[i] = size4D[i];
		index3D[i] = index4D[i];
		spacing3D[i] = spacing4D[i];
		origin3D[i] = origin4D[i];
	}
	image3D->SetSpacing(spacing3D);
	image3D->SetOrigin(origin3D);
	region3D.SetIndex(index3D);
	region3D.SetSize(size3D);
	image3D->SetRegions(region3D);
	image3D->Allocate();

	typedef itk::NumericTraits<PixelType>::AccumulateType SumType;
	typedef itk::NumericTraits<SumType>::RealType MeanType;
	const unsigned int timeLength = region4D.GetSize()[3];

	typedef itk::ImageLinearConstIteratorWithIndex<Image4DType> IteratorType;
	IteratorType it(image4D, region4D);
	it.SetDirection(3);
	it.GoToBegin();
	while (!it.IsAtEnd()) {
		SumType sum = itk::NumericTraits<SumType>::ZeroValue();
		it.GoToBeginOfLine();
		index4D = it.GetIndex();
		while (!it.IsAtEndOfLine()) {
			sum += it.Get();
			++it;
		}
		MeanType mean = static_cast<MeanType>(sum);
		static_cast<MeanType>(timeLength);
		index3D[0] = index4D[0];
		index3D[1] = index4D[1];
		index3D[2] = index4D[2];
		image3D->SetPixel(index3D, static_cast<PixelType>(mean));
		it.NextLine();
	}

	Writer3DType::Pointer writer3D = Writer3DType::New();
	writer3D->SetFileName("G:/itkOutput/ImageData/out633_2.png");
	writer3D->SetInput(image3D);
	try {
		writer3D->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void ImageSliceIteratorWithIndex634()
{
	typedef unsigned short              PixelType;
	typedef itk::Image< PixelType, 2 >  ImageType2D;
	typedef itk::Image< PixelType, 3 >  ImageType3D;
	typedef itk::ImageLinearIteratorWithIndex< ImageType2D >
													LinearIteratorType;
	typedef itk::ImageSliceConstIteratorWithIndex< ImageType3D
													> SliceIteratorType;
	
	typedef itk::ImageFileReader< ImageType3D > ReaderType;
	typedef itk::ImageFileWriter< ImageType2D > WriterType;

	ImageType3D::ConstPointer inputImage;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainProtonDensity3Slices.mha");
	try {
		reader->Update();
		inputImage = reader->GetOutput();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	unsigned int projectionDirection =
		static_cast<unsigned int>(2);

	unsigned int i, j;
	unsigned int direction[2];
	for (i = 0, j = 0; i < 3; ++i) {
		if (i != projectionDirection) {
			direction[j] = i;
			j++;
		}
	}
	ImageType2D::RegionType region;
	ImageType2D::RegionType::SizeType size;
	ImageType2D::RegionType::IndexType index;

	ImageType3D::RegionType requestedRegion = inputImage->GetRequestedRegion();

	index[direction[0]] = requestedRegion.GetIndex()[direction[0]];
	index[1 - direction[0]] = requestedRegion.GetIndex()[direction[1]];
	size[direction[0]] = requestedRegion.GetSize()[direction[0]];
	size[1 - direction[0]] = requestedRegion.GetSize()[direction[1]];
	
	region.SetSize(size);
	region.SetIndex(index);

	ImageType2D::Pointer outputImage = ImageType2D::New();

	outputImage->SetRegions(region);
	outputImage->Allocate();

	SliceIteratorType  inputIt(inputImage, inputImage->GetRequestedRegion());
	LinearIteratorType outputIt(outputImage,
								outputImage->GetRequestedRegion());

	inputIt.SetFirstDirection(direction[1]);
	inputIt.SetSecondDirection(direction[0]);

	outputIt.SetDirection(1 - direction[0]);
	
	outputIt.GoToBegin();
	while (!outputIt.IsAtEnd()) {
		while (!outputIt.IsAtEndOfLine()) {
			outputIt.Set(itk::NumericTraits<unsigned short>::NonpositiveMin());
			++outputIt;
		}
		outputIt.NextLine();
	}

	inputIt.GoToBegin();
	outputIt.GoToBegin();

	while (!inputIt.IsAtEnd()) {
		while (!inputIt.IsAtEndOfSlice()) {
			while (!inputIt.IsAtEndOfLine()) {
				PixelType p0 = outputIt.Get();
				PixelType p1 = inputIt.Get();
				PixelType m = std::max(p0, p1);
				outputIt.Set(m);
				++inputIt;
				++outputIt;
			}
			outputIt.NextLine();
			inputIt.NextLine();

		}
		outputIt.GoToBegin();
		inputIt.NextSlice();
	}

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out634.png");
	writer->SetInput(outputImage);
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void ImageRandomConstIteratorWithIndex635()
{
	const unsigned int Dimension = 2;

	typedef unsigned short                                      PixelType;
	typedef itk::Image< PixelType, Dimension >                  ImageType;
	typedef itk::ImageRandomConstIteratorWithIndex<
											ImageType > ConstIteratorType;

	typedef itk::ImageFileReader< ImageType > ReaderType;

	ImageType::ConstPointer inputImage;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\RatLungSlice1.mha");
	try {
		reader->Update();
		inputImage = reader->GetOutput();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return ;
	}
	const int sampleCount = 100;
	ConstIteratorType inputIt(inputImage, inputImage->GetRequestedRegion());
	inputIt.SetNumberOfSamples(sampleCount);
	inputIt.ReinitializeSeed();

	float mean = 0.0f;
	for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt) {
		mean += static_cast<float>(inputIt.Get());
	}
	mean = mean / (double)sampleCount;
	std::cout << "Mean estimate with " << sampleCount << " samples is " << mean << std::endl;
}
void NeighborhoodIterator641_1()
{
	typedef float                             PixelType;
	typedef itk::Image< PixelType, 2 >        ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ImageType>        IteratorType;
	
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainT1Slice.png");
	try {
		reader->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

	NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);
	NeighborhoodIteratorType it(radius, reader->GetOutput(),
								reader->GetOutput()->GetRequestedRegion());
	
	ImageType::Pointer output = ImageType::New();
	output->SetRegions(reader->GetOutput()->GetRequestedRegion());
	output->Allocate();

	IteratorType out(output, reader->GetOutput()->GetRequestedRegion());
	// x
// 	NeighborhoodIteratorType::OffsetType offset1 = { { -1, -1 } };
// 	NeighborhoodIteratorType::OffsetType offset2 = { { 1, -1 } };
// 	NeighborhoodIteratorType::OffsetType offset3 = { { -1, 0 } };
// 	NeighborhoodIteratorType::OffsetType offset4 = { { 1, 0 } };
// 	NeighborhoodIteratorType::OffsetType offset5 = { { -1, 1 } };
// 	NeighborhoodIteratorType::OffsetType offset6 = { { 1, 1 } };
	//y
	NeighborhoodIteratorType::OffsetType offset1 = { { -1, -1 } };
	NeighborhoodIteratorType::OffsetType offset2 = { { -1, 1 } };
	NeighborhoodIteratorType::OffsetType offset3 = { { 0, -1 } };
	NeighborhoodIteratorType::OffsetType offset4 = { { 0, 1 } };
	NeighborhoodIteratorType::OffsetType offset5 = { { 1, -1 } };
	NeighborhoodIteratorType::OffsetType offset6 = { { 1, 1 } };

	for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out) {
		float sum;
		sum = it.GetPixel(offset2) - it.GetPixel(offset1);
		sum += 2.0 * it.GetPixel(offset4) - 2.0 * it.GetPixel(offset3);
		sum += it.GetPixel(offset6) - it.GetPixel(offset5);
		out.Set(sum);
	}
	typedef unsigned char                          WritePixelType;
	typedef itk::Image< WritePixelType, 2 >        WriteImageType;
	typedef itk::ImageFileWriter< WriteImageType > WriterType;

	typedef itk::RescaleIntensityImageFilter<
		ImageType, WriteImageType > RescaleFilterType;

	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(output);

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out6411.png");
	writer->SetInput(rescaler->GetOutput());
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void NeighborhoodIterator641_2()
{
	typedef float                             PixelType;
	typedef itk::Image< PixelType, 2 >        ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ImageType>        IteratorType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainT1Slice.png");
	try {
		reader->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	ImageType::Pointer output = ImageType::New();
	output->SetRegions(reader->GetOutput()->GetRequestedRegion());
	output->Allocate();

	IteratorType out(output, reader->GetOutput()->GetRequestedRegion());

	const unsigned int direction = 1;
	itk::SobelOperator<PixelType, 2> sobelOperator;
	sobelOperator.SetDirection(direction);
	sobelOperator.CreateDirectional();

	NeighborhoodIteratorType::RadiusType radius = sobelOperator.GetRadius();
	NeighborhoodIteratorType it(radius, reader->GetOutput(),
								reader->GetOutput()->GetRequestedRegion());
	itk::NeighborhoodInnerProduct<ImageType> innerProduct;
	for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++ out) {
		out.Set(innerProduct(it, sobelOperator));
	}
	typedef unsigned char                          WritePixelType;
	typedef itk::Image< WritePixelType, 2 >        WriteImageType;
	typedef itk::ImageFileWriter< WriteImageType > WriterType;

	typedef itk::RescaleIntensityImageFilter<
		ImageType, WriteImageType > RescaleFilterType;

	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(output);

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out6412.png");
	writer->SetInput(rescaler->GetOutput());
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void NeighborhoodIterator641_3()
{
	typedef float                             PixelType;
	typedef itk::Image< PixelType, 2 >        ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ImageType>        IteratorType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainT1Slice.png");
	try {
		reader->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	ImageType::Pointer output = ImageType::New();
	output->SetRegions(reader->GetOutput()->GetRequestedRegion());
	output->Allocate();

	const unsigned int direction = 1;
	itk::SobelOperator<PixelType, 2> sobelOperator;
	sobelOperator.SetDirection(direction);
	sobelOperator.CreateDirectional();

	itk::NeighborhoodInnerProduct<ImageType> innerProduct;
	typedef itk::NeighborhoodAlgorithm
		::ImageBoundaryFacesCalculator<ImageType> FaceCalculatorType;
	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;

	faceList = faceCalculator(reader->GetOutput(), output->GetRequestedRegion(),
							  sobelOperator.GetRadius());
	FaceCalculatorType::FaceListType::iterator fit;
	IteratorType out;
	NeighborhoodIteratorType it;
	for (fit = faceList.begin(); fit != faceList.end(); ++ fit) {
		it = NeighborhoodIteratorType(sobelOperator.GetRadius(),
									  reader->GetOutput(), *fit);
		out = IteratorType(output, *fit);
		for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++ it, ++ out) {
			out.Set(innerProduct(it, sobelOperator));
		}
	}


	typedef unsigned char                          WritePixelType;
	typedef itk::Image< WritePixelType, 2 >        WriteImageType;
	typedef itk::ImageFileWriter< WriteImageType > WriterType;

	typedef itk::RescaleIntensityImageFilter<
		ImageType, WriteImageType > RescaleFilterType;

	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(output);

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out6413.png");
	writer->SetInput(rescaler->GetOutput());
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void NeighborhoodIterator641_4()
{
	typedef float                             PixelType;
	typedef itk::Image< PixelType, 2 >        ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ImageType>        IteratorType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainT1Slice.png");
	try {
		reader->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
	ImageType::Pointer output = ImageType::New();
	output->SetRegions(reader->GetOutput()->GetRequestedRegion());
	output->Allocate();

	const float standardDiff = 0.9;
	itk::GaussianOperator<PixelType, 2> gaussianOperator;
	gaussianOperator.SetVariance(standardDiff*standardDiff);

	itk::NeighborhoodInnerProduct<ImageType> innerProduct;
	typedef itk::NeighborhoodAlgorithm
		::ImageBoundaryFacesCalculator<ImageType> FaceCalculatorType;
	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;

	faceList = faceCalculator(reader->GetOutput(), output->GetRequestedRegion(),
							  gaussianOperator.GetRadius());
	FaceCalculatorType::FaceListType::iterator fit;
	IteratorType out;
	NeighborhoodIteratorType it;
	ImageType::Pointer input = reader->GetOutput();
	for (unsigned int i = 0; i < ImageType::ImageDimension; ++i) {
		gaussianOperator.SetDirection(i);
		gaussianOperator.CreateDirectional();

		faceList = faceCalculator(input, output->GetRequestedRegion(),
								  gaussianOperator.GetRadius());
		for (fit = faceList.begin(); fit != faceList.end(); ++fit) {
			it = NeighborhoodIteratorType(gaussianOperator.GetRadius(),
										  input, *fit);
			out = IteratorType(output, *fit);
			for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out) {
				out.Set(innerProduct(it, gaussianOperator));
			}
		}

		if (i != ImageType::ImageDimension - 1) {
			ImageType::Pointer tmp = input;
			input = output;
			output = tmp;
		}
	}


	typedef unsigned char                          WritePixelType;
	typedef itk::Image< WritePixelType, 2 >        WriteImageType;
	typedef itk::ImageFileWriter< WriteImageType > WriterType;

	typedef itk::RescaleIntensityImageFilter<
		ImageType, WriteImageType > RescaleFilterType;

	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(output);

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out6414.png");
	writer->SetInput(rescaler->GetOutput());
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void NeighborhoodIterator641_5()
{
	typedef float                             PixelType;
	typedef itk::Image< PixelType, 2 >        ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ImageType>        IteratorType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainT1Slice.png");
	try {
		reader->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

	ImageType::Pointer output = ImageType::New();
	output->SetRegions(reader->GetOutput()->GetRequestedRegion());
	output->Allocate();

	itk::NeighborhoodInnerProduct<ImageType> innerProduct;

	typedef itk::NeighborhoodAlgorithm
		::ImageBoundaryFacesCalculator< ImageType > FaceCalculatorType;

	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;
	FaceCalculatorType::FaceListType::iterator fit;

	IteratorType out;
	NeighborhoodIteratorType it;

	const float standardDiff = 0.9;
	itk::GaussianOperator< PixelType, 2 > gaussianOperator;
	gaussianOperator.SetDirection(0);
	gaussianOperator.SetVariance(standardDiff * standardDiff);
	gaussianOperator.CreateDirectional();

	NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(gaussianOperator.GetRadius()[0]);

	ImageType::Pointer input = reader->GetOutput();
	faceList = faceCalculator(input, output->GetRequestedRegion(), radius);

	for (unsigned int i = 0; i < ImageType::ImageDimension; ++i) {
		for (fit = faceList.begin(); fit != faceList.end(); ++fit) {
			it = NeighborhoodIteratorType(radius, input, *fit);
			out = IteratorType(output, *fit);
			for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out) {
				out.Set(innerProduct(it.GetSlice(i), it, gaussianOperator));
			}
		}

		// Swap the input and output buffers
		if (i != ImageType::ImageDimension - 1) {
			ImageType::Pointer tmp = input;
			input = output;
			output = tmp;
		}
	}

	typedef unsigned char                          WritePixelType;
	typedef itk::Image< WritePixelType, 2 >        WriteImageType;
	typedef itk::ImageFileWriter< WriteImageType > WriterType;

	typedef itk::RescaleIntensityImageFilter< ImageType,
		WriteImageType > RescaleFilterType;

	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(output);

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out6415.png");
	writer->SetInput(rescaler->GetOutput());
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void NeighborhoodIterator641_6()
{
	typedef float                             PixelType;
	typedef itk::Image< PixelType, 2 >        ImageType;
	typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
	typedef itk::ImageFileReader< ImageType > ReaderType;

	typedef itk::FastMarchingImageFilter< ImageType, ImageType> FastMarchingFilterType;

	FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();

	typedef FastMarchingFilterType::NodeContainer NodeContainer;
	typedef FastMarchingFilterType::NodeType NodeType;

	NodeContainer::Pointer seeds = NodeContainer::New();

	ImageType::IndexType seedPosition;
	seedPosition[0] = 128;
	seedPosition[1] = 128;
	const double initialDistance = 1.0;

	NodeType node;

	const double seedValue = -initialDistance;
	ImageType::SizeType size = { { 256, 256 } };

	node.SetValue(seedValue);
	node.SetIndex(seedPosition);
	seeds->Initialize();
	seeds->InsertElement(0, node);

	fastMarching->SetTrialPoints(seeds);
	fastMarching->SetSpeedConstant(1.0);

	itk::AddImageFilter<ImageType, ImageType, ImageType>::Pointer adder
		= itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
	itk::RandomImageSource<ImageType>::Pointer noise
		= itk::RandomImageSource<ImageType>::New();

	noise->SetSize(size.m_Size);
	noise->SetMin(-.7);
	noise->SetMax(.8);
	adder->SetInput1(noise->GetOutput());
	adder->SetInput2(fastMarching->GetOutput());

	try {
		fastMarching->SetOutputSize(size);
		fastMarching->Update();

		adder->Update();

	} catch (itk::ExceptionObject & excep) {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
	}

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainT1Slice.png");
	try {
		reader->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

	ImageType::Pointer input = adder->GetOutput();
	ImageType::IndexType index;
	index[0] = 10;
	index[1] = 10;
	NeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);
	NeighborhoodIteratorType it(radius, input, input->GetRequestedRegion());

	it.SetLocation(index);

	bool flag = true;
	while (flag == true) {
		NeighborhoodIteratorType::OffsetType nextMove;
		nextMove.Fill(0);

		flag = false;

		PixelType min = it.GetCenterPixel();
		for (unsigned i = 0; i < it.Size(); i++) {
			if (it.GetPixel(i) < min) {
				min = it.GetPixel(i);
				nextMove = it.GetOffset(i);
				flag = true;
			}
		}
		it.SetCenterPixel(255.0);
		it += nextMove;
	}

	typedef unsigned char                          WritePixelType;
	typedef itk::Image< WritePixelType, 2 >        WriteImageType;
	typedef itk::ImageFileWriter< WriteImageType > WriterType;

	typedef itk::RescaleIntensityImageFilter< ImageType,
		WriteImageType > RescaleFilterType;

	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(input);

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out6416.png");
	writer->SetInput(rescaler->GetOutput());
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}
}
void NeighborhoodIterator642()
{
	typedef unsigned char PixelType;
	typedef itk::Image<PixelType, 2> ImageType;
	typedef itk::ConstShapedNeighborhoodIterator<ImageType> ShapedNeighborhoodIteratorType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	typedef itk::ImageFileReader< ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BinaryImage.png");
	try {
		reader->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}

	ImageType::Pointer output = ImageType::New();
	output->SetRegions(reader->GetOutput()->GetRequestedRegion());
	output->Allocate();

	unsigned int element_radius = 1;
	ShapedNeighborhoodIteratorType::RadiusType radius;
	radius.Fill(element_radius);

	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<
		ImageType> FaceCalculatorType;

	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;
	FaceCalculatorType::FaceListType::iterator fit;
	faceList = faceCalculator(reader->GetOutput(),
							  output->GetRequestedRegion(),
							  radius);

	IteratorType out;

	const PixelType background_value = 0;
	const PixelType foreground_value = 255;
	const float rad = static_cast<float>(element_radius);

	for (fit = faceList.begin(); fit != faceList.end(); ++fit) {
		ShapedNeighborhoodIteratorType it(radius, reader->GetOutput(), *fit);
		out = IteratorType(output, *fit);
		for (float y = -rad; y <= rad; y ++) {
			for (float x = -rad; x <= rad; x++) {
				ShapedNeighborhoodIteratorType::OffsetType off;
				float dis = std::sqrt(x*x + y*y);
				if (dis <= rad) {
					off[0] = static_cast<int>(x);
					off[1] = static_cast<int>(y);
					it.ActivateOffset(off);
				}
			}
		}
		for (it.GoToBegin(), out.GoToBegin(); !it.IsAtEnd(); ++it, ++out) {
			ShapedNeighborhoodIteratorType::ConstIterator ci;
			bool flag = true;

			for (ci = it.Begin(); ci != it.End(); ci ++) {
				if (ci.Get() == background_value) {
					flag = false;
					break;
				}
			}
			if (true == flag) {
				out.Set(foreground_value);
			} else {
				out.Set(background_value);
			}
		}
	}
	typedef itk::ImageFileWriter< ImageType > WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("G:\\itkOutput\\ImageData\\out642.png");
	writer->SetInput(output);
	try {
		writer->Update();
	} catch (itk::ExceptionObject &err) {
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return;
	}




}






