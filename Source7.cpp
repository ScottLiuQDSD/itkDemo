
#include "ImageAdaptors\include\itkImageAdaptor.h"
#include "Common\include\itkImageRegionIteratorWithIndex.h"
#include "ImageBase\include\itkImageFileReader.h"
#include "ImageBase\include\itkImageFileWriter.h"
#include "Filtering\ImageIntensity\include\itkRescaleIntensityImageFilter.h"
#include "Filtering\ImageGradient\include\itkGradientRecursiveGaussianImageFilter.h"

class CastPixelAccessor
{
public:
	typedef unsigned char		InternalType;
	typedef float				ExternalType;
	static void Set(InternalType & output, const ExternalType & input)
	{
		output = static_cast<InternalType>(input);
	}

	static ExternalType Get(const InternalType & input)
	{
		return static_cast<ExternalType>(input);
	}

};
void ImageAdaptor7_1_0()
{
	typedef unsigned char InputPiexlType;
	const unsigned int Dimension = 2;
	typedef itk::Image<InputPiexlType, Dimension> ImageType;

	typedef itk::ImageAdaptor<ImageType, CastPixelAccessor> ImageAdaptorType;
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BinaryImage.png");
	reader->Update();

	adaptor->SetImage(reader->GetOutput());
	typedef itk::ImageRegionIteratorWithIndex<ImageAdaptorType> IteratorType;
	IteratorType it(adaptor, adaptor->GetBufferedRegion());
	
	double sum = 0.0;
	it.GoToBegin();
	while (!it.IsAtEnd()) {
		float value = it.Get();
		sum += value;
		++it;
	}
	std::cout << "Sum of pixels is: " << sum << std::endl;


}

class RedChannelPixelAccessor
{
public:
	typedef itk::RGBPixel<float> InternalType;
	typedef float ExternalType;
	static ExternalType Get(const InternalType &input)
	{
		return static_cast<ExternalType>(input.GetRed());
	}
};
class BlueChannelPixelAccessor {
public:
	typedef itk::RGBPixel<float> InternalType;
	typedef float ExternalType;
	static ExternalType Get(const InternalType &input)
	{
		return static_cast<ExternalType>(input.GetBlue());
	}
};
void ImageAdaptor7_1_1()
{
// 	typedef RedChannelPixelAccessor::InternalType  InputPixelType;
	typedef BlueChannelPixelAccessor::InternalType  InputPixelType;
	const   unsigned int   Dimension = 2;
	typedef itk::Image< InputPixelType, Dimension >   ImageType;

// 	typedef itk::ImageAdaptor<  ImageType,
// 		RedChannelPixelAccessor > ImageAdaptorType;
	typedef itk::ImageAdaptor<  ImageType,
		BlueChannelPixelAccessor > ImageAdaptorType;
	
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();

	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\VisibleWomanHeadSlice.png");

	reader->Update();

	adaptor->SetImage(reader->GetOutput());

	typedef itk::Image<unsigned char, Dimension> OutputImageType;
	typedef itk::RescaleIntensityImageFilter< ImageAdaptorType,
		OutputImageType > RescalerType;

	RescalerType::Pointer rescaler = RescalerType::New();
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName("G:\\itkOutput\\ImageData\\out711.png");

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(adaptor);
	writer->SetInput(rescaler->GetOutput());
	try {
		writer->Update();
	} catch (itk::ExceptionObject & excp) {
		std::cerr << "Exception caught " << excp << std::endl;
		return;
	}

}

class VectorPixelAccessor {
public:
	typedef itk::CovariantVector<float, 2>   InternalType;
	typedef                      float      ExternalType;

	VectorPixelAccessor() : m_Index(0) {}

	VectorPixelAccessor & operator=(const VectorPixelAccessor & vpa)
	{
		m_Index = vpa.m_Index;
		return *this;
	}
	ExternalType Get(const InternalType & input) const
	{
		return static_cast<ExternalType>(input[m_Index]);
	}
	void SetIndex(unsigned int index)
	{
		m_Index = index;
	}

private:
	unsigned int m_Index;
};
void VectorPixelAccessor7_3()
{
	typedef unsigned char  InputPixelType;
	const   unsigned int   Dimension = 2;
	typedef itk::Image< InputPixelType, Dimension  >   InputImageType;
	typedef itk::CovariantVector< float, Dimension  >   VectorPixelType;
	typedef itk::Image< VectorPixelType, Dimension  >   VectorImageType;
	typedef itk::GradientRecursiveGaussianImageFilter< InputImageType,
		VectorImageType> GradientFilterType;

	GradientFilterType::Pointer gradient = GradientFilterType::New();
	typedef itk::ImageAdaptor<  VectorImageType,
		VectorPixelAccessor > ImageAdaptorType;

	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	VectorPixelAccessor  accessor;
	accessor.SetIndex(1);
	adaptor->SetPixelAccessor(accessor);

	typedef itk::ImageFileReader< InputImageType >   ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	gradient->SetInput(reader->GetOutput());

	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainProtonDensitySlice.png");
	gradient->Update();

	adaptor->SetImage(gradient->GetOutput());

	typedef itk::Image< unsigned char, Dimension >   OutputImageType;
	typedef itk::RescaleIntensityImageFilter< ImageAdaptorType, OutputImageType>
		RescalerType;
	RescalerType::Pointer rescaler = RescalerType::New();
	typedef itk::ImageFileWriter< OutputImageType >   WriterType;
	WriterType::Pointer writer = WriterType::New();

	writer->SetFileName("G:\\itkOutput\\ImageData\\out73.png");

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);

	rescaler->SetInput(adaptor);
	writer->SetInput(rescaler->GetOutput());
	writer->Update();

}

class ThresholdingPixelAccessor
{
public:
	typedef unsigned char InternalType;
	typedef unsigned char ExternalType;
	ThresholdingPixelAccessor() : m_Threshold(0) {};
	ExternalType Get(const InternalType & input) const
	{
		return (input > m_Threshold) ? 1 : 0;
	}
	void SetThreshold(const InternalType threshold) 
	{
		m_Threshold = threshold;
	}
	ThresholdingPixelAccessor &
		operator=(const ThresholdingPixelAccessor & vpa)
	{
		m_Threshold = vpa.m_Threshold;
		return *this;
	}

private:
	InternalType m_Threshold;
};

void ThresholdingAdaptor7_4()
{
	typedef ThresholdingPixelAccessor::InternalType     PixelType;
	const   unsigned int   Dimension = 2;
	typedef itk::Image< PixelType, Dimension >   ImageType;

	typedef itk::ImageAdaptor< ImageType,
		ThresholdingPixelAccessor > ImageAdaptorType;
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	ThresholdingPixelAccessor accessor;
	const int threshold = 880;
	accessor.SetThreshold(threshold);
	adaptor->SetPixelAccessor(accessor);

	typedef itk::ImageFileReader< ImageType >   ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainProtonDensitySlice.png");
	reader->Update();

	adaptor->SetImage(reader->GetOutput());

	typedef itk::RescaleIntensityImageFilter< ImageAdaptorType,
		ImageType > RescalerType;

	RescalerType::Pointer rescaler = RescalerType::New();
	typedef itk::ImageFileWriter< ImageType >   WriterType;
	WriterType::Pointer writer = WriterType::New();


	writer->SetFileName("G:\\itkOutput\\ImageData\\out74.png");

	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);

	rescaler->SetInput(adaptor);
	writer->SetInput(rescaler->GetOutput());
	writer->Update();

}



