#include "Filtering\ImageGradient\include\itkGradientMagnitudeImageFilter.h"
#include "Filtering\Thresholding\include\itkThresholdImageFilter.h"
#include "Filtering\ImageIntensity\include\itkRescaleIntensityImageFilter.h"
#include "Common\include\itkImageToImageFilter.h"
#include "ImageBase\include\itkImageFileReader.h"
#include <ImageBase\include\itkImageFileWriter.h>

namespace itk
{
	template <typename TImage>
	class CompositeExampleImageFilter : public ImageToImageFilter<
		TImage, TImage>
	{
	public:
		typedef CompositeExampleImageFilter Self;
		typedef ImageToImageFilter<TImage, TImage> Superclass;
		typedef SmartPointer<Self> Pointer;
		typedef SmartPointer<const Self> ConstPointer;
		
		itkNewMacro(Self);
		itkTypeMacro(CompositeExampleImageFilter, ImageToImageFilter);
		
		typedef TImage                        ImageType;
		typedef typename ImageType::PixelType PixelType;

		itkGetMacro(Threshold, PixelType);
		itkSetMacro(Threshold, PixelType);
	protected:

		CompositeExampleImageFilter();
	protected:

		typedef ThresholdImageFilter< ImageType >                    ThresholdType;
		typedef GradientMagnitudeImageFilter< ImageType, ImageType > GradientType;
		typedef RescaleIntensityImageFilter< ImageType, ImageType >  RescalerType;

		virtual void GenerateData() ITK_OVERRIDE;

		/** Display */
		void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
	private:
		ITK_DISALLOW_COPY_AND_ASSIGN(CompositeExampleImageFilter);
		typename GradientType::Pointer     m_GradientFilter;
		typename ThresholdType::Pointer    m_ThresholdFilter;
		typename RescalerType::Pointer     m_RescaleFilter;

		PixelType m_Threshold;
	};
	template< typename TImage>
	CompositeExampleImageFilter< TImage> ::CompositeExampleImageFilter()
	{
		m_Threshold = 1;
		m_GradientFilter = GradientType::New();
		m_ThresholdFilter = ThresholdType::New();
		m_ThresholdFilter->SetInput(m_GradientFilter->GetOutput());
		m_RescaleFilter = RescalerType::New();
		m_RescaleFilter->SetInput(m_ThresholdFilter->GetOutput());
		m_RescaleFilter->SetOutputMinimum(
			NumericTraits<PixelType>::NonpositiveMin());
		m_RescaleFilter->SetOutputMaximum(NumericTraits<PixelType>::max());
	}
	template< typename TImage >
	void
		CompositeExampleImageFilter< TImage >
		::GenerateData()
	{
		typename ImageType::Pointer input = ImageType::New();
		input->Graft(const_cast<ImageType*>(this->GetInput()));
		m_GradientFilter->SetInput(input);

		m_ThresholdFilter->ThresholdBelow(this->m_Threshold);

		m_RescaleFilter->GraftOutput(this->GetOutput());
		m_RescaleFilter->Update();
		this->GraftOutput(m_RescaleFilter->GetOutput());
	}
	template< typename TImage >
	void
		CompositeExampleImageFilter< TImage >
		::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);

		os << indent << "Threshold:" << this->m_Threshold
			<< std::endl;
	}


}

void CompositeExample8()
{
	typedef itk::Image< short, 2 >            ImageType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName("G:\\InsightToolkit-4.13.0\\Examples\\Data\\BrainProtonDensitySlice.png");

	typedef itk::CompositeExampleImageFilter<ImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(reader->GetOutput());
	filter->SetThreshold(20);

	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(filter->GetOutput());
	writer->SetFileName("G:\\itkOutput\\ImageData\\out85.png");

	try {
		writer->Update();
	} catch (itk::ExceptionObject & e) {
		std::cerr << "Error: " << e << std::endl;
	}
}




