/* main.cpp:  Bed Removal from Clinical CT Scans using ITK filters.
//////////////////////////////////////////////////////////////////////
BED REMOVAL FROM CT SCANS  v1.0.0 (12 Aug 2013)
//////////////////////////////////////////////////////////////////////
Center for Infectious Disease Imaging (CIDI),
Radiology and Imaging Sciences,
National Institutes of Health (NIH),
Bethesda MD
12 Aug 2013
Written by Awais Mansoor, PhD
Primary Contact: awais.mansoor@nih.gov
Secondary Contact: Ulas.bagci@nih.gov
and
Daniel.mollura@nih.gov

//////////////////////////////////////////////////////////////////////
This code removes the bed out of the abdominal CT scan using intensity
thresholding coupled with morphological operations.
*/

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h" 


typedef float PixelType;
typedef unsigned char OutPixelType;

typedef itk::Image<PixelType, 3>  ImageType;
typedef itk::Image< OutPixelType, 3 > OutImageType;
 

int main(int argc, char * argv[])
{
	if( argc < 3 )
    {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0] << " inputImageFile outputImageFile Output Type (0=Bed Removed Image, 1=Body Mask; def.=0)" << std::endl;
      return EXIT_FAILURE;
    } 


  bool Output=0;

  if (argc==4) Output=argv[3];


  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();  
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  
  writer->SetFileName( argv[2] );




  
  //Pipeline
  try
    {
      reader->Update();
	  
    }
  catch ( itk::ExceptionObject &err)
    {
      std::cout<<"Problems reading input image"<<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }

  
  //Get image specs
  ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing(); 
  ImageType::PointType origin = reader->GetOutput()->GetOrigin(); 
  ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
  ImageType::SizeType  insize = reader->GetOutput()->GetRequestedRegion().GetSize();
  int pRow, pCol, pSli;
  pRow = insize[0];
  pCol = insize[1];
  pSli = insize[2]; 
  ImageType::RegionType region;
  region.SetSize( insize );

  
  //SKin BoundaryExtraction Starts

  ImageType::SizeType outputSize;
  outputSize[0]=insize[0]/3;
  outputSize[1]=insize[1]/3;
  outputSize[2]=insize[2];
  ImageType::SpacingType outputSpacing;
  outputSpacing[0] = reader->GetOutput()->GetSpacing()[0] * (static_cast<double>(insize[0]) / outputSize[0]);
  outputSpacing[1] = reader->GetOutput()->GetSpacing()[1] * (static_cast<double>(insize[1]) / outputSize[1]);
  outputSpacing[2] = reader->GetOutput()->GetSpacing()[2];

  typedef itk::IdentityTransform<double, 3> TransformType; 
  
  
  //std::cout<<"Entering Skin-body Edxtraction..."<<std::endl;

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resampleDN = ResampleImageFilterType::New();
  resampleDN->SetInput(reader->GetOutput());
  resampleDN->SetSize(outputSize);
  resampleDN->SetOutputSpacing(outputSpacing);
  resampleDN->SetOutputOrigin( origin );  
  resampleDN->SetOutputDirection( direction );
  resampleDN->SetTransform(TransformType::New());
  //resampleDN->SetInterpolator(pInterpolator);
  resampleDN->UpdateLargestPossibleRegion();

  
  typedef itk::BinaryThresholdImageFilter<ImageType, OutImageType > BinaryThresholdImageFilterType;
  BinaryThresholdImageFilterType::Pointer binaryskinfilter = BinaryThresholdImageFilterType::New();
  binaryskinfilter->SetInput( resampleDN->GetOutput() );
  binaryskinfilter->SetOutsideValue( 0 );
  binaryskinfilter->SetInsideValue( 1 );
  binaryskinfilter->SetLowerThreshold( -300 );
  binaryskinfilter->SetUpperThreshold( 3071 );
  binaryskinfilter->Update();
  
  typedef itk::BinaryBallStructuringElement<OutImageType::PixelType, OutImageType::ImageDimension>
              StructuringElementSkinType;
  StructuringElementSkinType structuringElementSkin;
  structuringElementSkin.SetRadius(3);
  structuringElementSkin.CreateStructuringElement();

  typedef itk::BinaryMorphologicalOpeningImageFilter <OutImageType, OutImageType, StructuringElementSkinType>  BinaryMorphologicalOpeningImageFilter;
  BinaryMorphologicalOpeningImageFilter::Pointer openingFilter
          = BinaryMorphologicalOpeningImageFilter::New();
  openingFilter->SetInput(binaryskinfilter->GetOutput());
  openingFilter->SetKernel(structuringElementSkin);
  openingFilter->SetForegroundValue(1);
  openingFilter->Update();


  
  typedef itk::BinaryMorphologicalClosingImageFilter <OutImageType, ImageType, StructuringElementSkinType>
          BinaryMorphologicalClosingImageFilter;

  BinaryMorphologicalClosingImageFilter::Pointer closingFilter = BinaryMorphologicalClosingImageFilter::New();
  closingFilter->SetInput(openingFilter->GetOutput());
  structuringElementSkin.SetRadius(40);
  structuringElementSkin.CreateStructuringElement();
  closingFilter->SetKernel(structuringElementSkin);
  closingFilter->SetForegroundValue(1);
  closingFilter->Update();


  
  ResampleImageFilterType::Pointer resampleUP = ResampleImageFilterType::New();
  resampleUP->SetInput(closingFilter->GetOutput());
  resampleUP->SetSize(insize);
  resampleUP->SetOutputSpacing(reader->GetOutput()->GetSpacing());
  resampleUP->SetOutputOrigin( origin );  
  resampleUP->SetOutputDirection( direction );
  resampleUP->SetTransform(TransformType::New());
  //resampleUP->SetInterpolator(pInterpolator);
  resampleUP->UpdateLargestPossibleRegion();
  //Skin Boundary Extraction Ends

  //std::cout<<"Ending Skin-body Edxtraction..."<<std::endl;
  
  
  ImageType::IndexType pixelIndex;

  for (int k=0; k<pSli; k++)
	for(int i=0; i<pRow; i++)
		for(int j=0; j<pCol; j++)
		{
			pixelIndex[0]=i;
			pixelIndex[1]=j;
			pixelIndex[2]=k;

			if(resampleUP->GetOutput()->GetPixel(pixelIndex)==0 && resampleUP->GetOutput()->GetPixel(pixelIndex)>-3000)
				reader->GetOutput()->SetPixel(pixelIndex,-1024);
			
		
		}
 
  

	reader->Update();

	if(Output==1)
	 writer->SetInput(resampleUP->GetOutput());
	else
	 writer->SetInput(reader->GetOutput());

  try
    {
      writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
      std::cout<<"ExceptionObject caught !"<<std::endl;
      std::cout<< err <<std::endl;
      return EXIT_FAILURE;
    }
 
  return EXIT_SUCCESS;
}
 

