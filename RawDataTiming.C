#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <memory>

R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>

namespace
{
  // small class to easily save multple pages pdf

  class PdfDocument
  {

    public:

    //* constructor
    PdfDocument( TString filename ): m_filename( filename )
    {}

    //* destructor
    ~PdfDocument()
    {
      if( m_first || m_filename.IsNull() ) return;
      // we save a new empty canvas at the end of the pdf file so that it is properly closed
      TCanvas().SaveAs( Form( "%s)", m_filename.Data() ) );
    }

    //* add pad
    void Add( TVirtualPad* pad )
    {
      if( m_filename.IsNull() ) return;
      if( m_first )
      {
        pad->SaveAs( Form( "%s(", m_filename.Data() ) );
        m_first = false;
      } else {
        pad->SaveAs( Form( "%s", m_filename.Data() ) );
      }
    }

    private:

    //* filename
    TString m_filename;

    //* true if first page
    Bool_t m_first = true;

    //* root dictionary
    ClassDef( PdfDocument, 0 );

  };

}

//_____________________________________________________________________________
TString RawDataTiming( int runNumber = 29863 )
{

  MicromegasMapping mapping;

  const TString inputFile = Form( "MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  const TString pdfFile = Form( "RawDataTiming-%08i-0000.pdf", runNumber );

  std::cout << "RawDataTiming - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataTiming - pdfFile: " << pdfFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  // open TFile
  auto f = TFile::Open( inputFile );
  if( !f )
  {
    std::cout << "RawDataTiming - could not open " << inputFile << std::endl;
    return TString();
  }

  // get tree
  auto tree = dynamic_cast<TTree*>(f->Get("T"));

  // project adc vs sample vs strip in a 3D histogram
  // strips are grouped by chuncks of 64 corresponding to the 4 Resist region in each of the detectors
  const TString var( "samples.adc:samples.sample" );
  const TString var3d = Form( "%s:((samples.strip/64)+4*(samples.tile+8*(samples.layer-55)))", var.Data() );
  auto h3d = new TH3I( "h3d", "h3d", 64, 0, 64, 150, 0, 150, 1200, 0, 1200 );
  h3d->GetXaxis()->SetTitle( "region_id" );
  h3d->GetYaxis()->SetTitle( "sample" );
  h3d->GetZaxis()->SetTitle( "adc" );

  // project
  tree->Project( h3d->GetName(), var3d, TCut() );
  std::cout << "projection done." << std::endl;

  // loop over layers, tiles and region
  for( int ilayer = 0; ilayer <2; ++ilayer )
  {

    // get actual layer and segmentation
    const int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;

    for( int itile = 0; itile <8; ++itile )
    {
      // generate hitset key from layer, tile and segmentation. This allows to retrieve the detector name from the Mapping
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );

      // create canvas and divide
      const auto cvname = Form( "cv_%i_%i", ilayer, itile );
      auto cv = new TCanvas( cvname, cvname, 900, 900 );
      cv->Divide(2,2);

      for( int iregion = 0; iregion < 4; ++iregion )
      {

        // due to mapping internals, region are named opposite to strips.
        const int region = 4-iregion;
        const auto hname = Form( "h_%i_%i_%i", ilayer, itile, iregion );
        const auto htitle = Form( "%s_R%i - %i,%i", name.c_str(), region, layer, itile );

        // get the bin matching layer, tile and region in the 3D histgoram
        const int bin = iregion + 4*(itile + 8*ilayer );
        h3d->GetXaxis()->SetRange( bin+1, bin+1 );

        // get the corresponding 2D histogram (adc vs sample)
        auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
        h2d->SetName( hname );
        h2d->SetTitle( htitle );

        // draw
        cv->cd(iregion+1);
        h2d->Draw("colz");
      }

      // add to pdfdocument
      pdfDocument.Add( cv );

    }
  }

  return pdfFile;

}
