
#include <TH1.h>
#include <TF1.h>

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
    PdfDocument( TString filename ): m_filename( filename ) {}

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
  };

}

//_____________________________________________________________________________
// TString RawDataSignal(int runNumber = 29863)
TString RawDataSignal(int runNumber = 25475)
{

  MicromegasMapping mapping;

  const TString inputFile = Form( "MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  const TString pdfFile = Form( "RawDataSignal-%08i-0000.pdf", runNumber );

  std::cout << "RawDataSignal - inputFile: " << inputFile << std::endl;
  std::cout << "RawDataSignal - pdfFile: " << pdfFile << std::endl;

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

  const TString var( "samples.adc:samples.strip" );
  const TString var3d = Form( "%s:(samples.tile+8*(samples.layer-55))", var.Data() );
  // auto h3d = new TH3I( "h3d", "h3d", 16, 0, 16, 256, 0, 256, 1024, 0, 1024 );
  auto h3d = new TH3F( "h3d", "h3d", 16, 0, 16, 256, 0, 256, 200, 0, 200 );
  h3d->GetXaxis()->SetTitle( "tile_id" );
  h3d->GetYaxis()->SetTitle( "strip" );
  h3d->GetZaxis()->SetTitle( "adc" );
  tree->Project( h3d->GetName(), var3d, TCut() );

  const TString var2d = ("samples.adc:samples.strip+256*(samples.tile+8*(samples.layer-55))" );
  auto h2d_all = new TH2F( "h2d_all", "h2d_all", 4096, 0, 4096, 200, 0, 200 );
  h2d_all->GetXaxis()->SetTitle( "strip" );
  h2d_all->GetYaxis()->SetTitle( "adc" );
  tree->Project( h2d_all->GetName(), var2d, TCut() );

  {
    const auto cvname = "cv_all";
    auto cv = new TCanvas( cvname, cvname, 900, 900 );
    h2d_all->Draw("colz");
    pdfDocument.Add( cv );
  }

  for( int ilayer = 0; ilayer <2; ++ilayer )
  {
    int layer = 55+ilayer;
    const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
    for( int tile = 0; tile <8; ++tile )
    {
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, tile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      const auto hname = Form( "h_%i_%i", ilayer, tile );
      const auto htitle = Form( "%s - %i,%i", name.c_str(), ilayer, tile );

      int bin = tile+8*ilayer;
      h3d->GetXaxis()->SetRange( bin+1, bin+1 );
      auto h2d = static_cast<TH2*>( h3d->Project3D( "zy" ) );
      h2d->SetName( hname );
      h2d->SetTitle( htitle );

      const auto cvname = Form( "cv_%i_%i", ilayer, tile );
      auto cv = new TCanvas( cvname, cvname, 900, 900 );
      h2d->Draw("colz");
      pdfDocument.Add( cv );
    }
  }

  return pdfFile;

}
