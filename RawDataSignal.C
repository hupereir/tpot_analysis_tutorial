
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
TString RawDataSignal(int runNumber = 9443)
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

  const TString var( "adc:strip" );
  const TString var3d = Form( "%s:(tile+8*(layer-55))", var.Data() );
  auto h3d = new TH3I( "h3d", "h3d", 16, 0, 16, 256, 0, 256, 1024, 0, 1024 );
  h3d->GetXaxis()->SetTitle( "tile_id" );
  h3d->GetYaxis()->SetTitle( "strip" );
  h3d->GetZaxis()->SetTitle( "adc" );
  tree->Project( h3d->GetName(), var3d, TCut() );  
    
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
