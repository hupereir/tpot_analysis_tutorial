#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasHotChannelMapData.h>

namespace
{
  MicromegasHotChannelMapData hot_channels;

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

  //__________________________________________________
  TLine* horizontal_line( TVirtualPad* pad, Double_t y )
  {
    Double_t xMin = pad->GetUxmin();
    Double_t xMax = pad->GetUxmax();

    if( pad->GetLogx() )
    {
      xMin = TMath::Power( 10, xMin );
      xMax = TMath::Power( 10, xMax );
    }

    TLine *line = new TLine( xMin, y, xMax, y );
    line->SetLineStyle( 2 );
    line->SetLineWidth( 1 );
    line->SetLineColor( 1 );
    return line;

  }



}

//_____________________________________________________________________________
void CreateHotChannelMap(int runNumber = 34567)
{
  MicromegasMapping mapping;

  hot_channels.clear();

  const TString inputFile = Form( "MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  const TString pdfFile = Form( "RawDataHitProfile-%08i-0000.pdf", runNumber );
  const TString calibrationFile = Form( "TPOT_HotChannels-%08i-0000.root", runNumber );

  std::cout << "CreateHotChannelMap - inputFile: " << inputFile << std::endl;
  std::cout << "CreateHotChannelMap - pdfFile: " << pdfFile << std::endl;
  std::cout << "CreateHotChannelMap - calibrationFile: " << calibrationFile << std::endl;

  PdfDocument pdfDocument( pdfFile );

  auto tfile = TFile::Open( inputFile );
  auto tree = static_cast<TTree*>(tfile->Get("T"));
  const auto entries = tree->GetEntries();

  const TString hname( "h2d" );
  const TString var2d = "waveforms[].strip:(waveforms[].tile+8*(waveforms[].layer-55))";
  const auto h2d = new TH2I( hname, "", 16, 0, 16, 256, 0, 256 );
  h2d->StatOverflows(true);
  h2d->GetYaxis()->SetTitle( "strip number" );
  h2d->GetZaxis()->SetTitle( "A.U." );

  const TCut signal_cut( "waveforms[].is_signal" );

  tree->Project(h2d->GetName(), var2d, signal_cut);

  // create canvas and divide
  const auto cvname = "cv";
  auto cv = new TCanvas( cvname, cvname, 900, 900 );
  cv->Divide(4,4);

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

      const auto hname = Form( "h_%i_%i", ilayer, itile );
      const auto htitle = Form( "%s - %i,%i", name.c_str(), layer, itile );

      // get the bin matching layer, tile and region in the 3D histgoram
      const int bin = itile + 8*ilayer;
      h2d->GetXaxis()->SetRange( bin+1, bin+1 );

      auto h = static_cast<TH1*>( h2d->ProjectionY() );
      h->Scale( 1./entries );
      h->SetName( hname );
      h->SetTitle( htitle );

      // draw
      cv->cd(bin+1);
      h->SetFillStyle(1001);
      h->SetFillColor(kYellow );
      h->Draw( "hist" );

      // get the mean (Hz)
      const double threshold = 0.2;

      gPad->Update();
      horizontal_line( gPad, threshold)->Draw();

      // loop over bins, store hot channels
      for( int i = 0; i < h->GetNbinsX(); ++i )
      {
        if( h->GetBinContent( i+1 ) > threshold )
        { hot_channels.add_hot_channel( layer, itile, i ); }
      }
    }
  }

  pdfDocument.Add( cv );

  std::cout << hot_channels << "; " << std::endl;
  hot_channels.write( calibrationFile.Data() );

  return pdfFile;

}
