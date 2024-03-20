#include <TH1.h>
#include <TF1.h>

#include <memory>

R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasCombinedDataEvaluation.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

const std::string status = "Internal";

namespace
{
  MicromegasMapping mapping;

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

  TLatex* put_text( Double_t x, Double_t y, TString value, Double_t fontSize )
  {
    auto text = new TLatex;
    text->SetTextSize( fontSize );
    text->DrawLatex( x, y, value );
    return text;
  }

  //__________________________________________________
  TLine* vertical_line( TVirtualPad* pad, Double_t x )
  {
    double yMin = pad->GetUymin();
    double yMax = pad->GetUymax();

    if( pad->GetLogy() )
    {
      yMin = TMath::Power( 10, yMin );
      yMax = TMath::Power( 10, yMax );
    }

    auto line = new TLine( x, yMin, x, yMax );
    line->SetLineStyle( 2 );
    line->SetLineWidth( 1 );
    line->SetLineColor( 1 );
    return line;

  }

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

  //______________________________________________________
  void draw_vertical_lines( TVirtualPad* pad )
  {
    pad->Update();
    for( int i = 1; i < 16; ++i )
    { vertical_line( pad, 256*i )->Draw(); }
  }

  //______________________________________________________
  std::vector<std::string> get_detector_names()
  {
    // get detector names that match tile/layer
    std::vector<std::string> detector_names;
    for( int ilayer = 0; ilayer < 2; ++ilayer )
      for( int itile = 0; itile < 8; ++itile )
    {
      const int layer = 55+ilayer;
      const auto segmentation = (layer==55) ? MicromegasDefs::SegmentationType::SEGMENTATION_PHI : MicromegasDefs::SegmentationType::SEGMENTATION_Z;
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation, itile );
      const auto name = mapping.get_detname_sphenix_from_hitsetkey( hitsetkey );
      detector_names.push_back( std::move( name ) );
    }
    return detector_names;
  }

  //______________________________________________________
  void draw_detector_names( TVirtualPad* pad, std::vector<std::string> detector_names, Double_t fontSize = 0.04 )
  {
    pad->Update();
    const Double_t xMin = pad->GetUxmin();
    const Double_t xMax = pad->GetUxmax();
    const Double_t yMax = pad->GetUymax();
    const auto dx = (xMax-xMin)/16;

    for( int i = 0; i < 16; ++i )
    {
      auto x = xMin + 0.1*dx + i*dx;
      put_text( x, 0.8*yMax, detector_names[i].c_str(), fontSize )->Draw();
    }
  }

  //______________________________________________________
  void draw_information( int runNumber )
  {

    if( false )
    {
      auto text = new TLatex;
      text->SetNDC( true );
      text->SetTextSize( 0.045 );
      text->DrawLatex( 0.77, 0.95, "#it{09/14/2023}" );
    }

    {
      auto text = new TPaveText(0.17,0.71,0.53,0.86, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( "#it{#bf{sPHENIX}}" );
      // text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
      /// text->AddText(Form( "Run %i (Cosmics)", runNumber));
      text->AddText("TPOT, Cosmic data");
      text->Draw();
    }
  }

}

//_____________________________________________________________________________
void NoiseEvaluation(
  int runNumber = 34567
)
{
  gStyle->SetOptStat(0);

  const TString inputFile = Form( "MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  const TString pdfFile = Form( "NoiseEvaluation-%08i-0000.pdf", runNumber );

  std::cout << "NoiseEvaluation - inputFile: " << inputFile << std::endl;
  std::cout << "NoiseEvaluation - pdfFile: " << pdfFile << std::endl;

  const auto detector_names = get_detector_names();

  PdfDocument pdfDocument( pdfFile );

  auto tfile = TFile::Open( inputFile );
  auto tree = static_cast<TTree*>( tfile->Get( "T" ) );

  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // get 2d histogram
  // auto h_adc_channel = new TH2I( "h_adc_channel", "", 4096, 0, 4096, 500, 0, 500 );
  auto h_adc_channel = new TH2I( "h_adc_channel", "", 4096, 0, 4096, 250, 0, 250 );
  h_adc_channel->GetXaxis()->SetTitle( "strip" );
  h_adc_channel->GetYaxis()->SetTitle( "adc" );
  h_adc_channel->GetYaxis()->SetTitleOffset( 1.3 );

  // also create 3D histogram for a per-detector view
  auto h_adc_channel_3d = new TH3F( "h_adc_channel_3d", "h_adc_channel_3d", 16, 0, 16, 256, 0, 256, 250, 0, 250 );
  h_adc_channel_3d->GetXaxis()->SetTitle( "tile_id" );
  h_adc_channel_3d->GetYaxis()->SetTitle( "strip" );
  h_adc_channel_3d->GetZaxis()->SetTitle( "adc" );

  auto p_adc_channel = new TProfile( "p_adc_channel", "", 4096, 0, 4096 );
  p_adc_channel->SetErrorOption("s");

  // loop over tree entries
  const int entries = std::min<int>( 1000, tree->GetEntries() );
  std::cout << "NoiseEvaluation - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    // some printout
    if( !(i%100) )
    { std::cout << "NoiseEvaluation - entry: " << i << std::endl; }

    tree->GetEntry(i);

    // loop over samples
    for( const auto& sample:container->samples )
    {
      if( !( sample.sample>=5 && sample.sample<=15 ) ) continue;
      const auto detid = sample.tile + 8*(sample.layer-55);
      const auto strip = sample.strip;
      const auto channel = sample.strip + 256*(sample.tile + 8*(sample.layer-55));
      const auto adc = sample.adc;

      h_adc_channel->Fill( channel, adc );
      p_adc_channel->Fill( channel, adc );
      h_adc_channel_3d->Fill( detid, strip, adc );
    }

  }

  std::vector<double> pedestal_array;
  std::vector<double> rms_array;

  // get mean and rms histogram
  auto h_pedestal = new TH1F( "h_pedestal", "", 4096, 0, 4096 );
  h_pedestal->GetXaxis()->SetTitle( "strip" );
  h_pedestal->GetYaxis()->SetTitle( "adc" );
  h_pedestal->GetYaxis()->SetTitleOffset( 1.3 );

  auto h_rms = new TH1F( "h_rms", "", 4096, 0, 4096 );
  h_rms->GetXaxis()->SetTitle( "strip" );
  h_rms->GetYaxis()->SetTitle( "adc" );
  h_rms->GetYaxis()->SetTitleOffset( 1.3 );

  for( int i = 0; i < h_adc_channel->GetNbinsX(); ++i )
  {
    const auto pedestal = p_adc_channel->GetBinContent(i+1);
    const auto rms = p_adc_channel->GetBinError(i+1);
    h_pedestal->SetBinContent( i+1,pedestal );
    h_rms->SetBinContent(i+1, rms );

    pedestal_array.push_back(pedestal);
    rms_array.push_back(rms);

  }

  // second loop to remove pedestal
  auto h_adc_channel_corrected = new TH2I( "h_adc_channel_corrected", "", 4096, 0, 4096, 400, -200, 200 );
  h_adc_channel_corrected->GetXaxis()->SetTitle( "strip" );
  h_adc_channel_corrected->GetYaxis()->SetTitle( "adc-pedestal" );
  h_adc_channel_corrected->GetYaxis()->SetTitleOffset( 1.3 );

  auto h_adc_corrected = new TH1F( "h_adc_corrected", "", 400, -200, 200 );
  h_adc_corrected->SetFillStyle(1001);
  h_adc_corrected->SetFillColor(kYellow );
  h_adc_corrected->GetXaxis()->SetTitle( "ADC" );

  for( int i = 0; i < entries; ++i )
  {
    // some printout
    if( !(i%100) )
    { std::cout << "NoiseEvaluation - entry: " << i << std::endl; }

    tree->GetEntry(i);

    // loop over samples
    for( const auto& sample:container->samples )
    {
      if( !( sample.sample>=5 && sample.sample<=15 ) ) continue;
      const auto detid = sample.tile + 8*(sample.layer-55);
      const auto channel = sample.strip + 256*(sample.tile + 8*(sample.layer-55));
      const double adc = sample.adc;
      const double adc_corrected = adc - pedestal_array[channel];
      h_adc_channel_corrected->Fill( channel, adc_corrected );

      if( detid != 8 ) h_adc_corrected->Fill( adc_corrected );

    }
  }

  // plots
  if( true )
  {
    auto cv = new TCanvas( "cv0", "cv0", 980, 900 );
    //  gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);
    h_adc_channel->Draw("col");

    draw_information( runNumber );

    // draw_vertical_lines( gPad );
    // draw_detector_names( gPad, detector_names, 0.02 );
    pdfDocument.Add(cv);
  }

  if( true )
  {
    auto cv = new TCanvas( "cv0", "cv0", 980, 900 );
    cv->Divide( 4, 4 );

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
        h_adc_channel_3d->GetXaxis()->SetRange( bin+1, bin+1 );
        auto h2d = static_cast<TH2*>( h_adc_channel_3d->Project3D( "zy" ) );
        h2d->SetName( hname );
        h2d->SetTitle( htitle );

        cv->cd( bin+1 );
        h2d->Draw("colz");
      }
    }

    pdfDocument.Add( cv );
  }

  if( true )
  {
    auto cv = new TCanvas( "cv0", "cv0", 980, 900 );
   //  gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);

    h_adc_channel_corrected->Draw("col");

    draw_information( runNumber );

    draw_vertical_lines( gPad );
    // draw_detector_names( gPad, detector_names, 0.02 );
    pdfDocument.Add(cv);
  }


  if( true )
  {
    auto cv = new TCanvas( "cv2", "cv2", 980, 900 );
   //  gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);

    h_adc_corrected->GetXaxis()->SetRangeUser(-100,100);
    // h_adc_corrected->Fit("gaus");
    h_adc_corrected->Draw();
    gPad->SetLogy(true);

    {
      auto text = new TPaveText(0.15,0.75,0.61,0.90, "NDC" );
      text->SetFillColor(0);
      text->SetFillStyle(0);
      text->SetBorderSize(0);
      text->SetTextAlign(11);
      text->AddText( "#it{#bf{sPHENIX}}" );
      // text->AddText( Form( "#it{#bf{sPHENIX}} %s", status.c_str() ));
      // text->AddText(Form( "Run %i (Cosmics)", runNumber));
      text->AddText("TPOT, Cosmic data");
      text->AddText(Form("#LTRMS#GT=%.1f ADC",  h_adc_corrected->GetRMS() ) );
      text->Draw();
    }

    pdfDocument.Add(cv);

    std::cout << "rms: " << h_adc_corrected->GetRMS() << std::endl;

  }

  if( true )
  {
    auto cv = new TCanvas( "cv1", "cv1", 980, 900 );
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);

    h_pedestal->SetLineColor(1);
    h_pedestal->Draw();


    draw_information( runNumber );

    draw_vertical_lines( gPad );
    draw_detector_names( gPad, detector_names, 0.02 );
    pdfDocument.Add(cv);

  }

  if( true )
  {
    auto cv = new TCanvas( "cv1", "cv1", 980, 900 );
    gPad->SetTopMargin( 0.07 );
    gPad->SetLeftMargin( 0.14);

    h_rms->SetLineColor(1);
    h_rms->SetMinimum(0);
    h_rms->SetMaximum(40);
    h_rms->Draw();

    draw_information( runNumber );

    draw_vertical_lines( gPad );
    draw_detector_names( gPad, detector_names, 0.02 );
    auto line = horizontal_line( gPad, 10 );
    line->SetLineColor(2);
    line->Draw();
    pdfDocument.Add(cv);
  }

}
