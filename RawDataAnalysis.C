#include <memory>
#include <array>

R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasMapping.h>
#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

#include <boost/container_hash/hash.hpp>

//_____________________________________________________________________________
void RawDataAnalysis( int runNumber = 20446 )
{
  const TString inputFile = Form( "MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  std::cout << "RawDataAnalysis - inputfilename: " << inputfilename << std::endl;

  bool verbosity = false;
  
  // filemanager
  auto tfile = TFile::Open( inputfilename, "READ" );
  if( !tfile )
  { 
    std::cout << "RawDataAnalysis - invalid file: " << inputfilename << std::endl;
    return;
  }
  
  auto tree = static_cast<TTree*>( tfile->Get( "T" ) );
  if( !tree )
  {
    std::cout << "RawDataAnalysis - invalid file: " << inputfilename << std::endl;
    return;
  }
  
  // here you can define all the relevant histograms you need to fill from the TTree content
  
  // create container and link to the TTree content
  // the actual definition of the container can be found in "$OFFLINE_MAIN/include/micromegas/MicromegasRawDataEvaluation.h"
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );
  
  // loop over tree entries
  const int entries = tree->GetEntries();
  std::cout << "RawDataAnalysis - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    tree->GetEntry(i);

    std::cout << "RawDataAnalysis -"
      << " entry: " << i 
      << " samples size: " << container->samples.size()
      << " waveforms size: " << container->waveforms.size()
      << std::endl;

    // loop over all samples for this event samples
    for( const auto& sample:container->samples )
    { 
      if( verbosity )
      {
        std::cout << "RawDataAnalysis -"
          << " sample layer: " << sample.layer
          << " tile: " << sample.tile
          << " strip: " << sample.strip 
          << " adc: " << sample.adc
          << " sample: " << sample.sample
          << std::endl;
      }
      /*
       * here you can put whatever criteria on the samples needed to fill or not
       * the histogram that you have define, e.g. timing information, 
       * number counts above threshold for a given sample, etc.
       */
    }

    // loop over waveforms
    /* 
     * waveforms already contains the information relevant from the sample corresponding to the max adc value
     * in the entire sample range 
     * it is therefore more convenient to fill hit-based histograms (like hit profile)
     */
    for( const auto& waveform:container->waveforms )
    { 
      if( verbosity )
      {
        std::cout << "RawDataAnalysis -"
          << " waveform layer: " << waveform.layer
          << " tile: " << waveform.tile
          << " strip: " << waveform.strip 
          << " adc_max: " << waveform.adc_max
          << " sample_max: " << waveform.sample_max
          << std::endl;
      }
      /*
       * here you can put whatever criteria on the samples needed to fill or not
       * the histogram that you have define, e.g. timing information, 
       * number counts above threshold for a given sample, etc.
       */
    }
  }

}
