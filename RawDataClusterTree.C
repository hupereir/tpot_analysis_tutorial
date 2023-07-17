#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>

#include <memory>
#include <array>

R__LOAD_LIBRARY(libmicromegas.so)

#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasRawDataEvaluation.h>

//____________________________________________________________________________
bool operator < (const MicromegasRawDataEvaluation::Waveform& lhs, const MicromegasRawDataEvaluation::Waveform& rhs )
{ return lhs.strip < rhs.strip; }

//____________________________________________________________________________
class FullEventContainer: public TObject
{
  
  public:
  
  void Reset()
  {
    lvl1_bco = 0;
    lvl1_counter = 0;
    n_waveforms_all = 0;
    n_waveforms_signal = 0;
    n_clusters = 0;
    n_phi_clusters = 0;
    n_z_clusters = 0;
    n_detector_clusters.clear();
  }
  
  /// ll1 bco
  uint64_t lvl1_bco = 0;

  /// ll1 counter
  uint32_t lvl1_counter = 0;
  
  // number of waveforms
  unsigned short n_waveforms_all = 0;
  
  // number of signal waveforms
  unsigned short n_waveforms_signal = 0;
  
  // number of clusters
  unsigned short n_clusters = 0;
  
  // number of clusters per detector
  std::vector<unsigned short> n_detector_clusters;
  
  // number of clusters
  unsigned short n_phi_clusters = 0;
  
  // number of clusters
  unsigned short n_z_clusters = 0;
    
  ClassDef(FullEventContainer,1)
  
};

namespace
{
  
  class cluster_t: public std::vector<MicromegasRawDataEvaluation::Waveform>
  {
    
    public:
    bool is_signal = true;
  };
  
  // sort waveforms per detid
  using cluster_list_t = std::vector<cluster_t>;
  using waveform_set_t = std::set<MicromegasRawDataEvaluation::Waveform>;
  using waveform_map_t = std::map<int,waveform_set_t>;
  
  // number of sigma above pedestal to remove noise
  static constexpr double n_sigma = 5.;
  
  // output tfile
  std::unique_ptr<TFile> tfile_out;
 
  // detector names
  // running container
  FullEventContainer* m_fullEventContainer = nullptr;

  // output tree
  TTree* tree_out = nullptr;
  
  // calibration data
  MicromegasCalibrationData calibration_data;
   
  // create output tree
  void create_tree( const TString& rootfilename )
  {
    // create output TFile
    tfile_out.reset(TFile::Open(rootfilename, "RECREATE"));
    tree_out = new TTree( "T", "T" );
    m_fullEventContainer = new FullEventContainer;
    tree_out->Branch( "Event", &m_fullEventContainer );
  }
  
  // save output tree
  void save_tree()
  {
    tfile_out->cd();
    tree_out->Write();
    tfile_out->Close();
  }
  
  // sample windows to determine signal and background clusters
  using sample_window_t = std::pair<unsigned short, unsigned short>;
  sample_window_t sample_window_signal = {20, 45};
  sample_window_t sample_window_background = {5, 15 };
      
  //_________________________________________________________________
  cluster_list_t get_clusters( const waveform_set_t& waveform_set )
  {
    // rolling cluster
    cluster_list_t clusters;
    cluster_t cluster;
    int previous_strip = -1;
    for( const auto& waveform:waveform_set )
    {
      if( previous_strip >= 0 && waveform.strip > previous_strip+1 )
      { 
        if( !cluster.empty() ) clusters.push_back(cluster);
        cluster.clear();
      }
      
      previous_strip = waveform.strip;
      cluster.push_back( waveform );
    }
    
    // store last cluster
    if( !cluster.empty() ) clusters.push_back( cluster );
    
    return clusters;
  }

  //_________________________________________________________________
  bool is_signal( const MicromegasRawDataEvaluation::Waveform& waveform, const sample_window_t& sample_window )
  {
    return 
      waveform.rms > 0 &&
      waveform.sample_max >= sample_window.first &&
      waveform.sample_max < sample_window.second &&
      waveform.adc_max > waveform.pedestal+n_sigma*waveform.rms;
  }
  
  //_________________________________________________________________
  bool is_signal( const MicromegasRawDataEvaluation::Waveform& waveform )
  { return is_signal( waveform, sample_window_signal ); }

  //_________________________________________________________________
  bool is_background( const MicromegasRawDataEvaluation::Waveform& waveform )
  { return is_signal( waveform, sample_window_background ); }
  
  //_________________________________________________________________
  void process_event( uint64_t lvl1_bco, uint32_t lvl1_counter, const MicromegasRawDataEvaluation::Waveform::List& waveforms )
  {

    std::cout
      << "process_event -"
      << " lvl1_bco: " << lvl1_bco 
      << " lvl1_counter: " << lvl1_counter 
      << " waveforms: " << waveforms.size() 
      << std::endl;

    // reset full event container
    m_fullEventContainer->Reset();
    m_fullEventContainer->lvl1_bco = lvl1_bco;
    m_fullEventContainer->lvl1_counter = lvl1_counter;
    m_fullEventContainer->n_waveforms_all = waveforms.size();
      
    waveform_map_t waveform_map;
    for( const auto& waveform:waveforms )
    {
      if( !waveform.layer ) continue;
      const int detid = waveform.tile + 8*(waveform.layer-55 );
      if( is_signal( waveform ) )
      {
        waveform_map[detid].insert( waveform );
        ++m_fullEventContainer->n_waveforms_signal;
      }
    }
        
    // store cluster multiplicity
    m_fullEventContainer->n_detector_clusters = std::vector<unsigned short>(16, 0);
    m_fullEventContainer->n_clusters = 0;
    for( const auto& [detid, waveform_set]:waveform_map )
    {
      auto clusters = get_clusters( waveform_set );
      for( auto&& cluster:clusters ) cluster.is_signal = true;
      
      // process clusters
      m_fullEventContainer->n_detector_clusters[detid] = clusters.size();
      m_fullEventContainer->n_clusters += clusters.size();
      
      // per view clusters
      if( detid < 8 ) m_fullEventContainer->n_phi_clusters += clusters.size();
      else m_fullEventContainer->n_z_clusters += clusters.size();
    } 
    
    tree_out->Fill();
    
  }
  
}

//_____________________________________________________________________________
TString RawDataClusterTree( int runNumber = 20445 )
{
  const TString inputfilename = Form( "MicromegasRawDataEvaluation-%08i-0000.root", runNumber );
  const TString rootfilename = Form( "RawDataClusterTree-%08i-0000.root", runNumber );
  
  std::cout << "RawDataClusterTree - inputfilename: " << inputfilename << std::endl;
  std::cout << "RawDataClusterTree - rootfilename: " << rootfilename << std::endl;

  // filemanager
  auto tfile = std::unique_ptr<TFile>( TFile::Open( inputfilename ));
  auto tree = static_cast<TTree*>( tfile->Get( "T" ) );
  if( !tree )
  {
    std::cout << "RawDataClusterTree - invalid file" << std::endl;
    return rootfilename;
  }
  
  auto container = new MicromegasRawDataEvaluation::Container;
  tree->SetBranchAddress( "Event", &container );

  // output tree
  create_tree( rootfilename );
  
  // keep track of all waveforms associated to a given bco 
  MicromegasRawDataEvaluation::Waveform::List waveforms;
  
  // previous bco and counter
  uint64_t prev_lvl1_bco = 0;
  uint32_t prev_lvl1_counter = 0;
  
  // loop over tree entries
  const int entries = tree->GetEntries();
  // const int entries = 1000;
  std::cout << "RawDataClusterTree - entries: " << entries << std::endl;
  for( int i = 0; i < entries; ++i )
  {
    // some printout
    if( !(i%100) )
    { std::cout << "RawDataClusterTree - entry: " << i << std::endl; }
    
    tree->GetEntry(i);    

    const bool consistent = (container->lvl1_count_list.size()==2 && container->lvl1_count_list[0]==container->lvl1_count_list[1]);
    if( !consistent ) 
    {
      std::cout 
        << "RawDataClusterTree -"
        << " entry: " << i 
        << " lvl1_count_list[0]: " << container->lvl1_count_list[0]
        << " lvl1_count_list[1]: " << container->lvl1_count_list[1]
        << " inconsistent" << std::endl;
      waveforms.clear();
      continue;
    }
    
    const auto& lvl1_bco = container->lvl1_bco_list[0];
    const auto& lvl1_counter = container->lvl1_count_list[0];
    
    if( lvl1_bco != prev_lvl1_bco )
    {
      if( !waveforms.empty() ) process_event( prev_lvl1_bco, prev_lvl1_counter, waveforms );
      waveforms.clear();
      prev_lvl1_bco = lvl1_bco;
      prev_lvl1_counter = lvl1_counter;
    }
    
    // copy current waveforms into container
    std::copy( container->waveforms.begin(), container->waveforms.end(), std::back_inserter( waveforms ) );
  }

  // process last event 
  if( !waveforms.empty() ) process_event( prev_lvl1_bco, prev_lvl1_counter, waveforms );

  // save output tree
  save_tree();
  std::cout << "done." << std::endl;
  
  return rootfilename;
}
