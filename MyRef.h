#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "TError.h"
#include "TRandom.h"
#include <vector>
#include <map>
#include "TString.h"

using namespace std ;

  namespace {
    typedef pair<Double_t, Int_t> keys;
  }

//______________________________________________________________________________
struct StRefMultCorr {
    // Specify the type of multiplicity (default is refmult)
    // "refmult"   - reference multiplicity defined in |eta|<0.5
    // "refmult2"  - reference multiplicity defined in 0.5<|eta|<1.0
    // "refmult3"  - reference multiplicity defined in |eta|<0.5 without protons
    // "toftray"   - TOF tray multiplicity
    StRefMultCorr(const TString name="refmult");
    virtual ~StRefMultCorr(); /// Default destructor

    // Bad run rejection
    Bool_t isBadRun(const Int_t RunId) ;

    // Event-by-event initialization. Call this function event-by-event
    //   * Default ZDC coincidence rate = 0 to make the function backward compatible 
    //   --> i.e. no correction will be applied unless users set the values for 200 GeV
    void initEvent(const UShort_t RefMult, const Double_t z,
        const Double_t zdcCoincidenceRate=0.0) ; // Set multiplicity, vz and zdc coincidence rate

    /// Get corrected multiplicity, correction as a function of primary z-vertex
    Double_t getRefMultCorr() const;

    // Corrected multiplity
    // flag=0:  Luminosity only
    // flag=1:  z-vertex only
    // flag=2:  full correction (default)
    Double_t getRefMultCorr(const UShort_t RefMult, const Double_t z,
       	const Double_t zdcCoincidenceRate, const UInt_t flag=2) const ;

    /// Get 16 centrality bins (5% increment, 0-5, 5-10, ..., 75-80)
    Int_t getCentralityBin16() const;

    /// Get 9 centrality bins (10% increment except for 0-5 and 5-10)
    Int_t getCentralityBin9() const;

    /// Re-weighting correction, correction is only applied up to mNormalize_step (energy dependent)
    Double_t getWeight() const;

    // Initialization of centrality bins etc
    void init(const Int_t RunId);

    // Return begin/end run from energy and year
    Int_t getBeginRun(const Double_t energy, const Int_t year) ;
    Int_t getEndRun(const Double_t energy, const Int_t year) ;

    // Print all parameters
    void print(const Option_t* option="") const ;

    const TString mName ; // refmult, refmult2, refmult3 or toftray (case insensitive)

    // Functions
    void read() ; /// Read input parameters from text file StRoot/StRefMultCorr/Centrality_def.txt
    void readBadRuns() ; /// Read bad run numbers
    void clear() ; /// Clear all arrays
    Bool_t isIndexOk() const ; /// 0 <= mParameterIndex < maxArraySize
    Bool_t isZvertexOk() const ; /// mStart_zvertex < z < mStop_zvertex
    Bool_t isRefMultOk() const ; /// 0-80%, (corrected multiplicity) > mCentrality_bins[0]
    Bool_t isCentralityOk(const Int_t icent) const ; /// centrality bin check
    Int_t setParameterIndex(const Int_t RunId) ; /// Parameter index from run id (return mParameterIndex)

    // Get table name based on the input multiplicity definition
    const Char_t* getTable() const ;

    // Data members
    enum {
        mNCentrality   = 16, /// 16 centrality bins, starting from 75-80% with 5% bin width
        mNPar_z_vertex = 8,
        mNPar_weight   = 6,
        mNPar_luminosity = 2
    };

    // Use these variables to avoid varying the corrected multiplicity
    // in the same event by random numbers
    UShort_t mRefMult ;     /// Current multiplicity
    Double_t mVz ;          /// Current primary z-vertex
    Double_t mZdcCoincidenceRate ; /// Current ZDC coincidence rate
    Double_t mRefMult_corr; /// Corrected refmult

    std::vector<Int_t> mYear              ; /// Year
    std::vector<Int_t> mStart_runId       ; /// Start run id
    std::vector<Int_t> mStop_runId        ; /// Stop run id
    std::vector<Double_t> mStart_zvertex  ; /// Start z-vertex (cm)
    std::vector<Double_t> mStop_zvertex   ; /// Stop z-vertex (cm)
    std::vector<Double_t> mNormalize_stop ; /// Normalization between MC and data (normalized in refmult>mNormalize_stop)
    std::vector<Int_t> mCentrality_bins[mNCentrality+1] ; /// Centrality bins (last value is set to 5000)
    std::vector<Double_t> mPar_z_vertex[mNPar_z_vertex] ; /// parameters for z-vertex correction
    std::vector<Double_t> mPar_weight[mNPar_weight] ; /// parameters for weight correction
    std::vector<Double_t> mPar_luminosity[mNPar_luminosity] ; /// parameters for luminosity correction (valid only for 200 GeV)
    Int_t mParameterIndex; /// Index of correction parameters

    std::multimap<std::pair<Double_t, Int_t>, Int_t> mBeginRun ; /// Begin run number for a given (energy, year)
    std::multimap<std::pair<Double_t, Int_t>, Int_t> mEndRun   ; /// End run number for a given (energy, year)
    std::vector<Int_t> mBadRun ; /// Bad run number list

};

//______________________________________________________________________________
// Default constructor
StRefMultCorr::StRefMultCorr(const TString name)
 : mName(name)
{
  mRefMult = 0 ;
  mVz = -9999. ;
  mRefMult_corr = -1.0 ;

  // Clear all data members
  clear() ;

  // Read parameters
  read() ;
  readBadRuns() ;
}

//______________________________________________________________________________
// Default destructor
StRefMultCorr::~StRefMultCorr()
{
}

//______________________________________________________________________________
Int_t StRefMultCorr::getBeginRun(const Double_t energy, const Int_t year)
{
  keys key(std::make_pair(energy, year));

  // Make sure key exists
  multimap<keys, Int_t>::iterator iterCheck = mBeginRun.find(key);
  if ( iterCheck == mBeginRun.end() ) {
    Error("StRefMultCorr::getBeginRun", "can't find energy = %1.1f, year = %d", energy, year);
    return -1;
  }

  pair<multimap<keys, Int_t>::iterator, multimap<keys, Int_t>::iterator> iterRange = mBeginRun.equal_range(key);

  return (*(iterRange.first)).second ;
}

//______________________________________________________________________________
Int_t StRefMultCorr::getEndRun(const Double_t energy, const Int_t year)
{
  keys key(std::make_pair(energy, year));

  // Make sure key exists
  multimap<keys, Int_t>::iterator iterCheck = mEndRun.find(key);
  if ( iterCheck == mEndRun.end() ) {
    Error("StRefMultCorr::getEndRun", "can't find energy = %1.1f, year = %d", energy, year);
    return -1;
  }

  pair<multimap<keys, Int_t>::iterator, multimap<keys, Int_t>::iterator> iterRange = mEndRun.equal_range(key);
  multimap<keys, Int_t>::iterator iter = iterRange.second ;
  iter--;

  return (*iter).second ;
}

//______________________________________________________________________________
void StRefMultCorr::clear()
{
  // Clear all arrays, and set parameter index = -1

  mYear.clear() ;
  mStart_runId.clear() ;
  mStop_runId.clear() ;
  mStart_zvertex.clear() ;
  mStop_zvertex.clear() ;
  mNormalize_stop.clear() ;

  for(Int_t i=0;i<mNCentrality;i++) {
    mCentrality_bins[i].clear() ;
  }
  mParameterIndex = -1 ;

  for(Int_t i=0;i<mNPar_z_vertex;i++) {
      mPar_z_vertex[i].clear() ;
  }

  for(Int_t i=0;i<mNPar_weight;i++) {
      mPar_weight[i].clear();
  }

  for(Int_t i=0;i<mNPar_luminosity;i++) {
      mPar_luminosity[i].clear();
  }

  mBeginRun.clear() ;
  mEndRun.clear() ;
  mBadRun.clear() ;
}

//______________________________________________________________________________
Bool_t StRefMultCorr::isBadRun(const Int_t RunId)
{
  // Return true if a given run id is bad run
  vector<Int_t>::iterator iter = std::find(mBadRun.begin(), mBadRun.end(), RunId);
#if 0
  if ( iter != mBadRun.end() ) {
    // QA
    cout << "StRefMultCorr::isBadRun  Find bad run = " << (*iter) << endl;
  }
#endif

  return ( iter != mBadRun.end() ) ;
}

//______________________________________________________________________________
void StRefMultCorr::initEvent(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate)
{
  // Set refmult, vz and corrected refmult if current (refmult,vz) are different from inputs
  // User must call this function event-by-event before 
  // calling any other public functions
  if ( mRefMult != RefMult || mVz != z || mZdcCoincidenceRate != zdcCoincidenceRate ) {
    mRefMult            = RefMult ;
    mVz                 = z ;
    mZdcCoincidenceRate = zdcCoincidenceRate ;
    mRefMult_corr       = getRefMultCorr(mRefMult, mVz, mZdcCoincidenceRate) ;
  }
}

//______________________________________________________________________________
Bool_t StRefMultCorr::isIndexOk() const
{
  // mParameterIndex not initialized (-1)
  if ( mParameterIndex == -1 ) {
    Error("StRefMultCorr::isIndexOk", "mParameterIndex = -1. Call init(const Int_t RunId) function to initialize centrality bins, corrections");
    Error("StRefMultCorr::isIndexOk", "mParameterIndex = -1. or use valid run numbers defined in Centrality_def_%s.txt", mName.Data());
    Error("StRefMultCorr::isIndexOk", "mParameterIndex = -1. exit");
    cout << endl;
    // Stop the process if invalid run number found
    exit(0);
  }

  // Out of bounds
  if ( mParameterIndex >= (Int_t)mStart_runId.size() ) {
    Error("StRefMultCorr::isIndexOk",
        Form("mParameterIndex = %d > max number of parameter set = %d. Make sure you put correct index for this energy",
          mParameterIndex, mStart_runId.size()));
    return kFALSE ;
  }

  return kTRUE ;
}

//______________________________________________________________________________
Bool_t StRefMultCorr::isZvertexOk() const
{
  // Primary z-vertex check
  return ( mVz > mStart_zvertex[mParameterIndex] && mVz < mStop_zvertex[mParameterIndex] ) ;
}

//______________________________________________________________________________
Bool_t StRefMultCorr::isRefMultOk() const
{
  // Invalid index
  if ( !isIndexOk() ) return kFALSE ;

  // select 0-80%
  return (mRefMult_corr > mCentrality_bins[0][mParameterIndex] && mRefMult_corr < mCentrality_bins[mNCentrality][mParameterIndex]);
}

//______________________________________________________________________________
Bool_t StRefMultCorr::isCentralityOk(const Int_t icent) const
{
  // Invalid centrality id
  if ( icent < -1 || icent >= mNCentrality+1 ) return kFALSE ;

  // Invalid index
  if ( !isIndexOk() ) return kFALSE ;

  // Special case
  // 1. 80-100% for icent=-1
  if ( icent == -1 ) return (mRefMult_corr <= mCentrality_bins[0][mParameterIndex]);

  // 2. icent = mNCentrality
  if ( icent == mNCentrality ) return (mRefMult_corr <= mCentrality_bins[mNCentrality][mParameterIndex]);

  const Bool_t ok = (mRefMult_corr > mCentrality_bins[icent][mParameterIndex] && mRefMult_corr <= mCentrality_bins[icent+1][mParameterIndex]);
//  if(ok){
//    cout << "StRefMultCorr::isCentralityOk  refmultcorr = " << mRefMult_corr
//      << "  min. bin = " << mCentrality_bins[icent][mParameterIndex]
//      << "  max. bin = " << mCentrality_bins[icent+1][mParameterIndex]
//      << endl;
//  }
  return ok ;
}

//______________________________________________________________________________
void StRefMultCorr::init(const Int_t RunId)
{
  // Reset mParameterIndex
  mParameterIndex = -1 ;

  // call setParameterIndex
  setParameterIndex(RunId) ;
}

//______________________________________________________________________________
Int_t StRefMultCorr::setParameterIndex(const Int_t RunId)
{
  // Determine the corresponding parameter set for the input RunId
  for(UInt_t npar = 0; npar < mStart_runId.size(); npar++)
  {
    if(RunId >= mStart_runId[npar] && RunId <= mStop_runId[npar])
    {
      mParameterIndex = npar ;
//      cout << "StRefMultCorr::setParameterIndex  Parameter set = " << mParameterIndex << " for RUN " << RunId << endl;
      break ;
    }
  }

  if(mParameterIndex == -1){
    Error("StRefMultCorr::setParameterIndex", "Parameter set does not exist for RUN %d", RunId);
  }
  //else cout << "Parameter set = " << npar_set << endl;

  return mParameterIndex ;
}

//______________________________________________________________________________
Double_t StRefMultCorr::getRefMultCorr() const
{
  // Call initEvent() first
  return mRefMult_corr ;
}

//______________________________________________________________________________
Double_t StRefMultCorr::getRefMultCorr(const UShort_t RefMult, const Double_t z,
    const Double_t zdcCoincidenceRate, const UInt_t flag) const
{
  // Apply correction if parameter index & z-vertex are ok
  if (!isIndexOk() || !isZvertexOk()) return RefMult ;

  // Correction function for RefMult, takes into account z_vertex dependence

  // Luminosity corrections
  // 200 GeV only. correction = 1 for all the other energies
  const Double_t par0l = mPar_luminosity[0][mParameterIndex] ;
  const Double_t par1l = mPar_luminosity[1][mParameterIndex] ;
  const Double_t correction_luminosity = (par0l==0.0) ? 1.0 : 1.0/(1.0 + par1l/par0l*zdcCoincidenceRate/1000.);

  // par0 to par5 define the parameters of a polynomial to parametrize z_vertex dependence of RefMult
  const Double_t par0 = mPar_z_vertex[0][mParameterIndex];
  const Double_t par1 = mPar_z_vertex[1][mParameterIndex];
  const Double_t par2 = mPar_z_vertex[2][mParameterIndex];
  const Double_t par3 = mPar_z_vertex[3][mParameterIndex];
  const Double_t par4 = mPar_z_vertex[4][mParameterIndex];
  const Double_t par5 = mPar_z_vertex[5][mParameterIndex];
  const Double_t par6 = mPar_z_vertex[6][mParameterIndex];
  const Double_t par7 = mPar_z_vertex[7][mParameterIndex]; // this parameter is usually 0, it takes care for an additional efficiency, usually difference between phase A and phase B parameter 0

  const Double_t  RefMult_ref = par0; // Reference mean RefMult at z=0
  const Double_t  RefMult_z = par0 + par1*z + par2*z*z + par3*z*z*z + par4*z*z*z*z + par5*z*z*z*z*z + par6*z*z*z*z*z*z; // Parametrization of mean RefMult vs. z_vertex position
  Double_t  Hovno = 1.0; // Correction factor for RefMult, takes into account z_vertex dependence

  if(RefMult_z > 0.0)
  {
    Hovno = (RefMult_ref + par7)/RefMult_z;
  }

  Double_t RefMult_d = (Double_t)(RefMult)+gRandom->Rndm(); // random sampling over bin width -> avoid peak structures in corrected distribution
  Double_t RefMult_corr  = -9999. ;
  switch ( flag ) {
    case 0: return RefMult_d*correction_luminosity;
    case 1: return RefMult_d*Hovno;
    case 2: return RefMult_d*Hovno*correction_luminosity;
    default:
      {
        Error("StRefMultCorr::getRefMultCorr", "invalid flag, flag=%d, should be 0,1 or 2", flag);
	return -9999.;
      }
  }
//  cout << "Input RefMult = " << RefMult << ", input z = " << z << ", RefMult_corr = " << RefMult_corr << endl;
  return RefMult_corr ;
}

//______________________________________________________________________________
Double_t StRefMultCorr::getWeight() const
{
  Double_t Weight = 1.0;

  // Invalid index
  if( !isIndexOk() ) return Weight ;

  // Invalid z-vertex
  if( !isZvertexOk() ) return Weight ;

  const Double_t par0 =   mPar_weight[0][mParameterIndex];
  const Double_t par1 =   mPar_weight[1][mParameterIndex];
  const Double_t par2 =   mPar_weight[2][mParameterIndex];
  const Double_t par3 =   mPar_weight[3][mParameterIndex];
  const Double_t par4 =   mPar_weight[4][mParameterIndex];
  const Double_t A    =   mPar_weight[5][mParameterIndex];

  // Additional z-vetex dependent correction
  //const Double_t A = ((1.27/1.21))/(30.0*30.0); // Don't ask...
  //const Double_t A = (0.05/0.21)/(30.0*30.0); // Don't ask...

  if(isRefMultOk() // 0-80%
      && mRefMult_corr < mNormalize_stop[mParameterIndex] // re-weighting only apply up to normalization point
      && mRefMult_corr != -(par3/par2) // avoid denominator = 0
    )
  {
    Weight = par0 + par1/(par2*mRefMult_corr + par3) + par4*(par2*mRefMult_corr + par3); // Parametrization of MC/data RefMult ratio
    Weight = Weight + (Weight-1.0)*(A*mVz*mVz); // z-dependent weight correction
  }

  return Weight;
}

//______________________________________________________________________________
Int_t StRefMultCorr::getCentralityBin16() const
{
  Int_t CentBin16 = -1;

  // Invalid index
  if( !isIndexOk() ) return CentBin16 ;

  while(CentBin16 < mNCentrality && !isCentralityOk(CentBin16) )
  {
    CentBin16++;
  }

  // return -1 if CentBin16 = 16 (very large refmult, refmult>5000)
  return (CentBin16==16) ? -1 : CentBin16;
}

//______________________________________________________________________________
Int_t StRefMultCorr::getCentralityBin9() const
{
  Int_t CentBin9 = -1;

  // Invalid index
  if ( !isIndexOk() ) return CentBin9 ;

  const Int_t CentBin16 = getCentralityBin16(); // Centrality bin 16
  const Bool_t isCentralityOk = CentBin16 >= 0 && CentBin16 < mNCentrality ;

  // No centrality is defined
  if (!isCentralityOk) return CentBin9 ;

  // First handle the exceptions
  if(mRefMult_corr > mCentrality_bins[15][mParameterIndex] && mRefMult_corr <= mCentrality_bins[16][mParameterIndex])
  {
    CentBin9 = 8; // most central 5%
  }
  else if(mRefMult_corr > mCentrality_bins[14][mParameterIndex] && mRefMult_corr <= mCentrality_bins[15][mParameterIndex])
  {
    CentBin9 = 7; // most central 5-10%
  }
  else
  {
    CentBin9 = (Int_t)(0.5*CentBin16);
  }

  return CentBin9;
}

//______________________________________________________________________________
const Char_t* StRefMultCorr::getTable() const
{
  if ( mName.CompareTo("refmult", TString::kIgnoreCase) == 0 ) {
    return "/home/gang/Analysis/Parity/200GeV_run11/StRefMultCorr/Centrality_def_refmult.txt";
  }
  else if ( mName.CompareTo("refmult2", TString::kIgnoreCase) == 0 ) {
    return "/home/gang/Analysis/Parity/200GeV_run11/StRefMultCorr/Centrality_def_refmult2.txt";
  }
  else if ( mName.CompareTo("refmult3", TString::kIgnoreCase) == 0 ) {
    return "/home/gang/Analysis/Parity/200GeV_run11/StRefMultCorr/Centrality_def_refmult3.txt";
  }
  else if ( mName.CompareTo("toftray", TString::kIgnoreCase) == 0 ) {
    return "/home/gang/Analysis/Parity/200GeV_run11/StRefMultCorr/Centrality_def_toftray.txt";
  }
  else{
    Error("StRefMultCorr::getTable", "No implementation for %s", mName.Data());
    cout << "Current available option is refmult or refmult2 or refmult3 or toftray" << endl;
    return "";
  }
}
//______________________________________________________________________________
void StRefMultCorr::read()
{
  // Open the parameter file and read the data
  const Char_t* inputFileName(getTable());
  ifstream ParamFile(inputFileName);
  if(!ParamFile){
    Error("StRefMultCorr::read", "cannot open %s", inputFileName);
    return;
  }
  cout << "StRefMultCorr::read  Open " << inputFileName << flush ;

  string line ;
  getline(ParamFile,line);

  if(line.find("Start_runId")!=string::npos)
  {
    while(ParamFile.good())
    {
      Int_t year;
      Double_t energy;
      ParamFile >> year >> energy ;

      Int_t startRunId=0, stopRunId=0 ;
      Double_t startZvertex=-9999., stopZvertex=-9999. ;
      ParamFile >> startRunId >> stopRunId >> startZvertex >> stopZvertex ;

      // Error check
      if(ParamFile.eof()) break;

      mYear.push_back(year) ;
      mBeginRun.insert(std::make_pair(std::make_pair(energy, year), startRunId));
      mEndRun.insert(std::make_pair(std::make_pair(energy, year), stopRunId));

      mStart_runId.push_back( startRunId ) ;
      mStop_runId.push_back( stopRunId ) ;
      mStart_zvertex.push_back( startZvertex ) ;
      mStop_zvertex.push_back( stopZvertex ) ;
      for(Int_t i=0;i<mNCentrality;i++) {
        Int_t centralitybins=-1;
        ParamFile >> centralitybins;
        mCentrality_bins[i].push_back( centralitybins );
      }
      Double_t normalize_stop=-1.0 ;
      ParamFile >> normalize_stop ;
      mNormalize_stop.push_back( normalize_stop );

      for(Int_t i=0;i<mNPar_z_vertex;i++) {
          Double_t param=-9999.;
          ParamFile >> param;
          mPar_z_vertex[i].push_back( param );
      }

      for(Int_t i=0;i<mNPar_weight;i++) {
          Double_t param=-9999.;
          ParamFile >> param;
          mPar_weight[i].push_back( param );
      }

      for(Int_t i=0;i<mNPar_luminosity;i++) {
          Double_t param=-9999.;
          ParamFile >> param;
          mPar_luminosity[i].push_back( param );
      }
      mCentrality_bins[mNCentrality].push_back( 5000 );
    }
  }
  else
  {
    cout << endl;
    Error("StRefMultCorr::read", "Input file is not correct! Wrong structure.");
    return;
  }
  ParamFile.close();

  cout << " [OK]" << endl;
}

//______________________________________________________________________________
void StRefMultCorr::readBadRuns()
{
  // Read bad run numbers
  //   - From year 2010 and 2011
  for(Int_t i=0; i<2; i++) {
    cout << "StRefMultCorr::readBadRuns  For " << mName << ": open " << flush ;
    const Int_t year = 2010 + i ;
    const Char_t* inputFileName(Form("/home/gang/Analysis/Parity/200GeV_run11/StRefMultCorr/bad_runs_refmult_year%d.txt", year));
    ifstream fin(inputFileName);
    if(!fin){
      Error("StRefMultCorr::readBadRuns", "can't open %s", inputFileName);
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t runId = 0 ;
    while( fin >> runId ) {
      mBadRun.push_back(runId);
    }
    cout << " [OK]" << endl;
  }
}

//______________________________________________________________________________
void StRefMultCorr::print(const Option_t* option) const
{
  cout << "StRefMultCorr::print  Print input parameters for " << mName << " ========================================" << endl << endl;
  // Option switched off, can be used to specify parameters
//  const TString opt(option);

//  Int_t input_counter = 0;
  for(UInt_t id=0; id<mStart_runId.size(); id++) {
    //cout << "Data line = " << input_counter << ", Start_runId = " << Start_runId[input_counter] << ", Stop_runId = " << Stop_runId[input_counter] << endl;
//    const UInt_t id = mStart_runId.size()-1;

    // Removed line break
    cout << "  Index=" << id;
    cout << Form(" Run=[%8d, %8d]", mStart_runId[id], mStop_runId[id]);
    cout << Form(" z-vertex=[%1.1f, %1.1f]", mStart_zvertex[id], mStop_zvertex[id]);
    cout << ", Normalize_stop=" << mNormalize_stop[id];
    cout << endl;

//    if(opt.IsWhitespace()){
//      continue ;
//    }

    cout << "Centrality:  ";
    for(Int_t i=0;i<mNCentrality;i++){
      cout << Form("  >%2d%%", 80-5*i);
    }
    cout << endl;
    cout << "RefMult:     ";
    for(Int_t i=0;i<mNCentrality;i++){
//      cout << Form("StRefMultCorr::read  Centrality %3d-%3d %%, refmult > %4d", 75-5*i, 80-5*i, mCentrality_bins[i][id]) << endl;
      const TString tmp(">");
      const TString centrality = tmp + Form("%d", mCentrality_bins[i][id]);
      cout << Form("%6s", centrality.Data());
    }
    cout << endl;

    for(Int_t i=0;i<mNPar_z_vertex;i++) {
      cout << "  mPar_z_vertex[" << i << "] = " << mPar_z_vertex[i][id];
    }
    cout << endl;
    for(Int_t i=0;i<mNPar_weight;i++) {
      cout << "  mPar_weight[" << i << "] = " << mPar_weight[i][id];
    }
    cout << endl;
    for(Int_t i=0;i<mNPar_luminosity;i++) {
      cout << "  mPar_luminosity[" << i << "] = " << mPar_luminosity[i][id];
    }
    cout << endl << endl;
  }
  cout << "=====================================================================================" << endl;
}

