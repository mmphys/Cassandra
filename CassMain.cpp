# include "Cassandra.h"
# include "CassMain.h"

using namespace std;
using namespace Cassandra;

Cassandra::DataManager::DataManager(ModelParams & mp, int iNumReplicas)
: m_Params{ mp }
{
  m_DataSets[BCDMS] = new Cassandra::BCDMSSet(mp, iNumReplicas);
  m_DataSets[JLAB] = new Cassandra::JLABSet(mp, iNumReplicas);
  m_DataSets[SLAC] = new Cassandra::SLACSet(mp, iNumReplicas);
  m_DataSets[BoNuS] = new Cassandra::BoNuSSet(mp, iNumReplicas);
}

Cassandra::DataManager::~DataManager()
{
  for( int i = 0 ; i < Unknown ; i++ )
    delete m_DataSets[i];
}

size_t Cassandra::DataManager::size(void) const noexcept
{
  size_t sz = 0;
  for( int i = 0 ; i < Unknown ; i++ )
    sz += m_DataSets[i]->size();
  return sz;
}

bool Cassandra::DataManager::InitPrediction( const DataSet &/*ds*/, const char * /*pszOutFilePrefix*/ )
{
  return true;
}

// Make a prediction

bool Cassandra::DataManager::MakeOnePrediction(DataSet & ds,
                                               const char * pPDFSet, const char * pMassScheme,
                                               const char * pszOutFilePrefix,
                                               int iSeq, double W2Min, double Q2Min, bool bStrictCut,
                                               size_t SampleSize, bool bShowCut,
                                               double /*xBinSize*/, double /*Q2BinSize*/)
{
  // Apply the cut. If the list hasn't changed, there's nothing to do
  if( !ds.ApplyCut( W2Min, Q2Min, bStrictCut, SampleSize, bShowCut, DataNodeSortFile ) )
  {
    LOG( Sometimes, "Skipping unchanged cut." << endl );
    return true;
  }
  
  if( !ds.m_TargetsNeeded )
  {
    LOG( Always, "Error - no targets required for prediction." << endl );
    return false;
  }
  
  // Load predictions from file in preference to creating them (to save time)
  bool bLoadError = false;
  int iReplica = ds.m_iNumReplicas;
  std::stringstream sFileName;
  sFileName << pszOutFilePrefix << '_' << ds.Name() << '_' << iSeq << "_Cassandra.F2";
  {
    std::ifstream is(sFileName.str().c_str());
    if( is.fail() )
    { LOG( Often, "Cached " << ds.Name() << " data in " << sFileName.str() << " unavailable" << endl ); }
    else
    {
      LOG( Often, "Loading cached " << ds.Name() << " data from " << sFileName.str() << endl );
      DataSet * pDS = ds.NewDataSet();
      if( pDS->Load(is, iReplica, pPDFSet, pMassScheme, ds.m_iNumReplicas,
                    static_cast<int>( ds.CutSize() ), W2Min, Q2Min, bStrictCut) )
      {
        // It's possible that the cached data is of no use
        // Treat this as an error so as not to overwrite the cache with less data
        if( iReplica >= ds.m_iNumReplicas )
        {
          bLoadError = true;
          LOG( Always, "Cache contains data for replicas >= " << iReplica
              << ", but I only need " << ds.m_iNumReplicas
              << ".\nAborting so as to ovoid overwriting cache." << endl );
        }
        else
        {
          // Copy predictions from the other DataSet
          int i = 0;
          for( list<DataNode *>::iterator it = pDS->m_list.begin()
              ; !bLoadError && it != pDS->m_list.end() ; it++, i++ )
          {
            DataNode * p = *it;
            // First, check that the data match
            if( p->m_x != ds.m_v[i]->m_x )
            {
              bLoadError = true;
              LOG( Always, "Cache " << i << " mismatch on x: " << p->m_x << "!="
                  << ds.m_v[i]->m_x << endl );
            }
            if( p->m_Q2 != ds.m_v[i]->m_Q2 )
            {
              bLoadError = true;
              LOG( Always, "Cache " << i << " mismatch on Q2: " << p->m_Q2 << "!="
                  << ds.m_v[i]->m_Q2 << endl );
            }
            if( p->m_eBeam != ds.m_v[i]->m_eBeam )
            {
              bLoadError = true;
              LOG( Always, "Cache " << i << " mismatch on eBeam: " << p->m_eBeam
                  << "!=" << ds.m_v[i]->m_eBeam << endl );
            }
            if( p->m_i != ds.m_v[i]->m_i )
            {
              bLoadError = true;
              LOG( Always, "Cache " << i << " mismatch on i: " << p->m_i << "!="
                  << ds.m_v[i]->m_i << endl );
            }
            if( p->m_Data.size() != ds.m_v[i]->m_Data.size() )
            {
              bLoadError = true;
              LOG( Always, "Cache " << i << " mismatch on data size: " << p->m_Data.size()
                  << "!="<< ds.m_v[i]->m_Data.size() << endl );
            }
            else
            {
              for(std::size_t j = p->m_Data.size(); j--; )
              {
                if( p->m_Data[j] != ds.m_v[i]->m_Data[j] )
                {
                  bLoadError = true;
                  LOG( Always, "Cache " << i << " mismatch on data[" << j << "]: "
                      << p->m_Data[j] << "!="<< ds.m_v[i]->m_Data[j] << endl );
                }
              }
            }
            if( !bLoadError )
            {
              for( int j = iReplica; j < ds.m_iNumReplicas; j++ )
                ds.m_v[i]->m_Theory[j] = p->m_Theory[j];
            }
          }
        }
      }
      is.close();
    }
  }
  if( bLoadError )
  {
    LOG( Always, "Error - no targets required for prediction." << endl );
    return false;
  }
  
  // References for Q0 = 1.65 GeV for NNPDF 3.1
  // ArXiv 1706.00428, pg 20
  // ArXiv 1605.06515, equation 2, pg 7
  double Q0 = 1.65;
  if( Q0 < ds.m_QMin )
    ds.m_QMin = Q0;
  if( Q0 > ds.m_QMax )
    ds.m_QMax = Q0;
  
  // No need to initialise APFEL if we don't need to run predictions
  if( iReplica )
  {
    IF_LOG(Mostly)
    ;
    else
      APFEL::EnableWelcomeMessage(false);
    APFEL::SetPDFSet(pPDFSet);
    
    // NNPDF sets are for five quark flavours
    int iNumQuarks = 5;
    APFEL::SetMaxFlavourPDFs(iNumQuarks);  // # flavours in DGLAP evolution
    APFEL::SetMaxFlavourAlpha(iNumQuarks); // # flavours in coupling equations
    //APFEL::SetMassScheme("FONLL-B");
    //APFEL::SetMassScheme("FONLL-C"); // Per Prof Ball 2018-06-17
    APFEL::SetMassScheme(pMassScheme);
    //APFEL::SetProcessDIS("EM");
    APFEL::SetQLimits(ds.m_QMin, ds.m_QMax);
    //APFEL::SetPolarizationDIS(0);  // Per Prof Ball 2018-06-08
    //APFEL::SetProjectileDIS("electron");
    //APFEL::SetTargetDIS("proton");
    APFEL::EnableTargetMassCorrections(true); // Per Prof Ball 2018-06-08
    APFEL::EnableDampingFONLL(false); // Per Prof Ball 2018-06-08
    //APFEL::SetFastEvolution(true);
    //APFEL::LockGrids(true);
    //APFEL::EnableEvolutionOperator(true);
    //APFEL::SetFFNS(3);
    //APFEL::SetTheory("QavDS");
    //APFEL::SetTheory("QED");
    //APFEL::SetTheory("QUniD");
    APFEL::EnableLeptonEvolution(true);
    //APFEL::SetTauMass(1e10);
    //  APFEL::SetPerturbativeOrder(2);  // NNLO
    //APFEL::SetPerturbativeOrder(0);
    //APFEL::SetPDFEvolution("exactalpha");
    APFEL::SetPDFEvolution("truncated"); // Per Prof Ball 2018-06-08
    //APFEL::SetPDFSet("NNPDF23_nlo_as_0119_qed.LHgrid");
    //APFEL::SetPDFSet("MRST2004qed.LHgrid");
    //APFEL::SetNumberOfGrids(1);
    //APFEL::SetGridParameters(1,30,3,1e-5);
    //APFEL::SetGridParameters(2,30,3,2e-1);
    //APFEL::SetGridParameters(3,30,3,8e-1);
    //APFEL::SetPDFSet("NNPDF30_nnlo_as_0118.LHgrid");
    //APFEL::SetAlphaQCDRef(0.118,91.2);
    //APFEL::SetAlphaEvolution("expanded");
    //APFEL::SetPDFEvolution("expandalpha");
    //APFEL::SetPoleMasses(1.275,4.18,173.03);
    //APFEL::SetMaxFlavourPDFs(5);
    //APFEL::SetMaxFlavourAlpha(5);
    APFEL::EnableIntrinsicCharm(true); // Per Prof Ball 2018-06-08
    APFEL::InitializeAPFEL_DIS();
  }

  // Perform predictions for requested PDF replicas as requested
  // Backwards so we fail early if there aren't enough replicas in the PDF
  while( iReplica-- )
  {
    APFEL::SetReplica( iReplica );
    // Perform predictions for all the targets required
    for( DataNode::Target t = DataNode::Target::Proton;
        t <= DataNode::Target::Deuteron;
        t = static_cast<DataNode::Target> (t << 1 ) )
    {
      if( t & ds.m_TargetsNeeded )
      {
        switch(t)
        {
          case DataNode::Target::Proton:
            cerr << "DataNode::Target::Proton" << std::endl;
            APFEL::SetTargetDIS("proton");
            break;
          case DataNode::Target::Neutron:
            cerr << "DataNode::Target::Neutron" << std::endl;
            APFEL::SetTargetDIS("neutron");
            break;
          case DataNode::Target::Deuteron:
            cerr << "DataNode::Target::Deuteron" << std::endl;
            APFEL::SetTargetDIS("isoscalar");
            break;
          case DataNode::Target::None:
            cerr << "DataNode::Target::None" << std::endl;
            break;
        }
        double LastQ = 0;
        bool   bFirst = true;
        for( DataNode * n : ds.m_l)
        {
          // If the node is interested in this target
          if( t & n->PredictionTargets() )
          {
            // As an optimisation, only compute structure functions for different energy
            if( bFirst || n->m_cQ != LastQ )
            {
              bFirst = false;
              LastQ = n->m_cQ;
              APFEL::ComputeStructureFunctionsAPFEL(Q0,LastQ);
            }
            n->MakePrediction(t, iReplica);
          }
        }
      }
    }
    
    // Now all the targets are done, we can perform any final calculations
    for( DataNode * n : ds.m_l)
      n->MakePrediction(DataNode::Target::None, iReplica);
    
    // Write this to disk in case we need to restart
    if( !ds.WritePrediction( pPDFSet, pMassScheme, sFileName.str(), iSeq,
                            W2Min, Q2Min, bStrictCut, iReplica) )
    {
      return false;
    }
  }
  return true;
}

bool Cassandra::DataManager::LoadData( const char * fileName )
{
  bool bRet = false;
  
  ifstream  fin(fileName);
  if(!fin.is_open())
  {LOG( Always, "Unable to open " << fileName << endl );}
  else
  {
    LOG( Sometimes, "Parsing data file " << fileName << endl );
    string line;
    Cassandra::DataManager::NodeType nodeType = Cassandra::DataManager::Unknown;
    bRet = true;
    
    // Read the next line from the file
    while( bRet && getline(fin,line) )
    {
      // Ignore blank lines
      if( line.find_first_not_of( " \t") != string::npos )
      {
        std::stringstream ss(line);
        
        // Try to identify the file type from the header
        if( nodeType == Cassandra::DataManager::NodeType::Unknown )
        {
          LOG( Rarely, line << endl ); // Print the header of each file until identified
          string Tokens[7];
          for( string &s : Tokens )
            ss >> s;
          for( nodeType = (Cassandra::DataManager::NodeType) 0;
              nodeType != Cassandra::DataManager::NodeType::Unknown
              && !m_DataSets[nodeType]->MyFileType( Tokens );
              ++(*(( int * ) &nodeType)) ) // a bit naughty, but it's defined to be an int
          {}
          if( nodeType != Cassandra::DataManager::NodeType::Unknown )
          {
            m_DataSets[nodeType]->NewFile();
          }
        }
        else
        {
          // Parse the next line from the file, print it if requested
          if( m_DataSets[nodeType]->ParseLine( ss, bRet ) )
            LOG( Rarely, line << endl );
        }
      }
    }
    
    //Close the file
    fin.close();
    
    if( nodeType != Cassandra::DataManager::NodeType::Unknown )
    {
      m_DataSets[nodeType]->ComputeOverallStats(Cassandra::LogLevel::Rarely, 10);
      LOG( Always, "Loaded " << m_DataSets[nodeType]->m_count << " records from " << fileName << endl );
    }
    else
    {
      LOG( Always, "Unknown file format: " << fileName << endl );
      bRet = false;
    }
  }
  return bRet;
}

void Cassandra::DataManager::ApplyCutoff( double W2, double Q2, bool bStrictCut,
                                         size_t SampleSize, bool bShowCut )
{
  // Apply absolute cut-off and dump lists again if we lost any records
  for( DataSet * pDataSet : m_DataSets )
  {
    if( pDataSet->size() )
    {
      const char * pszCut = Cassandra_StrictString( bStrictCut );
      LOG( Always, pDataSet->Name() << " applying upper cut-off at: W^2" << pszCut << W2
          << ", Q^2" << pszCut << Q2 << endl );
      if( pDataSet->ApplyCut( W2, Q2, bStrictCut, SampleSize, bShowCut, DataNodeSortFile ) )
        pDataSet->ParkCut( SampleSize );
    }
  }
}

void Cassandra::DataManager::ComputeOverallStats(size_t SampleSize)
{
  // Load complete - recompute statistics
  for( DataSet * pDataSet : m_DataSets )
  {
    pDataSet->m_list.sort( DataNodeSortFile );
    pDataSet->ComputeOverallStats(Cassandra::LogLevel::Mostly, SampleSize);
  }
}

bool Cassandra::DataManager::MakePrediction(const char * pPDFSet, const char * pMassScheme,
                                            const char * pszOutFilePrefix,
                                            bool bStrictFirstCut,
                                            int SampleSize, bool bParseOnly, bool bShowCut,
                                            const double PW2Min, const double PW2Max,
                                            const double PQ2Min, const double PQ2Max, const int PNumSteps,
                                            double xBinSize, double Q2BinSize)
{
  bool bMadePrediction = false;
  bool bNoError = true;
  
  // Apply the cut and dump lists again if we lost any records
  for( DataSet * pDataSet : m_DataSets )
  {
    if( pDataSet->size() )
    {
      // We're in a loop, so don't overwrite parameters for other loops
      double W2Min = PW2Min;
      double W2Max = PW2Max;
      double Q2Min = PQ2Min;
      double Q2Max = PQ2Max;
      int NumSteps = PNumSteps;
      
      // Work out where the cut(s) should be based on input parameters & each DataSet
      double W2, W2Step;
      double Q2, Q2Step;
      bool bStrictCut = bStrictFirstCut;
      
      if( W2Min == -1. )
      {
        W2Min = pDataSet->m_W2MinOverall;
        LOG( Mostly, "Taking minimum W2 from dataset = " << setprecision(7) << W2Min << endl );
      }
      else
      { LOG( Sometimes, "Minimum W2 specified >= " << setprecision(7) << W2Min << endl ); }
      
      if( Q2Min == -1. )
      {
        Q2Min = pDataSet->m_Q2MinOverall;
        LOG( Mostly, "Taking minimum Q2 from dataset = " << setprecision(7) << Q2Min << endl );
      }
      else
        LOG( Sometimes, "Minimum Q2 specified >= " << setprecision(7) << Q2Min << endl );
      
      if( NumSteps == 0 )
      {
        if( PW2Min == -1 || PQ2Min == -1 )
          bStrictCut = false;
        W2 = W2Min;
        W2Step = 0.;
        Q2 = Q2Min;
        Q2Step = 0.;
      }
      else
      {
        if( W2Max == -1. )
        {
          W2Max = pDataSet->m_W2MaxOverall;
          LOG( Mostly, "Taking maximum W2 from dataset > " << setprecision(7) << W2Max << endl );
          bStrictCut = true;
        }
        else
        {
          LOG( Sometimes, "Maximum W2 specified" << Cassandra_StrictString( bStrictCut )
              << W2Max << endl );
        }
        
        if( Q2Max == -1. )
        {
          Q2Max = pDataSet->m_Q2MaxOverall;
          LOG( Mostly, "Taking maximum Q2 from dataset > " << setprecision(7) << Q2Max << endl );
          bStrictCut = true;
        }
        else
        {
          LOG( Sometimes, "Maximum Q2 specified" << Cassandra_StrictString( bStrictCut )
              << Q2Max << endl );
        }
        
        W2 = W2Max;
        W2Step = (double) (((long double)(W2Min - W2Max)) / ((long double)NumSteps));
        Q2 = Q2Max;
        Q2Step = (double) (((long double)(Q2Min - Q2Max)) / ((long double)NumSteps));
      }
      
      // Step though each cut
      bool bOK = true; // Stop processing that data type if we get an error
      for( int i = 0 ; bOK && i <= NumSteps ; i++ )
      {
        if( i == 0 )
        {
          // Give the dataset an opportunity to initialise
          if( !bParseOnly && !InitPrediction( * pDataSet, pszOutFilePrefix ) )
            return false;
          LOG( Always, "================================================" << endl );
        }
        else
        {
          LOG( Always, "==========--------------------------------------" << endl );
          // If we're on our last step, adjust cut to avoid rounding errors
          if( NumSteps > 0 && i == NumSteps )
          {
            if( W2Step != 0. )
            {
              LOG( Sometimes, "Adusting W^2 cut from "
                  << setprecision(std::numeric_limits<double>::digits10)
                  << W2 << " to " << W2Min << endl );
              W2 = W2Min;
            }
            if( Q2Step != 0 )
            {
              LOG( Sometimes, "Adusting Q^2 cut from "
                  << setprecision(std::numeric_limits<double>::digits10)
                  << Q2 << " to " << Q2Min << endl );
              Q2 = Q2Min;
            }
          }
        }
        // Log information about this prediction
        const char * pszCut = Cassandra_StrictString( bStrictCut );
        LOG( Always, pDataSet->Name() << " " << i << ": W^2" << pszCut
            << setprecision(7) << W2 << setprecision(5) << ", Q^2" << pszCut << Q2 << endl );
        // Make the prediction
        if( bParseOnly )
        {
          if( pDataSet->ApplyCut( W2, Q2, bStrictCut, SampleSize, bShowCut, DataNodeSortFile ) )
            bMadePrediction = true;
        }
        else
        {
          bOK = MakeOnePrediction( * pDataSet, pPDFSet, pMassScheme, pszOutFilePrefix,
                                  i, W2, Q2, bStrictCut, SampleSize, bShowCut,
                                  xBinSize, Q2BinSize );
          if( bOK )
            bMadePrediction = true;
          else
            bNoError = false;
        }
        bStrictCut = false;
        W2 += W2Step;
        Q2 += Q2Step;
      }
    }
  }
  return bMadePrediction && bNoError;
}

int DISTest(const char * pPDFSet, const char * pMassScheme, int iReplica)
{
  LOG( Always, "Testing APFEL DIS with PDF Set \"" << pPDFSet << "\", Mass Scheme \"" << pMassScheme
      << "\", PDF Replica " << iReplica << endl );
  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2,
    1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  
  // Activate some options
  
  //
  // Settings
  //
  string proc = "NC";
  //APFEL::SetMassScheme("ZM-VFNS");
  APFEL::SetMassScheme(pMassScheme);
  APFEL::SetProcessDIS(proc);
  APFEL::SetQLimits(1.4,250);
  //APFEL::SetPolarizationDIS(0);
  //APFEL::SetProjectileDIS("electron");
  APFEL::SetTargetDIS("proton");
  //APFEL::EnableTargetMassCorrections(false);
  //APFEL::EnableDampingFONLL(true);
  //APFEL::SetFastEvolution(true);
  //APFEL::LockGrids(true);
  //APFEL::EnableEvolutionOperator(true);
  //APFEL::SetFFNS(3);
  //APFEL::SetTheory("QavDS");
  //APFEL::SetTheory("QED");
  //APFEL::SetTheory("QUniD");
  //APFEL::EnableLeptonEvolution(true);
  //APFEL::SetTauMass(1e10);
  //APFEL::SetPerturbativeOrder(0);
  //APFEL::SetPDFEvolution("exactalpha");
  APFEL::SetPDFSet(pPDFSet);
  //APFEL::SetPDFSet("MRST2004qed.LHgrid");
  //APFEL::SetNumberOfGrids(1);
  //APFEL::SetGridParameters(1,30,3,1e-5);
  //APFEL::SetGridParameters(2,30,3,2e-1);
  //APFEL::SetGridParameters(3,30,3,8e-1);
  //APFEL::SetPDFSet("NNPDF30_nnlo_as_0118.LHgrid");
  //APFEL::SetAlphaQCDRef(0.118,91.2);
  //APFEL::SetAlphaEvolution("expanded");
  //APFEL::SetPDFEvolution("expandalpha");
  //APFEL::SetPoleMasses(1.275,4.18,173.03);
  //APFEL::SetMaxFlavourPDFs(5);
  //APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetReplica( iReplica );
  
  // Initializes integrals on the grids
  APFEL::InitializeAPFEL_DIS();
  
  double Q02=2, Q2=100, eps = 1e-10;
  cout << "Initial scale = " << Q02 << " GeV^2. Final scale = " << Q2 << " GeV^2" << endl;
  
  // Load evolution
  double Q0 = sqrt(Q02) - eps;
  double Q  = sqrt(Q2);
  APFEL::ComputeStructureFunctionsAPFEL(Q0,Q);
  
  cout << scientific << setprecision(5) << endl;
  // Tabulate PDFs for the LHA x values
  cout << "alpha_QCD(mu2F) = " << APFEL::AlphaQCD(Q) << endl;
  cout << "alpha_QED(mu2F) = " << APFEL::AlphaQED(Q) << endl;
  cout << endl;
  
  cout << "   x   "
  << setw(11) << "   F2light   "
  << setw(11) << "   F2charm   "
  << setw(11) << "   F2bottom  "
  << setw(11) << "   F2total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4)
    << setw(11) << APFEL::F2light(xlha[i])  << "  "
    << setw(11) << APFEL::F2charm(xlha[i])  << "  "
    << setw(11) << APFEL::F2bottom(xlha[i]) << "  "
    << setw(11) << APFEL::F2total(xlha[i])	 << endl;
  cout << "      " << endl;
  
  cout << "   x   "
  << setw(11) << "   FLlight   "
  << setw(11) << "   FLcharm   "
  << setw(11) << "   FLbottom  "
  << setw(11) << "   FLtotal   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4)
    << setw(11) << APFEL::FLlight(xlha[i])  << "  "
    << setw(11) << APFEL::FLcharm(xlha[i])  << "  "
    << setw(11) << APFEL::FLbottom(xlha[i]) << "  "
    << setw(11) << APFEL::FLtotal(xlha[i])  << endl;
  cout << "      " << endl;
  
  cout << "   x   "
  << setw(11) << "   F3light   "
  << setw(11) << "   F3charm   "
  << setw(11) << "   F3bottom  "
  << setw(11) << "   F3total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4)
    << setw(11) << APFEL::F3light(xlha[i])  << "  "
    << setw(11) << APFEL::F3charm(xlha[i])  << "  "
    << setw(11) << APFEL::F3bottom(xlha[i]) << "  "
    << setw(11) << APFEL::F3total(xlha[i])	 << endl;
  cout << "      " << endl;
  //
  // Cache Structure functions
  //
  APFEL::CacheStructureFunctionsAPFEL(Q0);
  cout << "   x   "
  << setw(11) << "   F2light   "
  << setw(11) << "   F2charm   "
  << setw(11) << "   F2bottom  "
  << setw(11) << "   F2total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4)
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","light",xlha[i],Q)  << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","charm",xlha[i],Q)  << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","bottom",xlha[i],Q) << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F2","total",xlha[i],Q)  << endl;
  cout << "      " << endl;
  
  cout << "   x   "
  << setw(11) << "   FLlight   "
  << setw(11) << "   FLcharm   "
  << setw(11) << "   FLbottom  "
  << setw(11) << "   FLtotal   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4)
    << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","light",xlha[i],Q)  << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","charm",xlha[i],Q)  << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","bottom",xlha[i],Q) << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"FL","total",xlha[i],Q)  << endl;
  cout << "      " << endl;
  
  cout << "   x   "
  << setw(11) << "   F3light   "
  << setw(11) << "   F3charm   "
  << setw(11) << "   F3bottom  "
  << setw(11) << "   F3total   " << endl;
  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4)
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","light",xlha[i],Q)  << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","charm",xlha[i],Q)  << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","bottom",xlha[i],Q) << "  "
    << setw(11) << APFEL::StructureFunctionxQ(proc,"F3","total",xlha[i],Q)  << endl;
  cout << "      " << endl;
  
  return 0;
}

int SimpleTest(const char * pszTitle)
{
  if( pszTitle == nullptr || ! * pszTitle )
    pszTitle = "Use \"-zTitle_of_chart\"";
  TCanvas *c1 = new TCanvas("c1",pszTitle,200,10,700,500);
  
  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);
  
  const Int_t n = 10;
  Float_t x[n]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
  Float_t y[n]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
  //Float_t ex[n] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
  Float_t ey[n] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
  TGraphErrors *gr = new TGraphErrors(n,x,y,nullptr,ey);
  gr->SetTitle(pszTitle);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("ALP");
  
  c1->Update();
  c1->Print("test.pdf");
  /*gPad->Draw();
   cout << pszTitle << endl;
   std::string s;
   cin >> s;
   cout << s << endl;*/
  return 9;
}

int main(int argc, char * argv[])
{
  auto start_time = std::chrono::system_clock::now();
  
  const char * pszOutFilePrefix = "a_out";
  
  // Give instructions if no filename supplied
  if( argc < 2 )
  {
    cout << "usage: " << argv[0] << " [switches] <list-of-filenames-to-process>" << endl
    << "  -bx[n.n]  Bin size in X    (default=0.05)" << endl
    << "  -bq[n.n]  Bin size in Q**2 (default=0.25)" << endl
    << "  -c        Show cut records" << endl
    << "  -i        Input only, i.e. parse, but no predictions" << endl
    << "  -l[n]     Log level n, 0 = terse, 4 = detailed, default 1" << endl
    << "  -l[prefix]Log output to prefix.log, predictions to files starting with prefix" << endl
    << "  -m[mass]  Set mass scheme, default FONLL-C" << endl
    << "  -n[i]     number of intervals between cuts (default 0)" << endl
    << "  -pPDF     LHAPDF name of the PDF to apply, eg: -pNNPDF31_nnlo_as_0118_luxqed.LHgrid" << endl
    << "  -qc[n.n]  Throw away all data above Q^2 cut-off GeV^2 (default=3.5)" << endl
    << "  -ql[n.n]  lower limit for Q^2 (default=data min)" << endl
    << "  -qu[n.n]  Upper limit for Q^2 (default=data max)" << endl
    << "  -r[n]     How many PDF Replicas to make predictions for (0 to n-1)" << endl
    << "  -s[n]     Print n samples from each data set" << endl
    << "  -sc       Make the first cut a strict cut, i.e. > rather than >=" << endl
    << "  -t        Run a test of the APFEL::DIS library - don't process any data" << endl
    << "  -wc[n.n]  Throw away all data above W^2 cut-off GeV^2 (default=12.5)" << endl
    << "  -wl[n.n]  lower limit for W^2 (default=data min)" << endl
    << "  -wu[n.n]  Upper limit for W^2 (default=data max)" << endl
    << "  -z[file]  CERN ROOT chart test - output to file" << endl;
    return -1;
  }
  
  // Read command-line switches
  int SampleSize = 10;
  bool bParseOnly = false;
  bool bShowCut = false;
  double W2Max = -1.;
  double W2Min = -1.;
  double Q2Max = Q2_DEFAULT_CUT;
  double Q2Min = Q2_DEFAULT_CUT;
  int NumSteps = 0;
  bool bGotW2Cutoff = false;
  double   W2Cutoff = 12.5;
  bool bGotQ2Cutoff = false;
  double   Q2Cutoff = 3.5;
  
  int iSetPDFSet = 0;
  const char DefaultPDFSet[] = "NNPDF31_nnlo_as_0118.LHgrid";
  const char AlternatePDFSetM[] = "NNPDF31_nlo_as_0118.LHgrid";
  //const char AlternatePDFSetP[] = "NNPDF31_nnlo_as_0118_luxqed.LHgrid";
  const char * AlternatePDFSetP = AlternatePDFSetM;
  const char * pDefaultPDFSet = DefaultPDFSet;
  
  int iSetMassScheme = 0;
  const char DefaultMassScheme[] = "FONLL-C"; // Per Prof Ball 2018-06-17
  const char AlternateMassSchemeP[] = "FONLL-B";
  const char AlternateMassSchemeM[] = "ZM-VFNS";
  const char * pDefaultMassScheme = DefaultMassScheme;
  
  bool bStrictFirstCut = false;
  bool bTestOnly = false;
  int iNumReplicas = 1;
  
  double xBinSize = 0.05;
  double Q2BinSize = 0.25;
  
  int iReturnValue = 0;
  
  ModelParams Params;
  
  for( int i = 1 ; iReturnValue == 0 && i < argc ; i++ )
  {
    if( argv[i][0] == '-' )
    {
      switch( toupper( argv[i][1] ) )
      {
        case 'A':
        {
          int iWhich = toupper( argv[i][2] ) - 'A';
          if( iWhich >= 0 && iWhich <=3 )
          {
            char   c;
            double d[2];
            std::stringstream s( &argv[i][3] );
            s >> d[0] >> c >> d[1];
            if( !s.fail() && c == ',' )
            {
              int iCount;
              s >> c >> iCount;
              if( s.fail() || c != ',' )
                iCount = 10;
              Params[iWhich].FromTo( d[0], d[1], iCount );
            }
          }
        }
          break;
          
        case 'B':
        {
          char c = toupper( argv[i][2] );
          if( c != 0 && ( argv[i][3] == '-' || argv[i][3] == '.'
                         || ( argv[i][3] >= '0' && argv[i][3] <= '9' ) ) )
          {
            double d = atof( &argv[i][3] );
            if( c == 'X' )
            {
              xBinSize = d;
            }
            else if( c == 'Q' )
            {
              Q2BinSize = d;
            }
          }
        }
          break;
          
        case 'C':
          bShowCut = true;
          break;
          
        case 'I':
          bParseOnly = true;
          break;
          
        case 'L':
        {
          char c = argv[i][2];
          if( c == '-' || ( c >= '0' && c <= '9' ) )
          {
            Global_Log_Level = (Cassandra::LogLevel) atoi( &argv[i][2] );
            LOG( Always, "Log level set to " << ((int)Cassandra::Global_Log_Level) << endl );
          }
          else if( c != 0 )
          {
            pszOutFilePrefix = &argv[i][2];
            std::string sLogFile(pszOutFilePrefix);
            sLogFile += ".log";
            Cassandra::Global_Log_File.reset( new std::ofstream( sLogFile.c_str(),
                                                                ios_base::out | ios_base::trunc ) );
            LOG( Rarely, "Logging to file " << sLogFile << std::endl );
          }
        }
          break;
          
        case 'M':
          if( argv[i][2] != 0 )
          {
            pDefaultMassScheme = &argv[i][2];
            iSetMassScheme = 1;
          }
          else
          {
            pDefaultMassScheme = AlternateMassSchemeM;
            iSetMassScheme = -1;
          }
          break;
          
        case 'N':
        {
          if( argv[i][2] >= '0' && argv[i][2] <= '9' )
            NumSteps = atoi( &argv[i][2] );
        }
          break;
          
        case 'P':
          if( argv[i][2] != 0 )
          {
            pDefaultPDFSet = &argv[i][2];
            iSetPDFSet = 1;
          }
          else
          {
            pDefaultPDFSet = AlternatePDFSetP;
            iSetPDFSet = -1;
          }
          break;
          
        case 'Q':
          if( toupper( argv[i][2] ) == 'C' )
          {
            bGotQ2Cutoff = true;
            double d = atof( &argv[i][3] );
            if( d != 0 )
              Q2Cutoff = d;
          }
          else
          {
            if( toupper( argv[i][2] ) == 'L' && ( argv[i][3] == '-' || argv[i][3] == '.'
                                                 || ( argv[i][3] >= '0' && argv[i][3] <= '9' ) ) )
            {
              Q2Min = atof( &argv[i][3] );
            }
            else if( toupper( argv[i][2] ) == 'U' && ( argv[i][3] == '-' || argv[i][3] == '.'
                                                      || ( argv[i][3] >= '0' && argv[i][3] <= '9' ) ) )
            {
              Q2Max = atof( &argv[i][3] );
            }
          }
          break;
          
        case 'R':
        {
          int iNum = atoi( &argv[i][2] );
          if( iNum > 0 )
            iNumReplicas = iNum;
        }
          break;
          
        case 'S':
          if( toupper( argv[i][2] ) == 'C' )
            bStrictFirstCut = true;
          else
            SampleSize = atoi( &argv[i][2] );
          break;
          
        case 'T':
          bTestOnly = true;
          break;
          
        case 'W':
          if( toupper( argv[i][2] ) == 'C' )
          {
            bGotW2Cutoff = true;
            double d = atof( &argv[i][3] );
            if( d != 0 )
              W2Cutoff = d;
          }
          else
          {
            if( toupper( argv[i][2] ) == 'L' && ( argv[i][3] == '-' || argv[i][3] == '.'
                                                 || ( argv[i][3] >= '0' && argv[i][3] <= '9' ) ) )
            {
              W2Min = atof( &argv[i][3] );
            }
            else if( toupper( argv[i][2] ) == 'U' && ( argv[i][3] == '-' || argv[i][3] == '.'
                                                      || ( argv[i][3] >= '0' && argv[i][3] <= '9' ) ) )
            {
              W2Max = atof( &argv[i][3] );
            }
          }
          break;
          
        case 'Z':
          return SimpleTest( &argv[i][2] );
          
        default:
          cerr << "Ignoring unrecognised switch: " << argv[i] << endl;
          break;
      }
      argv[i][0] = 0;
    }
  }
  
  // Make sense of combined PDF evolution parameters
  if( iSetPDFSet == -1 && iSetMassScheme == 0 )
    pDefaultMassScheme = AlternateMassSchemeP;
  else if( iSetMassScheme == -1 && iSetPDFSet == 0 )
    pDefaultPDFSet = AlternatePDFSetM;
  LOG( Always, "Running prediction using PDF set " << pDefaultPDFSet
      << ", mass scheme " << pDefaultMassScheme << endl );
  
  // If we only specify one limit, make sure we use the one specified
  if( NumSteps == 0 )
  {
    if( Q2Min == Q2_DEFAULT_CUT )
      Q2Min = Q2Max;
    if( W2Min == -1. )
      W2Min = W2Max;
  }
  // If we specify W & Q limits, check whether n makes sense
  else if(    W2Min != -1. && W2Max != -1. && W2Min == W2Max
          && Q2Min != -1. && Q2Max != -1. && Q2Min == Q2Max )
  {
    LOG( Always, "Can't take multiple steps across fixed Q^2 & W^2 limits" << endl );
    return -1;
  }
  
  // Now give a summary of what cuts we will make
  if( NumSteps > 0 )
  {
    LOG( Always, NumSteps << " intervals for Q^2 " << fixed << setprecision(5) );
    if( Q2Min != -1. && Q2Max != -1. && Q2Min == Q2Max )
    { LOG( Always, "= " << Q2Min ); }
    else
    {
      LOG( Always, "on [" );
      if( Q2Min == -1. )
      { LOG( Always, "min" ); }
      else
      { LOG( Always, Q2Min ); }
      LOG( Always, ", " );
      if( Q2Max == -1. )
      { LOG( Always, "max" ); }
      else
      { LOG( Always, Q2Max ); }
      LOG( Always, "]" );
    }
    LOG( Always, ", W^2 " );
    if( W2Min != -1. && W2Max != -1. && W2Min == W2Max )
    { LOG( Always, "= " << W2Min ); }
    else
    {
      LOG( Always, "on [" );
      if( W2Min == -1. )
      { LOG( Always, "min" ); }
      else
      { LOG( Always, W2Min ); }
      LOG( Always, ", " );
      if( W2Max == -1. )
      { LOG( Always, "max" ); }
      else
      { LOG( Always, W2Max ); }
      LOG( Always, "]" );
    }
  }
  else
  {
    LOG( Always, "Q^2 >= " );
    if( Q2Min == -1. )
    { LOG( Always, "min" ); }
    else
    { LOG( Always, Q2Min ); }
    LOG( Always, ", W^2 >= " );
    if( W2Min == -1. )
    { LOG( Always, "min" ); }
    else
    { LOG( Always, W2Min ); }
  }
  LOG( Always, " GeV^2." << endl );
  
  if( bTestOnly )
  {
    LOG( Always, "Testing APFEL DIS library." << endl );
    return DISTest(pDefaultPDFSet, pDefaultMassScheme, iNumReplicas - 1);
  }
  
  // Say what bin sizes we are using
  LOG( Always, "X bin size set to " << xBinSize << endl );
  LOG( Always, "Q**2 bin size set to " << Q2BinSize << endl );
  
  // Say what the model is
  LOG( Always, scientific << setprecision(std::numeric_limits<double>::digits10)
      << "A: Start=" << Params[0].Start << ", Step=" << Params[0].Step << " Num=" << Params[0].Num << endl
      << "B: Start=" << Params[1].Start << ", Step=" << Params[1].Step << " Num=" << Params[1].Num << endl
      << "C: Start=" << Params[2].Start << ", Step=" << Params[2].Step << " Num=" << Params[2].Num << endl );
  
  Cassandra::DataManager myPred( Params, iNumReplicas );
  for( int i = 1 && iReturnValue == 0 ; i < argc ; i++ )
  {
    if( argv[i][0] != 0 )
    {
      glob_t  myGlob;
      myGlob.gl_pathc = 0;
      myGlob.gl_pathv = NULL;
      myGlob.gl_offs = 0;
      if( glob( argv[i], 0, NULL, &myGlob ) == 0 )
      {
        for( size_t j = 0; j < myGlob.gl_pathc ; j++ )
          if( !myPred.LoadData( myGlob.gl_pathv[j] ) )
            iReturnValue = -2;
      }
      globfree( &myGlob );
    }
  }
  
  // Make sure load was error free
  if( iReturnValue == 0 )
  {
    if( bGotW2Cutoff || bGotQ2Cutoff )
    {
      // Apply cutoff NB: recomputes statistics
      if( !bGotW2Cutoff ) W2Cutoff = 0;
      if( !bGotQ2Cutoff ) Q2Cutoff = 0;
      myPred.ApplyCutoff( W2Cutoff, Q2Cutoff, bStrictFirstCut, 10, false );
      // Strictness applies only to the cutoff, if specified
      bStrictFirstCut = false;
    }
    else
    {
      // Load complete - recompute statistics
      myPred.ComputeOverallStats( SampleSize );
    }
    if( !myPred.size() )
    {
      LOG( Always, "No data for predictions." << endl );
      iReturnValue = -3;
    }
    else
    {
      // Make predictions
      if( !myPred.MakePrediction( pDefaultPDFSet, pDefaultMassScheme,
                                 pszOutFilePrefix, bStrictFirstCut, SampleSize,
                                 bParseOnly, bShowCut,
                                 W2Min, W2Max, Q2Min, Q2Max, NumSteps,
                                 xBinSize, Q2BinSize ) )
        iReturnValue = -4;
    }
  }
  
  auto end_time = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end_time - start_time;
  LOG( Always, "Total run time " << fixed << setprecision(1) << elapsed_seconds.count() << " seconds." << endl
      << "Return code " << iReturnValue << endl );
  return iReturnValue;
}
