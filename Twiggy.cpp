# include "Cassandra.h"
# include "Twiggy.h"

using namespace std;
using namespace Cassandra;

double Cassandra::Model::RawModel( const DataNode * p, int k ) const
{
  double c_i;
  switch( m_iModelNum )
    {
    case 1:
      // c * x^a * ( 1 - x )^b * ( 1 + x * d )
      c_i = m_Model[k].p[2] * pow( p->m_x, m_Model[k].p[0] ) * pow( 1. - p->m_x, m_Model[k].p[1] )
	    * ( 1. + p->m_x * m_Model[k].p[3] );
      break;

    case 2:
      // c * x^a * ( 1 + x * b ) + d / Q^2
      c_i = m_Model[k].p[2] * pow( p->m_x, m_Model[k].p[0] ) * ( 1. + p->m_x * m_Model[k].p[1] );
      if( m_Model[k].p[3] != 0. )
	c_i += m_Model[k].p[3] / p->m_Q2;
      break;

    default:
      // Original model
      // c * x^a * ( 1 - x )^b + d / Q^2
      c_i = m_Model[k].p[2] * pow( p->m_x, m_Model[k].p[0] ) * pow( 1. - p->m_x, m_Model[k].p[1] );
      if( m_Model[k].p[3] != 0. )
	c_i += m_Model[k].p[3] / p->m_Q2;
      break;
    }
  return c_i;
}

double Cassandra::Model::HigherTwist( const DataNode * p, int k, int iReplica ) const
{
  double dHT = RawModel( p, k ) / p->m_Q2;
  if( !m_bAdditive )
    dHT *= p->m_Theory[iReplica];
  return dHT;
}

const char * Cassandra::DataManager::m_pszChiSuffix = ".chi";
const char Cassandra::DataManager::m_szModelSuffix[] = "_model";
const char Cassandra::DataManager::m_szModelHeader[] = "Seq a b c d w0 Chi0";
const char Cassandra::DataManager::m_szPDFSuffix[] = ".pdfwt";
const char Cassandra::DataManager::m_szPDFHeader[] = "Cut Strict X W**2 Q**2 Ndat Neff";
const char Cassandra::DataManager::m_szOptionStdDev[] = "s";

Cassandra::DataManager::DataManager(ModelParams & mp, int iNumReplicas,
				    NodeSortFunc TwiggySort, int NumOutliers, const char * pszModelPrefix,
				    double XCutoff, int iModelNum, bool bAdditive, const char * pszPDFSet,
				    double dQChartScale, double xBinSize, double Q2BinSize )
  : m_ht( iModelNum, bAdditive ), m_Params{ mp }, m_TwiggySort { TwiggySort }, m_NumOutliers { NumOutliers },
    m_pszModelPrefix{ pszModelPrefix }, m_XCutoff{ XCutoff },
    m_pszPDFSet{ pszPDFSet }, m_dQChartScale{ dQChartScale }, m_xBinSize{ xBinSize }, m_Q2BinSize{ Q2BinSize }
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

double Cassandra::DataManager::SeriesValue( const DataNode &n, int iSeries ) const
{
  double dVal;
  switch( iSeries )
    {
    case 0:
      dVal = n.m_x;
      break;
    case 1:
      dVal = n.m_Q2;
      break;
    default:
      dVal = 0;
    }
  return dVal;
}

double Cassandra::DataManager::SeriesValueBinned( const DataNode &n, int iSeries ) const
{
  double dVal = SeriesValue( n, iSeries );
  switch( iSeries )
    {
    case 0:
      if( dVal < 0.085 )
	dVal = 0.07;
      else if( dVal < 0.12 )
	dVal = 0.1;
      else if( dVal < 0.16 )
	dVal = 0.14;
      else if( dVal < 0.20125 )
	dVal = 0.18;
      else if( dVal < 0.25 )
	dVal = 0.225;
      else if( dVal < 0.3125 )
	dVal = 0.275;
      else
	dVal = ( ( int ) ( ( dVal - 0.35 ) / 0.1 ) ) * 0.1 + 0.35;
      break;
    case 1:
      if( m_Q2BinSize != 0. )
	{
	  double dBinSize = m_Q2BinSize;
	  if( ( int ) ( dVal / dBinSize ) < 5 )
	    dBinSize /=2;
	  dVal = ( ( int ) ( dVal / dBinSize ) ) * dBinSize;
	}
      break;
    }
  return dVal;
}

bool Cassandra::DataManager::InitPrediction( const char * pszOutFilePrefix, int iNumReplicas )
{
  bool bOK = true;
  for( int i = 0 ; bOK && i < 2 ; i++ )
    {
      bool bMakeFile = true;
      string sFileName( pszOutFilePrefix );
      if( i == 0 )
	sFileName.append( m_pszChiSuffix );
      else
	{
	  // Only create PDF weights if we're pre-loading the model and more than one PDF replica
	  if( m_pszModelPrefix && iNumReplicas > 1 )
	    sFileName.append( m_szPDFSuffix );
	  else
	    bMakeFile = false;
	}
      if( bMakeFile )
	{
	  // Create the file
	  std::ofstream os( sFileName, ios_base::out | ios_base::trunc );
	  if( os.fail() )
	    {
	      LOG( Always, "Unable to create " << sFileName << endl );
	      bOK = false;
	    }
	  else
	    {
	      // Write header
	      if( i == 0 )
		os << "Cut Strict X W**2 Q**2 Chi NDat Chi/NDat Chi_HT/NDat Chi_HT_Reweight/NDat Neff" << endl;
	      else
		{
		  os << m_szPDFHeader;
		  // All the replica weights first
		  for( int iRep = 1 ; iRep < iNumReplicas ; iRep++ )
		    os << " w" << iRep;
		  // Followed by the chi squared
		  for( int iRep = 1 ; iRep < iNumReplicas ; iRep++ )
		    os << " c" << iRep;
		  os << endl;
		}
	      os.close();
	    }
	}
    }
  return bOK;
}

Int_t Cassandra::DataManager::ValidateMatrices( const TMatrixDSym &MData, const TMatrixDSym &MTheory,
						const std::vector<DataNode *> & vNodes ) const
{
  Int_t ret_val = 0;
  Int_t MSize = MData.GetNrows();
  static const char szBy[] = " x ";
  if( MSize < 1 )
    {
      LOG( Always, "Error: Data matrix only has " << MSize << " rows" << endl );
    }
  else if( MData.GetNcols() != MSize )
    {
      LOG( Always, "Error: Data matrix not symmetric, " << MSize << " x " << MData.GetNcols() << endl );
    }
  else if( MTheory.GetNrows() != MSize || MTheory.GetNcols() != MSize )
    {
      LOG( Always, "Error: Theory matrix (" << MTheory.GetNrows() << szBy << MTheory.GetNcols()
	   << ") doesn't match data (" << MSize << szBy << MSize << ")" << endl );
    }
  else if( ( Int_t ) vNodes.size() != MSize )
    {
      LOG( Always, "Error: " << MSize << szBy << MSize << " matrices for " << vNodes.size() <<" nodes"<<endl);
    }
  else
    ret_val = MSize;
  return ret_val;
}

void Cassandra::DataManager::PlotHT( const DataSet &ds,
				     const std::string &sFileNameSeq, const std::string &sCut,
				     double * HTMean, const TMatrixDSym & MTheoryCovar ) const
{
  std::string sFileName( sFileNameSeq );
  sFileName.append( "_HT.pdf" );

  std::string sTitle( "Higher Twist " );
  sTitle.append( ds.Name() );
  std::string sTitleX( "X (" );
  sTitleX.append( std::to_string( ds.CutSize() ) );
  sTitleX.append( " data points, " );
  sTitleX.append( sCut );
  sTitleX.append( ")" );

  // Make a list of the x values in the data
  std::vector<double> x_vals;
  for( DataNode * p : ds.m_v )
    {
      size_t i = 0;
      while( i < x_vals.size() && x_vals[i] != p->m_x )
	i++;
      if( i == x_vals.size() )
	x_vals.push_back( p->m_x );
    }
  std::sort(x_vals.begin(), x_vals.end());
  int iNumX = ( int ) x_vals.size();
  double ex[iNumX];

  // Calculate error in each x ... by using half distance to nearest neighbour
  if( iNumX == 1 )
    ex[0] = 0.;
  else
    {
      // First bin
      ex[0] = ( x_vals[1] - x_vals[0] ) / 2;
      // Middle bins
      int i = 1;
      for( ; i < iNumX - 1 ; i++ )
	{
	  double dLeft = ( x_vals[i] - x_vals[i - 1] ) / 2;
	  double dRight = ( x_vals[i + 1] - x_vals[i] ) / 2;
	  ex[i] = ( dLeft < dRight ) ? dLeft : dRight;
	}
      // Last bin
      ex[i] = ( x_vals[i] - x_vals[i - 1] ) / 2;
    }

  // Now calculate the mean of the y-values and the error in the y-values
  int    iCount[iNumX];
  double y[iNumX];
  double ey[iNumX];
  for( int j = 0 ; j < iNumX ; j++ )
    {
      iCount[j] = 0;
      y[j] = 0;
      ey[j] = 0;
    }
  for( int i = 0 ; i < ( int ) ds.CutSize() ; i++ )
    {
      int j = 0;
      while( x_vals[j] != ds.m_v[i]->m_x )
	j++;
      iCount[j]++;
      y[j] += HTMean[i];
      ey[j] += MTheoryCovar[i][i];
    }
  for( int j = 0 ; j < iNumX ; j++ )
    {
      y[j] /= iCount[j];
      ey[j] = sqrt( ey[j] / iCount[j] );
    }

  // Now plot it
  //( "HTStats", sTitle.c_str(), iNumX, x_bins, m_szOptionStdDev );
  TGraphErrors * ge = new TGraphErrors( iNumX, &x_vals[0], y, nullptr, ey );

  auto c1 = new TCanvas( "c1", "", 200, 10, 700, 500 );
  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);

  ge->SetLineColor( kBlue );
  ge->SetFillColor( kRed );
  ge->SetFillStyle( 3010 );
  ge->SetTitle( sTitle.c_str() );
  ge->GetXaxis()->SetTitle( sTitleX.c_str() );
  auto t = ge->GetYaxis();
  t->SetTitle( "Higher Twist / Gev**2" );
  t->SetTitleOffset( 1.5 );
  ge->Draw( "a3 L" );
  //tp->Draw( "E1 X0 P L" );
  c1->Update();
  /*TPaveStats * st = (TPaveStats*) tp->FindObject("stats");
  if( st )
    {
      //double a = st->GetX2NDC() - st->GetX1NDC();
      //st->SetX1NDC( 0.7 * a );
      //st->SetX2NDC( 1.7 * a );
      double y2 = st->GetY2NDC();
      double y1 = st->GetY1NDC();
      double a = 2.2 * ( y2 - y1 );
      st->SetY2NDC( y2 - a );
      st->SetY1NDC( y1 - a );
      }*/
  auto OldErrLvl = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  c1->Print( sFileName.c_str() );
  gErrorIgnoreLevel = OldErrLvl;
  //if( st != nullptr )
    //delete st;
  delete c1;
  delete ge;
}

void Cassandra::DataManager::CreateGraphs( const TMatrixDSym & MData, const TMatrixDSym & MTheory,
					   const std::vector<DataNode *> & vNodes, double /*dNeff*/,
					   const std::string &sNamePrefix, const std::string &sFileNameSeq,
					   const std::string &sCut, const std::string &/*sCounts*/,
					   double * HTMean ) const
{
  Int_t MSize = ValidateMatrices( MData, MTheory, vNodes );
  if( MSize == 0 )
    return;

  // c=0 : x-axis=x,  series=Q2
  // c=1 : x-axis=Q2, series=x
  std::vector<double> Unbinned[2];
  std::vector<double> Binned[2];
  bool bSLACPresent = false;
  for( int c = 0 ; c < 2 ; c++ )
    {
      // Make a list of the x-axis values in the data
      for( DataNode * p : vNodes )
	{
	  // Do the unbinned, x_axis values
	  size_t i = 0;
	  double thisVal;
	  if( p->Type() == pszBCDMS )
	    thisVal = SeriesValue( *p, c );
	  else
	    thisVal = SeriesValueBinned( *p, c );
	  while( i < Unbinned[c].size() && Unbinned[c][i] != thisVal )
	    i++;
	  if( i == Unbinned[c].size() )
	    Unbinned[c].push_back( thisVal );
	  // Do the binned, x_axis values
	  i = 0;
	  thisVal = SeriesValueBinned( *p, c );
	  while( i < Binned[c].size() && Binned[c][i] != thisVal )
	    i++;
	  if( i == Binned[c].size() )
	    Binned[c].push_back( thisVal );
	  // Is there any SLAC data present?
	  if( p->Type() == pszSLAC )
	    bSLACPresent = true;
	}
      std::sort( Unbinned[c].begin(), Unbinned[c].end() );
      std::sort( Binned[c].begin(), Binned[c].end() );
    }

  // Create each chart
  for( int c = 0 ; c < 2 ; c++ )
    {
      std::vector<double> & xb = Binned[c];
      std::vector<double> & xu = Unbinned[c];
      std::vector<double> & Series = Binned[1 - c];
      const int iNumSeries = ( c == 0 && Series.size() > 11 ) ? 11 : ( int ) Series.size();

      std::string sFileName( sFileNameSeq );
      sFileName.append( "_F2DTvs" );
      sFileName.append( ( c == 0 ) ? "X.pdf" : "Q2.pdf" );

      std::string sTitle( "F2 Data/Theory vs " );
      sTitle.append( ( c == 0 ) ? "X, " : "Q**2, " );
      sTitle.append( sNamePrefix );
      sTitle.append( ( c == 0 ) ? "; X (" : "; Q**2 / GeV**2 (" );
      sTitle.append( std::to_string( MSize ) );
      sTitle.append( " data points, " );
      sTitle.append( sCut );
      sTitle.append( "); F2 (b&w=LT, colour=HT)" );

      LOG( Often, "Creating " << sFileName << endl );

      // Create one graph per bin, with and without HT
      auto mgLegend = new TMultiGraph();    // Add graphs that appear in legend here
      auto mgNoLegend = new TMultiGraph();  // Add graphs that do not appear in legend here
      mgLegend->SetTitle( sTitle.c_str() );
      const double dScale[] = { 4.8, 3., 2., 1.5, 1.2 };
      const double dScaleSize = sizeof( dScale ) / sizeof( double );
      for( int iSeries = 0 ; iSeries < iNumSeries ; iSeries++ )
	{
	  double yScale;
	  if( c== 0 )
	    yScale = 1 << (iNumSeries - 1 - iSeries);
	  else
	    yScale = ( iSeries < dScaleSize ) ? dScale[iSeries] : 1.;
	  // Unbinned
	  int    iUCount[xu.size()];
	  double F2D[xu.size()];
	  double F2DVar[xu.size()];
	  for( size_t j = 0 ; j < xu.size() ; j++ )
	    {
	      iUCount[j] = 0;
	      F2D[j]     = 0;
	      F2DVar[j]  = 0;
	    }
	  // Binned
	  int    iBCount[xb.size()];
	  double F2T[xb.size()];
	  double F2THT[xb.size()];
	  double F2THTVar[xb.size()];
	  for( size_t j = 0 ; j < xb.size() ; j++ )
	    {
	      iBCount[j]  = 0;
	      F2T[j]      = 0;
	      F2THT[j]    = 0;
	      F2THTVar[j] = 0;
	    }
	  //LOG( Rarely, "check" << endl );
	  int iSeriesCount = 0;
	  for( int i = 0 ; i < MSize ; i++ )
	    {
	      if( !bSLACPresent || vNodes[i]->Type() != pszJLAB )
		{
		  // Only include data points in this series
		  double dNodeSeries = SeriesValueBinned( * vNodes[i], 1 - c );
		  if( dNodeSeries == Series[iSeries]
		     || (iSeries == iNumSeries - 1 && dNodeSeries > Series[iSeries]) )
		    {
		      iSeriesCount++;
		      // Unbinned
		      double dNodeAxis;
		      if( vNodes[i]->Type() == pszBCDMS )
			dNodeAxis = SeriesValue( * vNodes[i], c );
		      else
			dNodeAxis = SeriesValueBinned( * vNodes[i], c );
		      size_t j = 0;
		      while( j < xu.size() && xu[j] != dNodeAxis )
			j++;
		      F2D[j]      += vNodes[i]->ScaledValue();
		      F2DVar[j]   += MData[i][i];
		      iUCount[j]++;
		      // Binned
		      dNodeAxis = SeriesValueBinned( * vNodes[i], c );
		      j = 0;
		      while( j < xb.size() && xb[j] != dNodeAxis )
			j++;
		      double dF2T  = vNodes[i]->m_Theory[0];
		      F2T[j]      += dF2T;
		      F2THT[j]    += dF2T + HTMean[i];
		      F2THTVar[j] += MTheory[i][i];
		      iBCount[j]++;
		    }
		}
	    }
	  //LOG( Rarely, "Series " << iSeries << ": Making graph values" << endl );
	  // Unbinned graphs
	  double g_xu[xu.size()];
	  double g_F2D[xu.size()];
	  double g_F2DErr[xu.size()];
	  int g_UCount = 0;
	  for( size_t j = 0 ; j < xu.size() ; j++ )
	    {
	      if( iUCount[j] > 0 )
		{
		  // Now add this to the graphs, plotting error = sqrt( var )
		  g_xu[g_UCount]     = xu[j];
		  g_F2D[g_UCount]    = F2D[j] / iUCount[j] * yScale;
		  g_F2DErr[g_UCount] = sqrt( F2DVar[j] / iUCount[j] ) * yScale;
		  g_UCount++;
		}
	    }
	  // Binned graphs
	  double g_xb[xb.size()];
	  double g_F2T[xb.size()];
	  double g_F2THT[xb.size()];
	  double g_F2THTErr[xb.size()];
	  int g_BCount = 0;
	  for( size_t j = 0 ; j < xb.size() ; j++ )
	    {
	      if( iBCount[j] > 0 )
		{
		  // Now add this to the graphs, plotting error = sqrt( var )
		  g_xb[g_BCount]       = xb[j];
		  g_F2T[g_BCount]      = F2T[j] / iBCount[j] * yScale;
		  g_F2THT[g_BCount]    = F2THT[j] / iBCount[j] * yScale;
		  g_F2THTErr[g_BCount] = sqrt( F2THTVar[j] / iBCount[j] ) * yScale;
		  g_BCount++;
		}
	    }

	  //Make series name
	  std::stringstream ssTitle;
	  int iSeriesPrecision;
	  const char * pszSeries;
	  if( c == 0 )
	    {
	      static const char szSeriesQ2[] = "Q2";
	      pszSeries = szSeriesQ2;
	      iSeriesPrecision = 1;
	    }
	  else
	    {
	      static const char szSeriesX[] = "x";
	      pszSeries = szSeriesX;
	      iSeriesPrecision = 2;
	    }
	  static const char szCompGE[] = ">=";
	  const char * pszComp = &szCompGE[1];
	  if( iSeries == iNumSeries - 1 && iNumSeries != ( int ) Series.size() )
	    pszComp--;
	  ssTitle << pszSeries << pszComp << fixed << setprecision( iSeriesPrecision ) << Series[iSeries];
	  if( iSeriesPrecision != 1 )
	    ssTitle << setprecision(1);
	  ssTitle << ", x" <<  yScale;

	  // Make the F2T Plot
	  TGraph * g = new TGraph( g_BCount, g_xb, g_F2T );
	  g->SetLineStyle( 7 );
	  mgNoLegend->Add( g, "" );

	  // Make the F2THT Plot
	  //ge = new TGraph( g_Count, g_x, g_F2THT );
	  auto ge = new TGraphErrors( g_BCount, g_xb, g_F2THT, nullptr, g_F2THTErr );
	  int iColor = iSeries % 8 + 2;
	  ge->SetLineColor( iColor );
	  ge->SetMarkerColor( iColor );
	  ge->SetFillColor( 29 );
	  ge->SetFillStyle( 3002 );
	  ge->SetTitle( ssTitle.str().c_str() );
	  mgLegend->Add( ge, "" );

	  // Make the F2D Plot
	  ge = new TGraphErrors( g_UCount, g_xu, g_F2D, nullptr, g_F2DErr );
	  ge->SetMarkerStyle( kOpenCircle );
	  ge->SetMarkerSize( 0.45 );
	  //ge->SetMarkerColor( iColor );
	  ge->SetMarkerColor( kBlack );
	  mgNoLegend->Add( ge, "P" );
	}
      //LOG( Rarely, "Printing graphs" << endl );
      TCanvas *c1 = new TCanvas("c1","Unused title - was ss1",200,10,700,500);
      c1->SetGrid();
      c1->GetFrame()->SetBorderSize(12);
      if( c == 1 )
	c1->SetLogx();
      c1->SetLogy();
      mgLegend->Draw( "A3 L" );
      TLegend * pLeg = c1->BuildLegend( 0.78, 0.71, 0.9, 0.94 );
      mgNoLegend->Draw("L");
      //if( pLeg ) pLeg->SetNColumns( 2 );
      auto OldErrLvl = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kWarning;
      c1->Print(sFileName.c_str());
      gErrorIgnoreLevel = OldErrLvl;
      delete mgNoLegend;
      delete mgLegend;
      if( pLeg )
	delete pLeg;
      delete c1;
    }
}

/*void Cassandra::DataManager::CreateGraphsOld( const TMatrixDSym & MData, const TMatrixDSym & MTheory,
					      const std::vector<DataNode *> & vNodes, double dNeff,
					      const std::string &sNamePrefix, const std::string &sFileNameSeq,
					      const std::string &sCut, const std::string &sCounts,
					      double * HTMean, double xBinSize, double Q2BinSize ) const
{
  Int_t MSize = ValidateMatrices( MData, MTheory, vNodes );
  if( MSize == 0 )
    return;

  // Do this twice, once for X, then for Q**2
  int iNumX[2];
  std::vector<double> x_axis[2];
  std::vector<double> x_axis_error[2];
  for( int c = 0 ; c < 2 ; c++ )
    {
      // Make a list of the x-axis values in the data
      for( DataNode * p : vNodes )
	{
	  int i = 0;
	  double thisVal;
	  if( c == 0 )
	    thisVal = BinnedValue( p->m_x, xBinSize );
	  else
	    thisVal = BinnedValue( p->m_Q2, Q2BinSize );
	  while( i < x_axis[c].size() && x_axis[c][i] != thisVal )
	    i++;
	  if( i == x_axis[c].size() )
	    x_axis[c].push_back( thisVal );
	}
      std::sort(x_axis[c].begin(), x_axis[c].end());

      // Calculate error in each x ... by using half distance to nearest neighbour
      iNumX[c] = x_axis[c].size();
      x_axis_error[c].resize( iNumX[c] );
      if( iNumX[c] == 1 )
	x_axis_error[c][0] = 0.;
      else
	{
	  // First bin
	  x_axis_error[c][0] = ( x_axis[c][1] - x_axis[c][0] ) / 2;
	  // Middle bins
	  int i = 1;
	  for( ; i < iNumX[c] - 1 ; i++ )
	    {
	      double dLeft = ( x_axis[c][i] - x_axis[c][i - 1] ) / 2;
	      double dRight = ( x_axis[c][i + 1] - x_axis[c][i] ) / 2;
	      x_axis_error[c][i] = ( dLeft < dRight ) ? dLeft : dRight;
	    }
	  // Last bin
	  x_axis_error[c][i] = ( x_axis[c][i] - x_axis[c][i - 1] ) / 2;
	}
    }

  for( int c = 0 ; c < 2 ; c++ )
    {
      std::vector<double> &x = x_axis[c];
      std::vector<double> &ex = x_axis_error[c];
      std::vector<double> &Series = x_axis[1 - c];
      // Limit the number of Q**2 series
      int iMaxSeries = Series.size();
      if( c == 0 && iMaxSeries > 11 )
	iMaxSeries = 11;

      std::string sFileName( sFileNameSeq );
      sFileName.append( "_F2DTvs" );
      sFileName.append( ( c == 0 ) ? "X.pdf" : "Q2.pdf" );

      std::string sTitle( "F2 D/T vs " );
      sTitle.append( ( c == 0 ) ? "X, " : "Q**2, " );
      sTitle.append( sNamePrefix );
      sTitle.append( ( c == 0 ) ? "; X (" : "; Q**2 / GeV**2 (" );
      sTitle.append( std::to_string( MSize ) );
      sTitle.append( " data points, " );
      sTitle.append( sCut );
      sTitle.append( "); F2 Data/Theory" );

      LOG( Often, "Creating " << sFileName << endl );

      // Create one graph per bin, with and without HT
      TGraphErrors * gF2DT[ Series.size() ];
      TGraphErrors * gF2DTHT[ Series.size() ];
      auto mgF2DT = new TMultiGraph();
      mgF2DT->SetTitle( sTitle.c_str() );
      double yScale = ( iMaxSeries + 1 ) * 0.5;
      for( int iSeries = iMaxSeries - 1 ; iSeries >= 0 ; iSeries--, yScale -= 0.5 )
	{
	  //LOG( Rarely, "Series " << iSeries << ": Making means and error for " \
	       << MSize << " data points using " << x.size() << " bins" <<endl);
	  // Now calculate the mean and error in F2, with and without the twist
	  int    iCount[x.size()];
	  double F2DT[x.size()];
	  double F2DTVar[x.size()];
	  double F2DTHT[x.size()];
	  double F2DTHTVar[x.size()];
	  for( int j = 0 ; j < x.size() ; j++ )
	    {
	      iCount[j]    = 0;
	      F2DT[j]      = 0;
	      F2DTVar[j]   = 0;
	      F2DTHT[j]    = 0;
	      F2DTHTVar[j] = 0;
	    }
	  //LOG( Rarely, "check" << endl );
	  int iSeriesCount = 0;
	  for( int i = 0 ; i < MSize ; i++ )
	    {
	      // Only include data points in this series
	      double dNodeSeries;
	      if( c == 1 )
		dNodeSeries = BinnedValue( vNodes[i]->m_x, xBinSize );
	      else
		dNodeSeries = BinnedValue( vNodes[i]->m_Q2, Q2BinSize );
	      //LOG( Rarely, iSeries << ": data point " << i << " is Series " << dNodeSeries << endl );
	      if( dNodeSeries == Series[iSeries] || iSeries == iMaxSeries-1 && dNodeSeries > Series[iSeries] )
		{
		  iSeriesCount++;
		  // Find bin for this x-axis value
		  double dNodeAxis;
		  if( c == 0 )
		    dNodeAxis = BinnedValue( vNodes[i]->m_x, xBinSize );
		  else
		    dNodeAxis = BinnedValue( vNodes[i]->m_Q2, Q2BinSize );
		  int j = 0;
		  while( j < x.size() && x[j] != dNodeAxis )
		    j++;
		  if( j == x.size() )
		    {
		      LOG( Always, iSeries << ": data point " << i << " unknown bin" << endl );
		      j = 0;
		    }
		  else
		    {
		      //LOG( Rarely, iSeries << ": data point " << i << " is bin " << j << endl );
		    }
		  // Raw values
		  double dF2D      = vNodes[i]->ScaledValue();
		  double dF2DVar   = MData[i][i];
		  double dF2T      = vNodes[i]->m_Theory[0];
		  double dF2THT    = dF2T + HTMean[i];
		  double dF2THTVar = MTheory[i][i];
		  // Calculate values/variance for current data point
		  double dF2DT      = dF2D / dF2T;
		  double dF2DTVar   = dF2DVar / ( dF2T * dF2T );
		  double dF2DTHT    = dF2D / dF2THT;
		  double dF2DTHTVAR = ( dF2DVar / ( dF2D * dF2D ) + dF2THTVar / ( dF2THT * dF2THT ) )
		                    * ( dF2DTHT * dF2DTHT );
		  // Now add these to statistics for bin
		  iCount[j]++;
		  F2DT[j]      += dF2DT;
		  F2DTVar[j]   += dF2DTVar;
		  F2DTHT[j]    += dF2DTHT;
		  F2DTHTVar[j] += dF2DTHTVAR;
		}
	    }
	  //LOG( Rarely, "Series " << iSeries << ": Making graph values" << endl );
	  // Only graph bins with values
	  double g_x[x.size()];
	  double g_xErr[x.size()];
	  double g_F2DT[x.size()];
	  double g_F2DTErr[x.size()];
	  double g_F2DTHT[x.size()];
	  double g_F2DTHTErr[x.size()];
	  int g_Count = 0;
	  for( int j = 0 ; j < x.size() ; j++ )
	    {
	      if( iCount[j] > 0 )
		{
		  // Normalise the means
		  F2DT[j]     /= iCount[j];
		  F2DTHT[j]   /= iCount[j];
		  // Normalise the variances, then square root to obtain error
		  F2DTVar[j]   /= iCount[j];
		  F2DTHTVar[j] /= iCount[j];
		  // Now add this to the graphs, plotting error = sqrt( var )
		  g_x[g_Count] = x[j];
		  g_xErr[g_Count] = ex[j];
		  g_F2DT[g_Count] = F2DT[j] * yScale;
		  g_F2DTErr[g_Count] = sqrt( F2DTVar[j] );
		  g_F2DTHT[g_Count] = F2DTHT[j] * yScale;
		  g_F2DTHTErr[g_Count] = sqrt( F2DTHTVar[j] );
		  g_Count++;
		}
	    }

	  // Now make plots
	  //LOG( Rarely, "Series " << iSeries << ": Making graphs" << endl );
	  static const char szDrawSingle[] = "APZ";
	  static const char szDrawMulti[] = "";
	  const char * pszOption;
	  bool bSingleton = ( g_Count < 2 );

	  // Make the F2D/T Plot
	  gF2DT[iSeries] = new TGraphErrors( g_Count, g_x, g_F2DT, nullptr, g_F2DTErr );
	  int iColor = iSeries % 8 + 2;
	  if( bSingleton )
	    {
	      pszOption = szDrawSingle;
	      gF2DT[iSeries]->SetMarkerStyle( 34 );
	      gF2DT[iSeries]->SetMarkerSize( 0.75 );
	      gF2DT[iSeries]->SetMarkerColor( kBlack );
	    }
	  else
	    {
	      pszOption = szDrawMulti;
	      gF2DT[iSeries]->SetLineStyle( 2 );
	      gF2DT[iSeries]->SetLineColor( kBlack );
	      gF2DT[iSeries]->SetFillColor( ( iSeries & 1 ) * 2 + 2 );
	      gF2DT[iSeries]->SetFillStyle( 3002 );
	    }
	  std::stringstream ssTitle;
	  int iSeriesPrecision;
	  static const char szSeriesQ2[] = "Q2";
	  static const char szSeriesX[] = "x";
	  const char * pszSeries;
	  if( c == 0 )
	    {
	      //ssTitle << "Q2=" << fixed << setprecision(1) << Series[iSeries];
	      iSeriesPrecision = 1;
	      pszSeries = szSeriesQ2;
	    }
	  else
	    {
	      iSeriesPrecision = 2;
	      pszSeries = szSeriesX;
	    }
	  static const char szCompGE[] = ">=";
	  const char * pszComp = &szCompGE[1];
	  if( iSeries == iMaxSeries - 1 && iMaxSeries != Series.size() )
	    pszComp--;
	  ssTitle << pszSeries << pszComp << fixed << setprecision( iSeriesPrecision ) << Series[iSeries];
	  if( iSeriesPrecision != 1 )
	    ssTitle << setprecision(1);
	  ssTitle << " (n=" << iSeriesCount << ", x" <<  yScale << ")";
	  gF2DT[iSeries]->SetTitle( ssTitle.str().c_str() );

	  // Make the F2D/T Plot with a higher twist
	  gF2DTHT[iSeries] = new TGraphErrors( g_Count, g_x, g_F2DTHT, nullptr, g_F2DTHTErr );
	  if( bSingleton )
	    {
	      pszOption = szDrawSingle;
	      gF2DTHT[iSeries]->SetMarkerStyle( 41 );
	      gF2DTHT[iSeries]->SetMarkerSize( 0.75 );
	      gF2DTHT[iSeries]->SetMarkerColor( iColor );
	      gF2DTHT[iSeries]->SetLineColor( iColor );
	    }
	  else
	    {
	      pszOption = szDrawMulti;
	      //gF2DTHT[iSeries]->SetLineWidth( 3 );
	      gF2DTHT[iSeries]->SetLineColor( iColor );
	      gF2DTHT[iSeries]->SetFillColor( ( 1 - ( iSeries & 1 ) ) * 2 + 2 );
	      gF2DTHT[iSeries]->SetFillStyle( 3004 + ( iSeries & 1 ) );
	    }
	  ssTitle << " HT";
	  gF2DTHT[iSeries]->SetTitle( ssTitle.str().c_str() );

	  // Now add the plots to the multi-graph ... in the order we want them to print
	  mgF2DT->Add( gF2DT[iSeries], pszOption );
	  mgF2DT->Add( gF2DTHT[iSeries], pszOption );
	}
      //LOG( Rarely, "Printing graphs" << endl );
      TCanvas *c1 = new TCanvas("c1","Unused title - was ss1",200,10,700,500);
      c1->SetGrid();
      c1->GetFrame()->SetBorderSize(12);
      //c1->SetLogx();
      mgF2DT->Draw( "A3 L" );
      //mgF2DT->Draw("APZ");
      TLegend * pLeg = c1->BuildLegend( 0.75, 0.1, 1, 0.35 );
      if( pLeg )
	pLeg->SetNColumns( 2 );
      auto OldErrLvl = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kWarning;
      c1->Print(sFileName.c_str());
      gErrorIgnoreLevel = OldErrLvl;
      delete mgF2DT;
      if( pLeg )
	delete pLeg;
      delete c1;
      //for( int iSeries = 0 ; iSeries < iMaxSeries ; iSeries++ )
        //delete gF2DT[iSeries];
    }
}*/

void Cassandra::DataManager::PlotDiagonals( const TMatrixDSym & MData, const TMatrixDSym & MTheory,
					    const std::vector<DataNode *> & vNodes, double dNeff,
					    const std::string &sNamePrefix, const std::string &sFileNameSeq,
					    const std::string &sCut, const std::string &sCounts ) const
{
  Int_t MSize = ValidateMatrices( MData, MTheory, vNodes );
  if( MSize == 0 )
    return;

  // Write the data to a file
  std::string sFileName( sFileNameSeq );
  sFileName.append( "_TheoryExpt" );
  std::ofstream os(sFileName, ios_base::out | ios_base::trunc);
  if( os.fail() )
    { LOG( Always, "Unable to create " << sFileName << endl ); }
  else
    {
      // Write file header
      const char sp = ' ';
      os << "Seq DataSet target x Q**2 W**2 Var_Data Var_Theory Var_Data-Var_Theory" << endl;
      for( int i = 0 ; i < MSize ; i++ )
	os << i << sp << vNodes[i]->Type() << sp << vNodes[i]->Name()
	   << sp << fixed << setprecision(5) << vNodes[i]->m_x << sp << setprecision(3) << vNodes[i]->m_Q2
	   << sp << vNodes[i]->m_cW2 << scientific << setprecision(std::numeric_limits<double>::digits10)
	   << sp << MData[i][i] << sp << MTheory[i][i] << sp << MData[i][i] - MTheory[i][i]
	   << endl;
      os.close();
    }

  // Now plot as graph
  sFileName.append( ".pdf" );

  std::string sTitle( "Theory vs Experiment " );
  sTitle.append( sNamePrefix );
  sTitle.append( "; " );
  sTitle.append( std::to_string( MSize ) );
  sTitle.append( " data points (" );
  sTitle.append( sCounts );
  sTitle.append( "); Error x10,000 (" );
  sTitle.append( sCut );
  sTitle.append( ", " );
  sTitle.append( std::to_string( ( int ) ( dNeff + 0.5 ) ) );
  sTitle.append( " eff reps)" );

  // Prepare my charts
  const int iNumCharts = 2;
  Matrix2D ChartData( iNumCharts + 1, MSize );
  TGraph * g[iNumCharts];
  auto mg = new TMultiGraph();
  // Build a list of outliers across all charts
  std::vector<double> dMax( m_NumOutliers + 1, 0 );
  for( int iPass = 0 ; iPass < 2 ; iPass++ )
    {
      for( int c = 0; c < iNumCharts ; c++ ) // c for chart
	{
	  const TMatrixDSym * M = nullptr;
	  switch( c )
	    {
	    case 0:
	      M = & MData;
	      break;
	    case 1:
	      M = & MTheory;
	      break;
	    }
	  for( int i = 0 ; i < MSize ; i++ )
	    {
	      double d = (*M)[i][i] * 10000; // Current value
	      switch( iPass )
		{
		case 0:
		  // Last row has (common) x-values
		  if( c == 0 )
		    ChartData[iNumCharts][i] = i;
		  // First pass, build list of requested number of outliers
		  if( dMax[0] < d )
		    {
		      dMax[0] = d;
		      if( m_NumOutliers > 0 )
			std::sort( dMax.begin(), dMax.end() );
		      IF_LOG( Rarely )
		      {
			LOG( Always, "Highest + Outliers:" );
			for( int j = 0 ; j <= m_NumOutliers ; j++ )
			  LOG( Always, ' ' << dMax[j] );
			LOG( Always, endl );
		      }
		    }
		  break;
		case 1:
		  // Second pass, fill chart data with outliers discarded (if enough outliers)
		  if( dMax[0] > 0. && d > dMax[0] )
		    d = dMax[0] * 1.02;
		  ChartData[c][i] = d;
		  break;
		}
	    }
	}
    }

  // Log outliers ... in case it's necessary to look them up
  IF_LOG( Mostly )
  {
    LOG( Always, "Highest + Outliers:" );
    for( int j = 0 ; j <= m_NumOutliers ; j++ )
      LOG( Always, ' ' << dMax[j] );
    LOG( Always, endl );
  }

  // Now create each graph
  for( int c = 0; c < iNumCharts ; c++ ) // c for chart
    {
      static const char * ChartNames[] = { "Data", "Theory" };
      g[c] = new TGraph( MSize, ChartData[iNumCharts], ChartData[c] );
      switch( c )
	{
	case 0:
	  g[c]->SetMarkerStyle( kPlus );
	  g[c]->SetMarkerColor( kBlack );
	  break;
	case 1:
	  g[c]->SetMarkerStyle( kCircle );
	  g[c]->SetMarkerColor( kRed );
	  break;
	}
      g[c]->SetLineColor( kWhite );
      g[c]->SetBit( TH1::kNoStats );
      g[c]->SetTitle( ChartNames[c] );
      mg->Add( g[c] );
    }
  mg->SetTitle( sTitle.c_str() );
  auto c1 = new TCanvas("c1","Unused title - was ss1",200,10,700,500);
  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);
  mg->Draw("APZ");
  TLegend * pLeg = c1->BuildLegend( 0.45, 0.82, 0.55, 0.9 );
  // Now draw vertical lines between datasets
  bool bFirst = true;
  int nDat = 0;
  for( DataSet * pDS : m_DataSets )
    {
      if( pDS->m_list.size() > 0 )
	{
	  if( bFirst )
	    bFirst = false;
	  else
	    {
	      double dPos = nDat;
	      dPos /= MSize;
	      dPos = dPos * 0.8 + 0.1;
	      TLine * l = new TLine( dPos, 0.1, dPos, 0.9 );
	      l->SetLineColor( 7 );
	      l->SetLineWidth( 1 );
	      l->Draw();
	      delete l;
	    }
	  nDat += pDS->m_list.size();
	}
    }

  // Now print
  auto OldErrLvl = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  c1->Print( sFileName.c_str() );
  gErrorIgnoreLevel = OldErrLvl;
  //delete mg;
  for( int c = 0; c < iNumCharts ; c++ )
    delete g[c];
  delete pLeg;
  delete c1;
}

void Cassandra::DataManager::PlotCorrelation( const TMatrixDSym & M, double dNeff,
					      const std::string & sChartTitle,
					      const std::string & sFileName, const std::string &sCut,
					      const std::string & sCounts ) const
{
  Int_t MSize = M.GetNrows();
  if( MSize < 1 )
    {
      LOG( Always, "Error: Matrix only has " << MSize << " rows" << endl );
      return;
    }

  if( M.GetNcols() != MSize )
    {
      LOG( Always, "Error: " << MSize << " x " << M.GetNcols() << " matrix not symmetric" << endl );
      return;
    }

  std::string sTitle( sChartTitle );
  sTitle.append( "; " );
  sTitle.append( std::to_string( MSize ) );
  sTitle.append( " data points (" );
  sTitle.append( sCut );
  sTitle.append( ", " );
  sTitle.append( std::to_string( ( int ) ( dNeff + 0.5 ) ) );
  sTitle.append( " effective replicas); " );
  sTitle.append( std::to_string( MSize ) );
  sTitle.append( " data points (" );
  sTitle.append( sCounts.c_str() );
  sTitle.append( ")" );

  auto hist = new TH2D( "TheoryCorrel", sTitle.c_str(), MSize, -0.5, MSize - 0.5, MSize, -(MSize - 0.5), 0.5);
  // Now update my histograms
  int iPos = 0;
  int iNeg = 0;
  int iZero = 0;
  for( int y = 0 ; y < MSize ; y++ )
    for( int x = 0 ; x < MSize ; x++ )
      {
	double d = M[y][x];
	if( d == 0. )
	  {
	    d = 0.00000000001; // So it shows up on the plot
	    iZero++;
	  }
	else if( d < 0. )
	  iNeg++;
	else
	  iPos++;
	hist->Fill( x, -y, d );
      }
  hist->SetBit( TH1::kNoStats );
  LOG( Mostly, sFileName << " has " << iPos << " positive, "
       << iZero << " zero and " << iNeg << " negative values" << endl );
  auto c1 = new TCanvas("c1","Unused title",200,10,700,500);
  //c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);
  hist->Draw("COLZ");
  c1->Update();
  auto OldErrLvl = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  c1->Print( sFileName.c_str() );
  gErrorIgnoreLevel = OldErrLvl;
  delete hist;
  delete c1;
}

// Open the numbered model, or any earlier model

bool Cassandra::DataManager::LoadHTModel( std::vector<HT> &ht, int iModelNum ) const
{
  bool bOK = false;
  bool bFirst = true;
  std::ifstream is;
  ht.clear();
  do
    {
      std::string sModelFileName( m_pszModelPrefix );
      if( bFirst )
	{
	  // See whether model number has been specified
	  bFirst = false;
	  iModelNum++;
	}
      else
	{
	  sModelFileName.append( 1, '_' );
	  sModelFileName.append( std::to_string( iModelNum ) );
	}
      sModelFileName.append( m_szModelSuffix );
      is.open( sModelFileName );
      if( is.is_open() )
	{
	  bOK = true;
	  LOG( Mostly, "Loaded model " << sModelFileName << endl );
	}
    }
  while( !bOK && iModelNum-- );
  if( !bOK )
    { LOG( Always, "Unable to load any version of model " << m_pszModelPrefix << endl ); }
  else
    {
      std::string sLine;
      if( !std::getline( is, sLine ) || sLine.compare( m_szModelHeader ) )
	{
	  bOK = false;
	  LOG( Always, "Model " << iModelNum << " header invalid" << endl );
	}
      else
	{
	  int iSeq = 0;
	  while( bOK && std::getline( is, sLine ) )
	    {
	      std::stringstream ss( sLine );
	      int iFileSeq;
	      HT htFile;
	      htFile.chisq = -9.9e99;
	      ss >> iFileSeq >> htFile.p[0] >> htFile.p[1] >> htFile.p[2] >> htFile.p[3]
		 >> htFile.weight >> htFile.chisq;
	      if( iFileSeq != iSeq || htFile.chisq == -9.9e99)
		{
		  LOG( Always, "Error reading model record " << iSeq << endl );
		  bOK = false;
		}
	      else
		{
		  ht.push_back( htFile );
		  iSeq++;
		}
	    }
	  if( bOK && iSeq == 0 )
	    {
	      LOG( Always, "Model " << iModelNum << " is empty" << endl );
	      bOK = false;
	    }
	  if( !bOK )
	    ht.clear();
	}
      is.close();
    }
  return bOK;
}

// Make a prediction

bool Cassandra::DataManager::MakeOnePrediction( DataSet & ds, const char * pszOutFilePrefix, CutChi & cc,
						size_t SampleSize, bool /*bShowCut*/ )
{
  const char sp = ' ';
  static const char szSep[] = ": ";

  const char * pszCut = Cassandra_StrictString( cc.bStrict );
  std::stringstream ssCut;
  ssCut << "Q**2" << pszCut << fixed << setprecision(1) << cc.Q2 << ", W**2" << pszCut << cc.W2;

  std::string sFileNameSeq( pszOutFilePrefix );
  sFileNameSeq.append( 1, '_' );
  sFileNameSeq.append( std::to_string( cc.iSeq ) );

  m_ht.m_Model.clear();

  // Useful to print this again just before we dump the final lists
  LOG( Sometimes, cc.iSeq << szSep << ds.Name() << " data set contains "
       << ds.size() << " data points." << endl);

  // Dump the lists
  ds.DumpList(SampleSize);

  // Clear the combined covariance matrix
  for( int i = 0 ; i < ( int ) ds.CutSize() ; i++ )
    for( int j = 0 ; j < ( int ) ds.CutSize() ; j++ )
      {
	ds.m_MCovar[i][j] = 0.;
	ds.m_MCovarInv[i][j] = 0.;
      }
  // Construct the covariance matrix from the block diagonals of the component matrices
  int iNumData = 0;
  int iNumDataSets = 0;
  std::string sCounts;
  bool bFirst = true;
  for( DataSet * pDS : m_DataSets )
    {
      if( pDS->CutSize() )
	{
	  // Construct a string describing components of DataSet
	  if( bFirst )
	    bFirst = false;
	  else
	    sCounts.append( ", " );
	  sCounts.append( pDS->Name() );
	  sCounts.append( "=" );
	  sCounts.append( std::to_string( pDS->CutSize() ) );
	  // Construct the covariance matrix from the block diagonals of the component matrices
	  pDS->ConstructCovar();
	  for( int i = 0 ; i < ( int ) pDS->CutSize() ; i++ )
	    for( int j = 0 ; j < ( int ) pDS->CutSize() ; j++ )
	      {
		ds.m_MCovar   [iNumData + i][iNumData + j] = pDS->m_MCovar[i][j];
		ds.m_MCovarInv[iNumData + i][iNumData + j] = pDS->m_MCovarInv[i][j];
	      }
	  iNumData += pDS->CutSize();
	  iNumDataSets++;
	}
    }
  if( iNumData != ( int ) ds.CutSize() )
    {
      LOG( Always, cc.iSeq << ": Error constructing combined " << ds.CutSize() << "x" << ds.CutSize()
	   << " covariance matrix - component data sets only have " << iNumData
	   << " elements" << endl );
      return false;
    }
  LOG( Often, cc.iSeq << ": Combined covariance matrix is " << iNumData << "x" << iNumData << endl );

  // Plot the data sets
  {
    std::string sFileName( sFileNameSeq );
    sFileName.append( "_DataSet.pdf" );
    LOG( Mostly, "Creating " << sFileName << endl );
    auto mg = new TMultiGraph();
    int c = 0;
    for( DataSet * pDS : m_DataSets )
      if( pDS->CutSize() )
	{
	  static const int iMarkerStyles[] = { 22, 23, 29, 33, 39, 41, 47, 48 }; 
	  double x[pDS->CutSize()];
	  double y[pDS->CutSize()];
	  for( int i = 0 ; i < ( int ) pDS->CutSize() ; i++ )
	    {
	      x[i] = pDS->m_v[i]->m_x + ( c - iNumDataSets / 2. ) * 0.002;
	      y[i] = pDS->m_v[i]->m_Q2;
	    }
	  auto gr = new TGraph( ( Int_t ) pDS->CutSize(), x, y );
	  gr->SetTitle( pDS->Name() );
	  gr->SetMarkerStyle( iMarkerStyles[c % 8] );
	  int iColor;
	  if( c <= 2 )
	    iColor = c * 2 + 2;
	  else if( c <= 6 )
	    iColor = c * 2 - 3;
	  else
	    iColor = c + 13;
	  gr->SetMarkerColor( iColor );
	  gr->SetMarkerSize( 0.5 );
	  mg->Add( gr );
	  c++;
	}
    std::stringstream ssTitle;
    ssTitle << m_pszTarget << " data surviving cuts at " << ssCut.str().c_str()
	    << "; X (" << ds.CutSize() << " data points, " << sCounts.c_str() << "); Q**2 / GeV**2";
    mg->SetTitle( ssTitle.str().c_str() );
    TCanvas *c1 = new TCanvas("c1","Unused title - was ss1",200,10,700,500);
    c1->SetGrid();
    c1->GetFrame()->SetBorderSize(12);
    c1->SetLogy();
    mg->Draw( "AP" );
    TLegend * pLeg = c1->BuildLegend( 0.1, 0.8, 0.25, 0.94 );
    auto OldErrLvl = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;
    c1->Print( sFileName.c_str() );
    gErrorIgnoreLevel = OldErrLvl;
    if( pLeg )
      delete pLeg;
    delete mg;
    delete c1;
  }

  // Calculate chi squared goodness of fit before applying HT model
  double chisq = 0.0;
  for( int j = 0 ; j < iNumData ; j++)
    {
      double d = 0.0;
      for( int i = 0 ; i < iNumData ; i++)
	d += ( ds.m_v[i]->ScaledValue() - ds.m_v[i]->m_Theory[0] ) * ds.m_MCovarInv[i][j];
      chisq += d * ( ds.m_v[j]->ScaledValue() - ds.m_v[j]->m_Theory[0] );
    }
  LOG( Always, cc.iSeq << ": Without higher twist, reduced chi squared = " << fixed << setprecision(5)
       << chisq/iNumData << ", with " << iNumData << " data points" << endl );

  // Either load, or create the higher twist model
  int iNumHTReplicas = 0;
  double dNeff = 0.;
  bool bOK;
  if( m_pszModelPrefix )
    {
      // Load the model from file
      bOK = LoadHTModel( m_ht.m_Model, cc.iSeq );
      if( bOK )
	{
	  iNumHTReplicas = ( int ) m_ht.size();
	  // Calculate the weights across the HT replicas
	  double lnReplicas = std::log( iNumHTReplicas );
	  for( int k = 0 ; k < iNumHTReplicas ; k++ )
	    if( m_ht[k].weight != 0. )
	      dNeff += m_ht[k].weight * ( lnReplicas - std::log( m_ht[k].weight ) );
	  dNeff = std::exp( dNeff / iNumHTReplicas );
	}
    }
  else
    {
      // Create the model based on loaded data
      std::string sModelName( sFileNameSeq );
      sModelName.append( m_szModelSuffix );
      // Write out the model to chi table
      std::ofstream os(sModelName, ios_base::out | ios_base::trunc);
      if( os.fail() )
	{
	  bOK = false;
	  LOG( Always, cc.iSeq << ": Unable to create " << sModelName << endl );
	}
      else
	{
	  bOK = true;
	  // Write file header
	  const char sp = ' ';
	  os << fixed << setprecision(2)
	     << m_szModelHeader
	     << scientific << setprecision(std::numeric_limits<double>::digits10) << endl;

	  // Now calculate statistics for each variant of parameters
	  //std::vector<double> vTheoryAdjust( iNumData );
	  std::vector<double> vNewTheory( iNumData );
	  iNumHTReplicas = m_Params[0].Num * m_Params[1].Num * m_Params[2].Num * m_Params[3].Num;
	  m_ht.m_Model.resize( iNumHTReplicas );
	  double dMinWeight = 0;
	  int k = 0;
	  int ai, bi, ci, di;
	  double a, b, c, d;
	  for(ai=0, a=m_Params[0].Start; ai < m_Params[0].Num; a+=m_Params[0].Step, ai++)
	    for(bi=0, b=m_Params[1].Start; bi < m_Params[1].Num; b+=m_Params[1].Step, bi++)
	      for(ci=0, c=m_Params[2].Start; ci < m_Params[2].Num; c+=m_Params[2].Step, ci++)
		for(di=0, d=m_Params[3].Start; di < m_Params[3].Num; d+=m_Params[3].Step, di++, k++)
		  {
		    // Remember parameters
		    m_ht[k].p[0] = a;
		    m_ht[k].p[1] = b;
		    m_ht[k].p[2] = c;
		    m_ht[k].p[3] = d;
		    // Calculate adjustments to (multipliers of) the new theoretical predictions
		    for( int i = 0 ; i < iNumData ; i++)
		      {
			//double x = ds.m_v[i]->m_x;
			//double Q2 = ds.m_v[i]->m_Q2;
			//double Q4 = Q2 * Q2;
			//vTheoryAdjust[i] = 1. + m_ht.Adjust( ds.m_v[i], k );
			vNewTheory[i] = m_ht.NewTheory( ds.m_v[i], k );
		      }
		    // Calculate chi squared goodness of fit
		    m_ht[k].chisq = 0.;
		    for( int j = 0 ; j < iNumData ; j++ )
		      {
			double dTemp = 0.0;
			for( int i = 0 ; i < iNumData ; i++ )
			  dTemp += ( ds.m_v[i]->ScaledValue() - vNewTheory[i] ) * ds.m_MCovarInv[i][j];
			m_ht[k].chisq += dTemp * ( ds.m_v[j]->ScaledValue() - vNewTheory[j] );
		      }
		    // Calculate the weights across the HT replicas
		    double chi = sqrt( m_ht[k].chisq );
		    m_ht[k].weight = ( iNumData - 1 ) * std::log( chi ) - 0.5 * m_ht[k].chisq;
		    if( k == 0 || dMinWeight > m_ht[k].weight )
		      dMinWeight = m_ht[k].weight;
		  }

	  // Divide all the weights by the minimum weight
	  double dTotalWeight = 0.;
	  for( int k = 0 ; bOK && k < iNumHTReplicas ; k++ )
	    {
	      std::feclearexcept( FE_ALL_EXCEPT );
	      m_ht[k].weight = std::exp( m_ht[k].weight - dMinWeight - 735 );
	      if( std::fetestexcept( FE_OVERFLOW ) )
		{
		  bOK = false;
		  LOG( Always, cc.iSeq << ": HT Replica " << k
		       << " overflow error while calculating chi squared" << endl );
		}
	      else
		dTotalWeight += m_ht[k].weight;
	    }
	  if( bOK )
	    {
	      // Calculate the weights across the HT replicas
	      dTotalWeight /= iNumHTReplicas;
	      double lnReplicas = std::log( iNumHTReplicas );
	      double dNewTotalWeight = 0;
	      for( int k = 0 ; k < iNumHTReplicas ; k++ )
		{
		  if( dTotalWeight != 0. )
		    m_ht[k].weight /= dTotalWeight;
		  //if( m_ht[k].weight < 0.1 )
		  //m_ht[k].weight = 0.;  // Discard low-weight replicas
		  //else
		  if( m_ht[k].weight != 0. )
		    {
		      dNewTotalWeight += m_ht[k].weight;
		      dNeff += m_ht[k].weight * ( lnReplicas - std::log( m_ht[k].weight ) );
		    }
		  // Now write out the chi squared and the weights
		  os << k << sp << m_ht[k].p[0] << sp << m_ht[k].p[1] << sp << m_ht[k].p[2] << sp
		     << m_ht[k].p[3] << sp << m_ht[k].weight << sp << m_ht[k].chisq / iNumData << endl;
		}
	      dNeff = std::exp( dNeff / iNumHTReplicas );

	      // If we're doing a detailed log, show the higher twists
	      IF_LOG( Rarely )
	      {
		LOG( Always, "Building model" << endl );
		double HTMean[iNumData];
		for( int i = 0 ; i < iNumData ; i++ )
		  {
		    HTMean[i] = 0.; // Mean for this datapoint
		    for( int k = 0 ; k < iNumHTReplicas ; k++ )
		      HTMean[i] += m_ht[k].weight * m_ht.HigherTwist( ds.m_v[i], k );
		    HTMean[i] /= dNewTotalWeight;
		    LOG( Always, cc.iSeq << szSep << i << sp << ds.m_v[i]->Type()
			 << " x=" << fixed << setprecision(3) << ds.m_v[i]->m_x
			 << ", Q**2=" << setprecision(1) << ds.m_v[i]->m_Q2 << setprecision(6)
			 << ", F2D=" << ds.m_v[i]->ScaledValue() << ", F2T[central]=" << ds.m_v[i]->m_Theory[0]
			 << setprecision(std::numeric_limits<double>::digits10)
			 << ", HT " << HTMean[i] << endl );
		  }
	      }

	      // I might as well graph the parameter weights
	      TH1D * tp[4];
	      for( int i = 0; i < 4; i++ )
		{
		  if( m_Params[i].Num > 1 )
		    {
		      static char szStats[] = "Stats ?";
		      static char szParam[] = "Parameter ?";
		      char c = 'A' + i;
		      szStats[6] = c;
		      szParam[10] = c;
		      tp[i] = new TH1D( szStats, szParam, m_Params[i].Num,
					m_Params[i].Start - m_Params[i].Step * 0.5,
					m_Params[i].Start + m_Params[i].Step * ( m_Params[i].Num - 0.5 ) );
		      // Now update my histograms
		      for( int k = 0 ; k < iNumHTReplicas ; k++ )
			tp[i]->Fill( m_ht[k].p[i], m_ht[k].weight );
		      auto c1 = new TCanvas("c1","Unused title - was ss1",200,10,700,500);
		      c1->SetGrid();
		      c1->GetFrame()->SetBorderSize(12);
		      std::string sName(sModelName);
		      static char szSuffix[] = "_?.pdf";
		      szSuffix[1] = 'A' + i;
		      sName.append( szSuffix );
		      std::stringstream ss;
		      ss << fixed << setprecision(0)
			 << tp[i]->GetTitle() << " (" << ssCut.str() << ", "
			 << dNeff << " effective replicas)";
		      tp[i]->SetTitle( ss.str().c_str() );
		      tp[i]->Draw("hist");
		      c1->Update();
		      TPaveStats * st = (TPaveStats*) tp[i]->FindObject("stats");
		      if( st )
			{
			  //double a = st->GetX2NDC() - st->GetX1NDC();
			  //st->SetX1NDC( 0.7 * a );
			  //st->SetX2NDC( 1.7 * a );
			  double a = st->GetY2NDC() - st->GetY1NDC();
			  st->SetY1NDC( 0.8 * a );
			  st->SetY2NDC( 1.8 * a );
			}
		      auto OldErrLvl = gErrorIgnoreLevel;
		      gErrorIgnoreLevel = kWarning;
		      c1->Print( sName.c_str() );
		      gErrorIgnoreLevel = OldErrLvl;
		      if( st != nullptr )
			delete st;
		      delete c1;
		      delete tp[i];
		    }
		}
	    }
	  os.close();
	}
    }

  // Error loading / creating model?
  if( !bOK )
    return false;

  // Calculate total weight for normalisation. Make sure it's not zero
  double dTotalWeight = 0.;
  for( int k = 0 ; k < iNumHTReplicas ; k++ )
    dTotalWeight += m_ht[k].weight;
  if( dTotalWeight <= 0. )
    {
      LOG( Always, cc.iSeq << ": Error: model invalid, total replica weight " << dTotalWeight << endl );
      return false;
    }

  // Create matrix of the higher twist for each data point by every model
  Matrix2D MDataTwist( iNumData, iNumHTReplicas );
  double HTMean[ iNumData ]; // Higher twist mean for each data point
  LOG( Rarely, "Recomputed model" << endl );
  for( int i = 0 ; i < iNumData ; i++ )
    {
      HTMean[i] = 0.; // Mean for this datapoint
      for( int k = 0 ; k < iNumHTReplicas ; k++ )
	{
	  MDataTwist[i][k] = m_ht.HigherTwist( ds.m_v[i], k );
	  HTMean[i] += m_ht[k].weight * MDataTwist[i][k];
	  //cout << sp << MDataTwist[i][k];
	}
      HTMean[i] /= dTotalWeight;
      LOG( Rarely, cc.iSeq << szSep << i << sp << ds.m_v[i]->Type()
	   << " x=" << fixed << setprecision(3) << ds.m_v[i]->m_x
	   << ", Q**2=" << setprecision(1) << ds.m_v[i]->m_Q2 << setprecision(6)
	   << ", F2D=" << ds.m_v[i]->ScaledValue() << ", F2T[central]=" << ds.m_v[i]->m_Theory[0]
	   << setprecision(std::numeric_limits<double>::digits10)
	   << ", HT " << HTMean[i] << endl );
    }

  // Create theory-covariance matrix
  TMatrixDSym MTheoryCovar( iNumData );
  for( int i = 0 ; i < iNumData ; i++ )
    for( int j = 0 ; j <= i ; j++ )
      {
	double d = 0.;
	for( int k = 0 ; k < iNumHTReplicas ; k++ )
	  d += m_ht[k].weight * ( MDataTwist[i][k] - HTMean[i] ) * ( MDataTwist[j][k] - HTMean[j] );
	d /= dTotalWeight;
	MTheoryCovar[i][j] = d;
	if( i != j )
	  MTheoryCovar[j][i] = d;
      }

  // Create and plot theory-correlation matrix
  TMatrixDSym MTheoryCorrel( iNumData );
  for( int i = 0 ; i < iNumData ; i++ )
    for( int j = 0 ; j <= i ; j++ )
      {
	double d = MTheoryCovar[i][j] / sqrt( MTheoryCovar[i][i] * MTheoryCovar[j][j] );
	MTheoryCorrel[i][j] = d;
	if( i != j )
	  MTheoryCorrel[j][i] = d;
      }
  std::string sChartTitle( "Theory Correlation " );
  sChartTitle.append( ds.Name() );
  std::string sChartName( sFileNameSeq );
  sChartName.append( "_CorrelTheory.pdf" );
  PlotCorrelation( MTheoryCorrel, dNeff, sChartTitle, sChartName, ssCut.str(), sCounts );

  // Create and plot experiment-correlation matrix
  TMatrixDSym MExptCorrel( iNumData );
  for( int i = 0 ; i < iNumData ; i++ )
    for( int j = 0 ; j <= i ; j++ )
      {
	double d = ds.m_MCovar[i][j] / sqrt( ds.m_MCovar[i][i] * ds.m_MCovar[j][j] );
	MExptCorrel[i][j] = d;
	if( i != j )
	  MExptCorrel[j][i] = d;
      }
  sChartTitle = "Experiment Correlation ";
  sChartTitle.append( ds.Name() );
  sChartName = sFileNameSeq;
  sChartName.append( "_CorrelExpt.pdf" );
  PlotCorrelation( MExptCorrel, dNeff, sChartTitle, sChartName, ssCut.str(), sCounts );

  // Create and plot theory-experiment-correlation matrix
  TMatrixDSym MExptTheoryCorrel( iNumData );
  for( int i = 0 ; i < iNumData ; i++ )
    for( int j = 0 ; j <= i ; j++ )
      {
	double d = ( MTheoryCovar[i][j] + ds.m_MCovar[i][j] )
	  / sqrt(( MTheoryCovar[i][i]+ds.m_MCovar[i][i] )*( MTheoryCovar[j][j]+ds.m_MCovar[j][j] ));
	MExptTheoryCorrel[i][j] = d;
	if( i != j )
	  MExptTheoryCorrel[j][i] = d;
      }
  sChartTitle = "Experiment+Theory Correlation ";
  sChartTitle.append( ds.Name() );
  sChartName = sFileNameSeq;
  sChartName.append( "_CorrelExptTheory.pdf" );
  PlotCorrelation( MExptTheoryCorrel, dNeff, sChartTitle, sChartName, ssCut.str(), sCounts );

  // Now calculate chi squared goodness of fit for HT model
  double chisqHT = 0.;
  double chisqHTReweight = 0.; // We'll calculate this later if we reweight PDFs
  for( int j = 0 ; j < iNumData ; j++)
    {
      double d = 0.0;
      for( int i = 0 ; i < iNumData ; i++)
	d += ( ds.m_v[i]->ScaledValue() - ( ds.m_v[i]->m_Theory[0] + HTMean[i] ) )
	  * ds.m_MCovarInv[i][j];
      chisqHT += d * ( ds.m_v[j]->ScaledValue() - ( ds.m_v[j]->m_Theory[0] + HTMean[j] ) );
    }
  chisqHT /= iNumData;

  // I will need to append this prediction to chi table ... but will do at end
  string sChiFile( pszOutFilePrefix );
  sChiFile += m_pszChiSuffix;
  std::ofstream os( sChiFile, ios_base::out | ios_base::app );
  if( os.fail() )
    {
      LOG( Always, cc.iSeq << ": Unable to append to " << sChiFile << endl );
      return false;
    }

  // If we're performing a detailed log, summarize higher twists
  IF_LOG( Rarely )
  {
    for( int i = 0 ; i < iNumData ; i++ )
      {
	int iPos = 0;
	int iNeg = 0;
	int iZero = 0;
	for( int k = 0 ; k < iNumHTReplicas ; k++ )
	  {
	    double d = MDataTwist[i][k] - HTMean[i];
	    if( d < 0 )
	      iNeg++;
	    else if( d == 0 )
	      iZero++;
	    else
	      iPos++;
	  }
	LOG( Always, cc.iSeq << szSep << i << szSep << iPos << " positive, "
	     << iNeg << " negative, " << iZero << "=0" << endl );
      }
  }

  // Plot the higher twist
  PlotHT( ds, sFileNameSeq, ssCut.str(), HTMean, MTheoryCovar );

  // Create graphs of the data - this time using the model
  CreateGraphs( ds.m_MCovar, MTheoryCovar, ds.m_v, dNeff,
		ds.Name(), sFileNameSeq, ssCut.str(), sCounts, HTMean );

  // Plot the diagonal experimental and theoretical errors
  PlotDiagonals( ds.m_MCovar, MTheoryCovar, ds.m_v, dNeff,
		 ds.Name(), sFileNameSeq, ssCut.str(), sCounts );

  // If I'm running a pre-existing model with more data, determine effect on PDFs
  if( m_pszModelPrefix && ds.m_iNumReplicas > 1 )
    {
      // Calculate PDF reweighting
      double vWeight[ds.m_iNumReplicas];
      std::string sPDFWeightName( pszOutFilePrefix );
      sPDFWeightName.append( m_szPDFSuffix );
      std::ofstream wos( sPDFWeightName, ios_base::out | ios_base::app );
      if( wos.fail() )
	{
	  LOG( Always, "Unable to create " << sPDFWeightName << endl );
	  bOK = false;
	}
      else
	{
	  LOG( Often, "Creating " << sPDFWeightName << endl );

	  // Start by inverting the theory + experiment covariance matrix
	  LOG( Often, "Inverting " << iNumData << " by " << iNumData
	       << " theory + experiment covariance matrix" << endl );
	  TMatrixDSym MExptTheoryCovar( iNumData );
	  TMatrixDSym MExptTheoryCovarInv( iNumData );
	  for( int i = 0 ; i < iNumData ; i++ )
	    for( int j = 0 ; j < iNumData ; j++ )
	      {
		double d = MTheoryCovar[i][j] + ds.m_MCovar[i][j];
		MExptTheoryCovar[i][j] = d;
		MExptTheoryCovarInv[i][j] = d;
	      }
	  Double_t det=-9.9e99;
	  auto tStart = std::chrono::high_resolution_clock::now();
	  MExptTheoryCovarInv.Invert(&det);
	  auto tEnd = std::chrono::high_resolution_clock::now();
	  std::chrono::duration<double> tElapsed = tEnd - tStart;
	  LOG( Often, "Matrix inversion took " << fixed << tElapsed.count() << " seconds." << endl );

	  double vNewTheory[iNumData];
	  double vChi[ds.m_iNumReplicas];
	  double dMinWeight = 0;
	  // Calculate chi squared for each replica
	  for( int iRep = 1 ; iRep < ds.m_iNumReplicas ; iRep++ )
	    {
	      // Calculate the new theoretical predictions
	      for( int i = 0 ; i < iNumData ; i++)
		vNewTheory[i] = ds.m_v[i]->m_Theory[iRep] + HTMean[i];

	      // Calculate chi squared goodness of fit
	      vChi[iRep] = 0;
	      for( int j = 0 ; j < iNumData ; j++ )
		{
		  double d = 0.0;
		  for( int i = 0 ; i < iNumData ; i++ )
		    d += ( ds.m_v[i]->ScaledValue() - ( ds.m_v[i]->m_Theory[iRep] + HTMean[i] ) )
		      * MExptTheoryCovarInv[i][j];
		  vChi[iRep] += d * ( ds.m_v[j]->ScaledValue() - ( ds.m_v[j]->m_Theory[iRep] + HTMean[j] ) );
		}
	      // Now calculate the weight
	      double chi = sqrt( vChi[iRep] );
	      LOG( Rarely, scientific << setprecision(std::numeric_limits<double>::digits10)
		   << iRep << ": vChi[iRep]=" << vChi[iRep]
		   << ", chi=" << chi
		   << ", log(chi)=" << std::log(chi)
		   << ", (iNumData-1)=" << iNumData - 1
		   << endl );
	      //vWeight[iRep] = std::exp( ( iNumData - 1 ) * std::log( chi ) - 0.5 * vChi[iRep] );
	      vWeight[iRep] = ( iNumData - 1 ) * std::log( chi ) - 0.5 * vChi[iRep];
	      if( iRep == 1 || dMinWeight > vWeight[iRep] )
		dMinWeight = vWeight[iRep];
	    }
	  dNeff = 0.;
	  dTotalWeight = 0.;
	  // Divide all the weights by the minimum weight
	  for( int iRep = 1 ; iRep < ds.m_iNumReplicas ; iRep++ )
	    {
	      std::feclearexcept( FE_ALL_EXCEPT );
	      vWeight[iRep] = std::exp( vWeight[iRep] - dMinWeight -735 );
	      if( std::fetestexcept( FE_OVERFLOW ) )
		{
		  if( bOK )
		    {
		      bOK = false;
		      LOG( Always, cc.iSeq << ": Replica " << iRep
			   << " overflow error while calculating chi squared" << endl );
		    }
		  vWeight[iRep] = 0.;
		}
	      else
		dTotalWeight += vWeight[iRep];
	    }
	  if( bOK )
	    {
	      // Calculate the weights and number of effective replicas
	      dTotalWeight /= ds.m_iNumReplicas - 1.;
	      double dNewTotalWeight = 0;
	      // We'll repeatedly need the log of the number of replicas
	      // Ignoring replica 0, i.e. central value
	      double lnReplicas = std::log( ds.m_iNumReplicas - 1 );
	      for( int iRep = 1 ; iRep < ds.m_iNumReplicas ; iRep++ )
		{
		  vWeight[iRep] /= dTotalWeight;
		  if( vWeight[iRep] != 0. )
		    dNeff += vWeight[iRep] * ( lnReplicas - std::log( vWeight[iRep] ) );
		  dNewTotalWeight += vWeight[iRep];
		  chisqHTReweight += vChi[iRep] * vWeight[iRep];
		}
	      dNeff = std::exp( dNeff / ( ds.m_iNumReplicas - 1 ) );
	      dTotalWeight = dNewTotalWeight;
	      chisqHTReweight /= dTotalWeight * iNumData;
	    }

	  // Write out what we know whether successful or not - useful for debugging
	  // Now write out the chi squared and the weights
	  wos << cc.iSeq << Cassandra_StrictString( cc.bStrict )
	     << fixed << setprecision(5) << m_XCutoff << sp << cc.W2 << sp << cc.Q2
	     << sp << iNumData << sp << dNeff
	     << scientific << setprecision(std::numeric_limits<double>::digits10);
	  for( int iRep = 1 ; iRep < ds.m_iNumReplicas ; iRep++ )
	    wos << sp << vWeight[iRep];
	  for( int iRep = 1 ; iRep < ds.m_iNumReplicas ; iRep++ )
	    wos << sp << vChi[iRep] / iNumData;
	  wos << endl;
	  wos.close();

	  // If there's been an error, or a PDF set hasn't been specified, that's all I need
	  if( bOK && m_pszPDFSet != nullptr )
	    {
	      // Now let's plot the impact on PDFs
	      // References for Q0 = 1.65 GeV for NNPDF 3.1
	      // ArXiv 1706.00428, pg 20
	      // ArXiv 1605.06515, equation 2, pg 7
	      double eps = 1e-10;
	      double Q0 = 1.65;
	      double QMin = m_dQChartScale;
	      double QMax = m_dQChartScale;
	      if( QMin >= Q0 )
		QMin = Q0 - eps;
	      if( QMax <= Q0 )
		QMax = Q0 + eps;

	      // Build a list of x values we are interested in
	      std::vector<double> x;
	      for( int i = -25 ; i <= 0 ; i++ )
		{
		  if( i == 0 )
		    x.push_back( 0.9 );
		  else
		    x.push_back( std::pow( 10, i / 5. ) );
		}

	      // These are the partons we're interested in
	      Parton Partons[] = { {2, "up", "Up quark"},
				   {1, "down", "Down quark"},
				   {0, "gluon", "Gluon"},
				   {4, "charm", "Charm quark"},
				   {-2, "ubar", "Anti-Up quark"},
				   {-1, "dbar", "Anti-Down quark"},
				   {-4, "cbar", "Anti-Charm quark"}
	      };
	      const int iNumPartons = sizeof( Partons ) / sizeof( Parton );

	      // We'll need space for every parton's data
	      Matrix2D MPartonData[iNumPartons];
	      for( int i = 0 ; i < iNumPartons ; i++ )
		MPartonData[i].Resize( ( int ) x.size(), ds.m_iNumReplicas );

	      // Load evolution
	      IF_LOG(Mostly)
		;
	      else
		APFEL::EnableWelcomeMessage( false );
	      APFEL::SetSmallxResummation(0, "NLL"); // Do I need this? Added because sample had it
	      APFEL::SetPDFSet( m_pszPDFSet );
	      // NNPDF sets are for five quark flavours
	      int iNumQuarks = 5;
	      APFEL::SetMaxFlavourPDFs(iNumQuarks);  // # flavours in DGLAP evolution
	      APFEL::SetMaxFlavourAlpha(iNumQuarks); // # flavours in coupling equations
	      APFEL::SetMassScheme("FONLL-C"); // Per Prof Ball 2018-06-17
	      APFEL::SetQLimits( QMin, QMax );
	      APFEL::EnableTargetMassCorrections(true); // Per Prof Ball 2018-06-08
	      APFEL::EnableDampingFONLL(false); // Per Prof Ball 2018-06-08
	      APFEL::EnableLeptonEvolution(true);
	      APFEL::SetPDFEvolution("truncated"); // Per Prof Ball 2018-06-08
	      APFEL::EnableIntrinsicCharm(true); // Per Prof Ball 2018-06-08
	      APFEL::InitializeAPFEL();

	      // Perform predictions for requested PDF replicas
	      // Backwards so we fail early if there aren't enough replicas in the PDF
	      LOG( Often, "Evaluating impact to PDF at Q=" << m_dQChartScale << " GeV" << endl );
	      dTotalWeight = 0.;
	      for( int iRep = ds.m_iNumReplicas - 1 ; iRep > 0 ; iRep-- )
		{
		  APFEL::SetReplica( iRep );
		  APFEL::EvolveAPFEL( Q0, m_dQChartScale );
		  // Walk the list of x-values
		  for( int i = 0 ; i < ( int ) x.size() ; i++ )
		    {
		      double xf[13];
		      APFEL::xPDFall( x[i], xf );
		      // Save data for all the partons we're interested in
		      for( int j = 0 ; j < iNumPartons ; j++ )
			MPartonData[j][i][iRep] = xf[Partons[j].m_APFELIndex + 6];
		    }
		  dTotalWeight += vWeight[iRep];
		}

	      // Now create charts of impact to PDF
	      for( int c = 0 ; c < iNumPartons ; c++ )
		{
		  // First, compute statistics
		  double y_Orig[x.size()];
		  double y_OrigErr[x.size()];
		  double y_Mod[x.size()];
		  double y_ModErr[x.size()];
		  double y_ErrRatio[x.size()];
		  for( int i = 0 ; i < ( int ) x.size() ; i++ )
		    {
		      y_Orig[i] = 0.;
		      y_OrigErr[i] = 0.;
		      y_Mod[i] = 0.;
		      y_ModErr[i] = 0.;
		      // Calculate means first ... to avoid overflow
		      for( int iRep = 1 ; iRep < ds.m_iNumReplicas ; iRep++ )
			{
			  double d = MPartonData[c][i][iRep];
			  y_Orig[i] += d;
			  y_Mod[i] += vWeight[iRep] * d;
			}
		      y_Orig[i] /= ( ds.m_iNumReplicas - 1 );
		      y_Mod[i] /= dTotalWeight;
		      // Now calculate error
		      for( int iRep = 1 ; iRep < ds.m_iNumReplicas ; iRep++ )
			{
			  double d = MPartonData[c][i][iRep];
			  y_OrigErr[i] += ( d - y_Orig[i] ) * ( d - y_Orig[i] );
			  y_ModErr[i] += vWeight[iRep] * ( d - y_Mod[i] ) * ( d - y_Mod[i] );
			}
		      y_OrigErr[i] = std::sqrt( y_OrigErr[i] / ( ds.m_iNumReplicas - 1 ) );
		      y_ModErr[i] = std::sqrt( y_ModErr[i] / dTotalWeight );
		      // cout << cc.iSeq << szSep << Partons[c].m_szName << " raw  " << i
		      // << " x=" << x[i]
		      // << " Orig=" << y_Orig[i] << " +/- " << y_OrigErr[i]
		      // << " Mod=" << y_Mod[i] << " +/- " << y_ModErr[i]
		      // << endl;
		      // Calculate ratio of new error to old
		      if( y_OrigErr[i] == 0. )
			y_ErrRatio[i] = 1E300;
		      else
			y_ErrRatio[i] = y_ModErr[i] / y_OrigErr[i];
		      // Now normalise to original values, keeping relative errors
		      if( y_Orig[i] != 0. )
			{
			  y_Mod[i] /= y_Orig[i];
			  y_ModErr[i] /= y_Orig[i];
			  y_OrigErr[i] /= y_Orig[i];
			  y_Orig[i] = 1.;
			}
		      // Debug
		      LOG( Rarely, cc.iSeq << szSep << Partons[c].m_szName << " norm " << i
			   << " x=" << x[i]
			   << " Orig=" << y_Orig[i] << " +/- " << y_OrigErr[i]
			   << " Mod=" << y_Mod[i] << " +/- " << y_ModErr[i]
			   << endl );
		      //y_OrigErr[i] = 0.1 + 0.1 * ( ( i + 1 ) % 2 );
		      //y_ModErr[i] = 0.1 + 0.1 * ( i % 2 );
		      //y_Mod[i] = 1.1;
		    }
		  // Now create charts
		  auto g_Orig = new TGraphErrors( ( Int_t ) x.size(), &x[0], y_Orig, nullptr, y_OrigErr );
		  //g_Orig->SetLineStyle( 2 );
		  g_Orig->SetLineColor( kBlack );
		  g_Orig->SetFillColor( 2 );
		  g_Orig->SetFillStyle( 3002 );
		  g_Orig->SetTitle( "NNPDF 3.1" );
		  auto g_Mod = new TGraphErrors( ( Int_t ) x.size(), &x[0], y_Mod, nullptr, y_ModErr );
		  g_Mod->SetLineColor( 2 );
		  g_Mod->SetFillColor( 4 );
		  g_Mod->SetFillStyle( 3004 );
		  g_Mod->SetTitle( "+HT" );
		  auto g_ErrRatio = new TGraph( ( Int_t ) x.size(), &x[0], y_ErrRatio );
		  g_ErrRatio->SetLineColor( 7 ); // acqua
		  g_ErrRatio->SetMarkerColor( 7 ); // acqua
		  g_ErrRatio->SetTitle( "Rel. Error" );
		  auto mg = new TMultiGraph();
		  mg->Add( g_Orig, "" );
		  mg->Add( g_Mod, "" );
		  mg->Add( g_ErrRatio, "C" );

		  TCanvas *c1 = new TCanvas("c1","Unused title - was ss1",200,10,700,500);
		  c1->SetGrid();
		  c1->GetFrame()->SetBorderSize(12);
		  c1->SetLogx();
		  //g_Orig->Draw( "A3 L" );
		  //g_Mod->Draw( "A3 L SAME" );
		  auto OldErrLvl = gErrorIgnoreLevel;
		  gErrorIgnoreLevel = kWarning;
		  {
		    std::stringstream sTitle;
		    sTitle << Partons[c].m_szChartName << ' ' << ds.Name()
			   << " (Q=" << fixed << setprecision(1) << m_dQChartScale << " GeV); X ("
			   << ds.CutSize() << " data points, " << ssCut.str().c_str()
			   << "); x * PDF (NNPDF 3.1=1)";
		    mg->SetTitle( sTitle.str().c_str() );
		    mg->Draw( "A3 L" );
		    c1->Update();
		    mg->SetMinimum( 0.8 );
		    mg->SetMaximum( 1.3 );
		    //TLegend * pLeg = c1->BuildLegend( 0.4, 0.74, 0.6, 0.9 );

		    //sChartTitle = Partons[c].m_szChartName;
		    //sChartTitle.append( ds.Name() );
		    //g_Orig->SetTitle( sChartTitle.c_str() );
		    //g_Orig->GetXaxis()->SetTitle( "NNPDF 3.1=1" );
		    //auto t = g_Orig->GetYaxis();
		    //t->SetTitle( "NNPDF 3.1 = 1" );
		    //t->SetTitleOffset( 1.5 );
		    std::string sFileName( sFileNameSeq );
		    sFileName.append( "_q_" );
		    sFileName.append( Partons[c].m_szName );
		    sFileName.append( ".pdf" );
		    c1->Print(sFileName.c_str());
		  }
		  gErrorIgnoreLevel = OldErrLvl;
		  //delete g_Mod;
		  //delete g_Orig;
		  delete mg;
		  delete c1;
		}
	    }
	}
    }

  // Append this prediction to chi table
  os << cc.iSeq << Cassandra_StrictString( cc.bStrict )
     << fixed << setprecision(3) << m_XCutoff
     << setprecision(2) << sp << cc.W2 << sp << cc.Q2
     << setprecision(std::numeric_limits<double>::digits10)
     << sp << chisq << sp << iNumData << sp << chisq / iNumData
     << sp << chisqHT << sp << chisqHTReweight
     << setprecision(1) << sp << dNeff
     << endl;
  os.close();
  cc.Chi = chisq / iNumData;
  cc.ChiTwist = chisqHT;
  cc.ChiTwistReweight = chisqHTReweight;
  return bOK;
}

// Sample code from https://github.com/scarrazza/apfel/blob/master/examples/TabulationCxx.cc

void APFELTest( const char * pszPDFSet )
{
  double xlha[] = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 
		1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  
  // Activate some options
  //APFEL::SetFastEvolution(true);
  //APFEL::SetPerturbativeOrder(0);
  //APFEL::SetPDFSet("MRST2004qed");
  //APFEL::EnableEvolutionOperator(true);
  // Initializes integrals on the grids
  APFEL::SetTheory("QCD");
  APFEL::SetSmallxResummation(0, "NLL");
  APFEL::SetPDFSet( pszPDFSet );
  APFEL::InitializeAPFEL();

  double Q02, Q2, eps = 1e-10;
  cout << "Enter initial and final scale in GeV^2" << endl;
  cin >> Q02 >> Q2;

  // Load evolution
  double Q0 = sqrt(Q02) - eps;
  double Q  = sqrt(Q2);
  APFEL::EvolveAPFEL(Q0,Q);

  cout << scientific << setprecision(5) << endl;
  // Tabulate PDFs for the LHA x values
  cout << "alpha_QCD(mu2F) = " << APFEL::AlphaQCD(Q) << endl;
  cout << "alpha_QED(mu2F) = " << APFEL::AlphaQED(Q) << endl;
  cout << endl;

  cout << "Standard evolution:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " 
       << setw(11) << "     photon  "
       << setw(11) << "    e^-+e^+  "
       << setw(11) << "   mu^-+mu^+ "
       << setw(11) <<"   tau^-+tau^+" << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::xPDFj(2,xlha[i]) - APFEL::xPDFj(-2,xlha[i]) << "  "
	 << setw(11) << APFEL::xPDFj(1,xlha[i]) - APFEL::xPDFj(-1,xlha[i]) << "  "
	 << setw(11) << 2*(APFEL::xPDFj(-1,xlha[i]) + APFEL::xPDFj(-2,xlha[i])) << "  "
	 << setw(11) << APFEL::xPDFj(4,xlha[i]) + APFEL::xPDFj(-4,xlha[i]) << "  "
	 << setw(11) << APFEL::xPDFj(0,xlha[i]) << "  "
	 << setw(11) << APFEL::xgammaj(xlha[i]) << "  "
	 << setw(11) << APFEL::xLeptonj(1,xlha[i]) + APFEL::xLeptonj(-1,xlha[i]) << "  "
	 << setw(11) << APFEL::xLeptonj(2,xlha[i]) + APFEL::xLeptonj(-2,xlha[i]) << "  "
	 << setw(11) << APFEL::xLeptonj(3,xlha[i]) + APFEL::xLeptonj(-3,xlha[i]) << "  "
	 << endl;
  cout << "      " << endl;

  cout << "Standard evolution using the xPDFall function:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " << endl;

  cout << scientific;
  double xf[13];
  for (int i = 2; i < 11; i++) {
    APFEL::xPDFall(xlha[i],xf);
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << xf[2+6] - xf[-2+6] << "  "
	 << setw(11) << xf[1+6] - xf[-1+6] << "  "
	 << setw(11) << 2*(xf[-1+6] + xf[-2+6]) << "  "
	 << setw(11) << xf[4+6] + xf[-4+6] << "  "
	 << setw(11) << xf[0+6] << "  "
	 << endl;
  }
  cout << "      " << endl;
  //
  // Cached PDFs
  //
  APFEL::CachePDFsAPFEL(Q0);

  cout << "Cached evolution:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " 
       << setw(11) << "     photon  "
       << setw(11) << "    e^-+e^+  "
       << setw(11) << "   mu^-+mu^+ "
       << setw(11) <<"   tau^-+tau^+" << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++)
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << APFEL::xPDFxQ(2,xlha[i],Q) - APFEL::xPDFxQ(-2,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(1,xlha[i],Q) - APFEL::xPDFxQ(-1,xlha[i],Q) << "  "
	 << setw(11) << 2*(APFEL::xPDFxQ(-1,xlha[i],Q) + APFEL::xPDFxQ(-2,xlha[i],Q)) << "  "
	 << setw(11) << APFEL::xPDFxQ(4,xlha[i],Q) + APFEL::xPDFxQ(-4,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(0,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(22,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(11,xlha[i],Q) + APFEL::xPDFxQ(-11,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(13,xlha[i],Q) + APFEL::xPDFxQ(-13,xlha[i],Q) << "  "
	 << setw(11) << APFEL::xPDFxQ(15,xlha[i],Q) + APFEL::xPDFxQ(-15,xlha[i],Q) << "  "
	 << endl;
  cout << "      " << endl;

  cout << "Cached evolution using the xPDFxQall function:" << endl;
  cout << "   x   " 
       << setw(11) << "    u-ubar    " 
       << setw(11) << "   d-dbar    " 
       << setw(11) << " 2(ubr+dbr)  " 
       << setw(11) << " c+cbar  " 
       << setw(11) << "     gluon   " << endl;

  cout << scientific;
  for (int i = 2; i < 11; i++) {
    APFEL::xPDFxQall(xlha[i],Q,xf);
    cout << setprecision(1) << xlha[i] << "\t" << setprecision(4) 
	 << setw(11) << xf[2+6] - xf[-2+6] << "  "
	 << setw(11) << xf[1+6] - xf[-1+6] << "  "
	 << setw(11) << 2*(xf[-1+6] + xf[-2+6]) << "  "
	 << setw(11) << xf[4+6] + xf[-4+6] << "  "
	 << setw(11) << xf[0+6] << "  "
	 << endl;
  }
  cout << "      " << endl;
}

bool Cassandra::DataManager::LoadData( const char * fileName )
{
  // Work out which DataSet to use from the name
  string sFileName( fileName );
  Cassandra::DataManager::NodeType nodeType;
  if( sFileName.find( pszBCDMS ) != std::string::npos )
    nodeType = NodeType::BCDMS;
  else if( sFileName.find( pszJLAB ) != std::string::npos )
    nodeType = NodeType::JLAB;
  else if( sFileName.find( pszSLAC ) != std::string::npos )
    nodeType = NodeType::SLAC;
  else if( sFileName.find( pszBoNuS ) != std::string::npos )
    nodeType = NodeType::BoNuS;
  else
    return false;

  bool bRet = false;

  ifstream  fin(fileName);
  if(!fin.is_open())
    {LOG( Always, "Unable to open " << fileName << endl );}
  else
    {
      LOG( Sometimes, "Parsing data file " << fileName << endl );
      string line;
      int iLowestValidReplica;
      // This last parameter means don't error check params;
      bRet = m_DataSets[nodeType]->Load(fin, iLowestValidReplica, nullptr, nullptr,
					m_DataSets[nodeType]->m_iNumReplicas,
					-1, -1., -1., false, false);
      if( bRet )
	{
	  if( iLowestValidReplica != 0 )
	    {
	      LOG( Always, "Loaded incomplete set, Lowest Valid Replica " << iLowestValidReplica
		   << " from " << fileName << endl );
	      bRet = false;
	    }
	}

      //Close the file
      fin.close();

      m_DataSets[nodeType]->ComputeOverallStats(Cassandra::LogLevel::Rarely, 10);
      LOG( Always, "Loaded " << m_DataSets[nodeType]->m_count << " records from " << fileName << endl );
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
	  if( pDataSet->ApplyCut( W2, Q2, bStrictCut, SampleSize, bShowCut, m_TwiggySort ) )
	      pDataSet->ParkCut( SampleSize );
	}
    }
}

void Cassandra::DataManager::ParkXBelow( double X )
{
  for( DataSet * pDataSet : m_DataSets )
    if( pDataSet->size() )
      pDataSet->ParkXBelow( X );
}

void Cassandra::DataManager::ComputeOverallStats(size_t SampleSize)
{
  // Load complete - recompute statistics
  bool bFirst = true;
  m_NDat = 0;
  for( DataSet * pDataSet : m_DataSets )
    {
      pDataSet->m_list.sort( DataNodeSortFile );
      pDataSet->ComputeOverallStats(Cassandra::LogLevel::Mostly, SampleSize);
      if( pDataSet->m_list.size() > 0 )
	{
	  m_NDat += pDataSet->m_list.size();
	  if( bFirst )
	    {
	      bFirst = false;
	      m_WMinOverall = pDataSet->m_WMinOverall;
	      m_WMaxOverall = pDataSet->m_WMaxOverall;
	      m_QMinOverall = pDataSet->m_QMinOverall;
	      m_QMaxOverall = pDataSet->m_QMaxOverall;
	      m_W2MinOverall = pDataSet->m_W2MinOverall;
	      m_W2MaxOverall = pDataSet->m_W2MaxOverall;
	      m_Q2MinOverall = pDataSet->m_Q2MinOverall;
	      m_Q2MaxOverall = pDataSet->m_Q2MaxOverall;
	    }
	  else
	    {
	      if( pDataSet->m_WMinOverall < m_WMinOverall )
		m_WMinOverall = pDataSet->m_WMinOverall;
	      if( pDataSet->m_WMaxOverall > m_WMaxOverall )
		m_WMaxOverall = pDataSet->m_WMaxOverall;
	      if( pDataSet->m_QMinOverall < m_QMinOverall )
		m_QMinOverall = pDataSet->m_QMinOverall;
	      if( pDataSet->m_QMaxOverall > m_QMaxOverall )
		m_QMaxOverall = pDataSet->m_QMaxOverall;
	      if( pDataSet->m_W2MinOverall < m_W2MinOverall )
		m_W2MinOverall = pDataSet->m_W2MinOverall;
	      if( pDataSet->m_W2MaxOverall > m_W2MaxOverall )
		m_W2MaxOverall = pDataSet->m_W2MaxOverall;
	      if( pDataSet->m_Q2MinOverall < m_Q2MinOverall )
		m_Q2MinOverall = pDataSet->m_Q2MinOverall;
	      if( pDataSet->m_Q2MaxOverall > m_Q2MaxOverall )
		m_Q2MaxOverall = pDataSet->m_Q2MaxOverall;
	    }
	}
    }
  if( bFirst )
    {
      m_WMinOverall = 0;
      m_WMaxOverall = 0;
      m_QMinOverall = 0;
      m_QMaxOverall = 0;
      m_W2MinOverall = 0;
      m_W2MaxOverall = 0;
      m_Q2MinOverall = 0;
      m_Q2MaxOverall = 0;
    }
  else
    {
      LOG( Mostly, "Combined data set contains " << m_NDat
	   << " data points with statistics:" << endl
	   <<"    Item         Min         Max       Min^2       Max^2"
	   << fixed << setprecision(5) << endl
	   << "       W" << setw(12) << m_WMinOverall << setw(12) << m_WMaxOverall
	   << setw(12) << m_W2MinOverall << setw(12) << m_W2MaxOverall << endl
	   << "       Q" << setw(12) << m_QMinOverall << setw(12) << m_QMaxOverall
	   << setw(12) << m_Q2MinOverall << setw(12) << m_Q2MaxOverall << endl );
    }
}

bool Cassandra::DataManager::MakePrediction(const char * pszOutFilePrefix,
					    bool bStrictFirstCut,
					    int SampleSize, bool bParseOnly, bool bShowCut,
					    const double PW2Min, double W2Max,
					    const double PQ2Min, double Q2Max, int NumSteps )
{
  bool bMadePrediction = false;
  bool bNoError = true;

  // Get name for combined DataSet
  bool bFirst = true;
  m_pszTarget = nullptr;
  for( DataSet * pDS : m_DataSets )
    {
      if( pDS->size() )
	{
	  for( DataNode * p : pDS->m_list )
	    {
	      if( bFirst )
		{
		  bFirst = false;
		  m_pszTarget = p->Name();
		}
	      else if( m_pszTarget != p->Name() )
		{
		  static const char szCombined[] = "Combined";
		  m_pszTarget = szCombined;
		}
	    }
	}
    }
  std::string sName( m_pszTarget );
  for( DataSet * pDS : m_DataSets )
    {
      if( pDS->size() )
	{
	  sName.append( 1, ' ' );
	  sName.append( pDS->Name() );
	}
    }

  // Make a combined DataSet ... but clear it's list before it's deleted
  Cassandra::BCDMSSet dsAll( m_Params, m_DataSets[0]->m_iNumReplicas, sName.c_str() );
  for( DataSet * pDS : m_DataSets )
    for( DataNode * n : pDS->m_list )
      dsAll.m_list.push_back( n );
  dsAll.ReleaseCutInternalOnly( false );

  // Append name of dataset to file prefix
  //std::string sFileName( pszOutFilePrefix );
  //sFileName.append( 1, '_' );
  //sFileName.append( sName );

  // We might need the original parameters
  double W2Min = PW2Min;
  double Q2Min = PQ2Min;

  // Work out where the cut(s) should be based on input parameters & DataSet parameters
  double W2, W2Step;
  double Q2, Q2Step;
  bool bStrictCut = bStrictFirstCut;

  if( W2Min == -1. )
    {
      W2Min = m_W2MinOverall;
      LOG( Mostly, "Taking minimum W2 from dataset = " << setprecision(7) << W2Min << endl );
    }
  else
    { LOG( Sometimes, "Minimum W2 specified >= " << setprecision(7) << W2Min << endl ); }

  if( Q2Min == -1. )
    {
      Q2Min = m_Q2MinOverall;
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
	  W2Max = m_W2MaxOverall;
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
	  Q2Max = m_Q2MaxOverall;
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
  std::vector<CutChi> vChi;
  for( int i = 0 ; bOK && i <= NumSteps ; i++ )
    {
      if( i == 0 )
	{
	  // Give the dataset an opportunity to initialise
	  if( !bParseOnly && !InitPrediction( pszOutFilePrefix, dsAll.m_iNumReplicas ) )
	    bOK = false;
	  else
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
      if( bOK )
	{
	  // Log information about this prediction
	  const char * pszCut = Cassandra_StrictString( bStrictCut );
	  LOG( Always, i << ": W^2" << pszCut
	       << setprecision(7) << W2 << setprecision(5) << ", Q^2" << pszCut << Q2 << endl );

	  // Apply the cut to each DataSet. If any of them change, do a full prediction on combined
	  bool bDataChanged = false;
	  for( DataSet * pDataSet : m_DataSets )
	    if( pDataSet->size() && pDataSet->ApplyCut( W2, Q2, bStrictCut, SampleSize,
							bShowCut, m_TwiggySort ) )
	      bDataChanged = true;

	  // Make the prediction
	  if( bParseOnly )
	    bMadePrediction = true;
	  else if( bDataChanged )
	    {
	      if( !dsAll.ApplyCut( W2, Q2, bStrictCut, SampleSize, bShowCut, m_TwiggySort ) )
		{
		  LOG( Always, "Error: Combined dataset didn't change, but individual sets did" << endl );
		  bOK = false;
		}
	      else
		{
		  int j = 0;
		  /* This test was really only necessary for debugging
		  // Make sure we ended up with the same DataNodes in the same order
		  for( DataSet * pDS : m_DataSets )
		    for( DataNode * p : pDS->m_v )
		      if( p != dsAll.m_v[j++] )
			{
			  bOK = false;
			  LOG( Always, "*** " );
			  p->Serialise( LOG_STREAM, 5 );
			}*/
		  if( !bOK )
		    {
		      LOG( Always, i << ": Error: DataSets are in different order" << endl );
		      bNoError = false;
		    }
		  /*else if( j != ( int ) dsAll.CutSize() )
		    {
		      LOG( Always, i << ": Error: DataSets are different size" << endl );
		      bOK = false;
		    }*/
		  else
		    {
		      LOG( Sometimes, i << ": DataSets are in same order" << endl );
		      CutChi thisChi;
		      thisChi.iSeq = i;
		      thisChi.W2 = W2;
		      thisChi.Q2 = Q2;
		      thisChi.bStrict = bStrictCut;
		      bOK = MakeOnePrediction( dsAll, pszOutFilePrefix, thisChi, SampleSize, bShowCut );
		      if( bOK )
			{
			  bMadePrediction = true;
			  vChi.push_back( thisChi );
			}
		      else
			bNoError = false;
		    }
		}
	    }
	  bStrictCut = false;
	  W2 += W2Step;
	  Q2 += Q2Step;
	}
    }
  // Plot the chi table
  int iNumChi = ( int ) vChi.size();
  if( iNumChi )
    {
      // Work out whether to plot the reweighted Chi as well
      int iNumSeries = 2;
      for( int i = 0 ; i < iNumChi && iNumSeries == 2 ; i++ )
	if( vChi[i].ChiTwistReweight != 0. )
	  iNumSeries = 3;
      // Now make plots
      auto mg = new TMultiGraph();
      for( int c = 0 ; c < iNumSeries ; c++ )
	{
	  double x[iNumChi];
	  double y[iNumChi];
	  for( int i = 0 ; i < iNumChi ; i++ )
	    {
	      x[i] = vChi[i].W2;
	      switch( c )
		{
		case 0:
		  y[i] = vChi[i].Chi;
		  break;
		case 1:
		  y[i] = vChi[i].ChiTwist;
		  break;
		case 2:
		  y[i] = vChi[i].ChiTwistReweight;
		  break;
		}
	    }
	  auto gr = new TGraph( iNumChi, x, y );
	  static const char * szSeries[] = { "No HT", "+ HT", "+ HT reweight" };
	  gr->SetTitle( szSeries[c] );
	  gr->SetMarkerStyle( 23 - c );
	  gr->SetMarkerColor( 2 + c * 2 );
	  gr->SetLineColor( 2 + c * 2 );
	  mg->Add( gr );
	}
      string sFileName( pszOutFilePrefix );
      sFileName.append( ".chi.pdf" );
      LOG( Often, "Creating " << sFileName << endl );
      std::stringstream ssTitle;
      ssTitle << "Chi**2/Ndat vs W**2 " << sName << "; W**2 / GeV**2; Chi**2/Ndat";
      mg->SetTitle( ssTitle.str().c_str() );
      TCanvas *c1 = new TCanvas("c1","Unused title - was ss1",200,10,700,500);
      c1->SetGrid();
      c1->GetFrame()->SetBorderSize(12);
      //c1->SetLogx();
      mg->Draw( "ACP" );
      TLegend * pLeg = c1->BuildLegend( 0.75, 0.82, 0.9, 0.9 );
      auto OldErrLvl = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kWarning;
      c1->Print( sFileName.c_str() );
      gErrorIgnoreLevel = OldErrLvl;
      if( pLeg )
	delete pLeg;
      delete mg;
      delete c1;
    }

  // Clean up
  dsAll.m_list.clear();
  dsAll.m_l.clear();
  dsAll.m_lParked.clear();
  return bMadePrediction && bNoError;
}

int main(int argc, char * argv[])
{
  auto start_time = std::chrono::system_clock::now();

  static const char * pszOutFilePrefix = "a_out";
  static const char szDefaultPDFSet[] = "NNPDF31_nnlo_as_0118.LHgrid";

  // Give instructions if no filename supplied
  if( argc < 2 )
    {
      cout << "usage: " << argv[0] << " [switches] <list-of-filenames-to-process>" << endl
	   << "  -apx,y,z  Set model parameters, p=parameter letter (a, b, c or d)" << endl
	   << "            x=start, y=stop, [optional] x=number of intervals" << endl
	   << "  -bx[n.n]  Bin size in X    (default=0.05)" << endl
	   << "  -bq[n.n]  Bin size in Q**2 (default=0.25)" << endl
	   << "  -c        Show cut records" << endl
	   << "  -i        Input only, i.e. parse, but no predictions" << endl
	   << "  -l[n]     Log level n, 0 = terse, 4 = detailed, default 1" << endl
	   << "  -l[prefix]Log output to prefix.log. Prediction files also start with prefix" << endl
	   << "  -m[prefix]Read model from files beginning with prefix (don't compute model)" << endl
	   << "  -n[i]     number of intervals between cuts (default 0)" << endl
	   << "  -o[n]     number of outliers to cap on correlation plots (default 0)" << endl
	   << "  -p[PDF]   Show HT model impact to LHAPDF set (don't forget .LHgrid suffix)" << endl
	   << "            Default " << szDefaultPDFSet << endl
	   << "  -qc[n.n]  Throw away all data above Q^2 cut-off GeV^2 (default=3.5)" << endl
	   << "  -ql[n.n]  lower limit for Q^2 (default=data min)" << endl
	   << "  -qu[n.n]  Upper limit for Q^2 (default=data max)" << endl
	   << "  -r[n]     How many PDF Replicas to make predictions for (0 to n-1)" << endl
	   << "  -s[n]     Print n samples from each data set" << endl
	   << "  -sc       Make the first cut a strict cut, i.e. > rather than >=" << endl
	   << "  -tt[n]    Model Type: t = 'a'dditive or 'm'ultiplicative, n=0, 1 or 2" << endl
	   << "  -wc[n.n]  Throw away all data above W^2 cut-off GeV^2 (default=12.5)" << endl
	   << "  -wl[n.n]  lower limit for W^2 (default=data min)" << endl
	   << "  -wu[n.n]  Upper limit for W^2 (default=data max)" << endl
	   << "  -xc[n.n]  Throw away all data below X cut-off (default 0)" << endl
	   << "  -x[n]     Sort covariance charts by: 0 = Type, X, Q**2; 1 = Type, Q**2, X;" << endl
	   << "            2 = DataFile; (default 0)" << endl;
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
  double   dxCutoff = 0.;

  bool bStrictFirstCut = false;
  int iNumReplicas = 1;

  double xBinSize = 0.01;
  double Q2BinSize = 2;

  int iReturnValue = 0;

  ModelParams Params;
  NodeSortFunc TwiggySort = DataNodeSortTypeXQ2;
  int iNumOutliers = 3;
  const char * pszModelPrefix = nullptr;
  int iModelNum = 0;
  bool bAdditive = false;
  const char * pszPDFSet = nullptr;
  double dQChartScale = 3.;

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
			       || (argv[i][3] >= '0' && argv[i][3] <= '9') ) )
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
	      if( argv[i][2] )
		pszModelPrefix = &argv[i][2];
	      break;

	    case 'N':
	      if( argv[i][2] >= '0' && argv[i][2] <= '9' ) 
		NumSteps = atoi( &argv[i][2] );
	      break;

	    case 'O':
	      if( argv[i][2] >= '0' && argv[i][2] <= '9' )
		iNumOutliers = atoi( &argv[i][2] );
	      break;

	    case 'P':
	      {
		char c = argv[i][2];
		if( c == 0 )
		  pszPDFSet = szDefaultPDFSet;
		else if( (c >= '0' && c <= '9') || c == '.' )
		  {
		    dQChartScale = atof( &argv[i][2] );
		    // If we've specified a chart scale, clearly we want a chart!
		    if( pszPDFSet == nullptr )
		      pszPDFSet = szDefaultPDFSet;
		  }
		else
		  pszPDFSet = &argv[i][2];
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
						       || (argv[i][3] >= '0' && argv[i][3] <= '9') ) )
		    {
		      Q2Min = atof( &argv[i][3] );
		    }
		  else if( toupper( argv[i][2] ) == 'U' && ( argv[i][3] == '-' || argv[i][3] == '.'
							    || (argv[i][3] >= '0' && argv[i][3] <= '9') ) )
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
	      {
		int iOffset = 2;
		char c = toupper( argv[i][iOffset] );
		if( c == 0 )
		  {
		    APFELTest( ( pszPDFSet == nullptr ) ? szDefaultPDFSet : pszPDFSet );
		    return -5;
		  }
		if( c == 'A' )
		  {
		    bAdditive = true;
		    c = argv[i][++iOffset];
		  }
		else if( c == 'M' )
		  {
		    bAdditive = false;
		    c = argv[i][++iOffset];
		  }
		if( c >= '0' && c <= '9' )
		  iModelNum = atoi( &argv[i][iOffset] );
	      }
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
						       || (argv[i][3] >= '0' && argv[i][3] <= '9') ) )
		    {
		      W2Min = atof( &argv[i][3] );
		    }
		  else if( toupper( argv[i][2] ) == 'U' && ( argv[i][3] == '-' || argv[i][3] == '.'
							    || (argv[i][3] >= '0' && argv[i][3] <= '9') ) )
		    {
		      W2Max = atof( &argv[i][3] );
		    }
		}
	      break;

	    case 'X':
	      switch( toupper( argv[i][2] ) )
		{
		case 'C':
		    if( (argv[i][3] >= '0' && argv[i][3] <= '9') || argv[i][3] == '.' )
		    dxCutoff = atof( &argv[i][3] );
		  break;
		case '1':
		  TwiggySort = DataNodeSortTypeQ2X;
		  break;
		case '2':
		  TwiggySort = DataNodeSortFile;
		  break;
		default:
		  TwiggySort = DataNodeSortTypeXQ2;
		  break;
		}
	      break;

	    default:
	      cerr << "Ignoring unrecognised switch: " << argv[i] << endl;
	      break;
	    }
	  argv[i][0] = 0;
	}
    }

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
  LOG( Always, " GeV^2. " << endl );
  if( dxCutoff == 0. )
    { LOG( Always, "No x cutoff" ); }
  else
    { LOG( Always, "X >= " << fixed << setprecision(2) << dxCutoff ); }
  LOG( Always, endl );

  // Say what bin sizes we are using
  LOG( Always, "X bin size set to " << xBinSize << endl );
  LOG( Always, "Q**2 bin size set to " << Q2BinSize << endl );

  // Say what the model is
  if( bAdditive )
    { LOG( Always, "Additive" ); }
  else
    { LOG( Always, "Multiplicative" ); }
  LOG( Always, " model " << iModelNum << scientific
       << setprecision( std::numeric_limits<double>::digits10 )  << endl );
  for( int i = 0 ; i < 4 ; i++ )
    {
      static char szPrefix[] = "?: Start=";
      szPrefix[0] = 'A' + i;
      LOG( Always, szPrefix << Params[i].Start << ", Step=" << Params[i].Step
	   << " Num=" << Params[i].Num << endl );
    }

  // Describe whether models will be loaded or created, and which PDF set used to show impact
  if( pszModelPrefix )
    {
      LOG( Always, "Higher Twist Models will be loaded from " << pszModelPrefix << endl );
      if( pszPDFSet )
	{ LOG( Always, "Impact to PDF Set " << pszPDFSet << " will be charted" << endl ); }
      else
	{ LOG( Always, "No PDF Set specified to show impact of HT (use \"-p[PDFSet]\" to specify)" << endl );}
    }
  else
    {
      LOG( Always, "Higher Twist Models will be created from input data" << endl );
      if( pszPDFSet )
	{ LOG( Always, "PDF Set " << pszPDFSet << " specified, but not needed" << endl ); }
    }

  // Log how heatmaps will be sorted
  IF_LOG( Mostly )
  {
    LOG( Mostly, "Covariance plots will be sorted by " );
    if( TwiggySort == DataNodeSortTypeXQ2 )
      { LOG( Mostly, "DataSet, then X then Q**2" << endl ); }
    else if( TwiggySort == DataNodeSortTypeQ2X )
      { LOG( Mostly, "DataSet, then Q**2 then X" << endl ); }
    else if( TwiggySort == DataNodeSortFile )
      { LOG( Mostly, "DataSet, then original source file order" << endl ); }
    else
      { LOG( Mostly, "Error: sort undefined" << endl ); }
  }

  // Okay, finally start running our model and making predictions
  Cassandra::DataManager myPred( Params, iNumReplicas, TwiggySort, iNumOutliers, pszModelPrefix,
				 dxCutoff, iModelNum, bAdditive, pszPDFSet, dQChartScale, xBinSize, Q2BinSize);
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
      else if( dxCutoff != 0. )
	myPred.ParkXBelow( dxCutoff );
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
	  if( !myPred.MakePrediction( pszOutFilePrefix, bStrictFirstCut, SampleSize,
				      bParseOnly, bShowCut,
				      W2Min, W2Max, Q2Min, Q2Max, NumSteps ) )
	    iReturnValue = -4;
	}
    }

  auto end_time = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end_time - start_time;
  LOG( Always, "Total run time " << fixed << setprecision(1) << elapsed_seconds.count() << " seconds." << endl
       << "Return code " << iReturnValue << endl );
  return iReturnValue;
}
