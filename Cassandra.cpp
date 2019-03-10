# include "Cassandra.h"

using namespace std;
using namespace Cassandra;

// Identifies type: SLAC; JLAB; BCDMS; etc
// NB: Pointer equality guarantees string equality
const char * Cassandra::pszProton = "Proton";
const char * Cassandra::pszDeuteron = "Deuteron";
const char * Cassandra::pszBCDMS = "BCDMS";
const char * Cassandra::pszSLAC = "SLAC";
const char * Cassandra::pszJLAB = "JLAB";
const char * Cassandra::pszBoNuS = "JLAB_BoNuS";

//static Cassandra::Logger lout;

std::unique_ptr<std::ofstream> Cassandra::Global_Log_File;
Cassandra::LogLevel Cassandra::Global_Log_Level = Cassandra::LogLevel::Mostly;

// Serialise a TMatrixD to an output stream

std::ostream& Cassandra::operator<<(std::ostream& os, const TMatrixD &m)
{
  auto rowLB = m.GetRowLwb();
  auto colLB = m.GetColLwb();
  auto rowCount = m.GetNrows();
  auto colCount = m.GetNcols();
  if( rowCount > 10 )
    rowCount = 10;
  if( colCount > 10 )
    colCount = 10;

  for( auto row = 0 ; row < rowCount ; row++ )
    {
      os << "( ";
      for( auto col = 0 ; col < colCount ; col++ )
	{
	  os << m[row + rowLB][col + colLB] << " ";
	}
      os << ")" << endl;
    }
  return os;
}

// Serialise a Datum to an output stream

std::ostream& Cassandra::operator<<(std::ostream& os, const Datum& d)
{
  os << std::setw(9) << d.m_Value << std::setw(9) << d.m_errStat << std::setw(9) << d.m_errSys;
  return os;
}

// String describing strictness of cut

const char * Cassandra::Cassandra_StrictString( bool bStrictCut )
{
  const char * pszCut;
  if( bStrictCut )
    pszCut = " > ";
  else
    pszCut = " >= ";
  return pszCut;
}

// Helper functions to sort lists of DataNode

// File order

bool Cassandra::DataNodeSortFile( const DataNode * p1, const DataNode * p2 )
{
  if( p1->Type() < p2->Type() )
    return true;
  if( p1->Type() > p2->Type() )
    return false;
  if( p1->PredictionTargets() < p2->PredictionTargets() )
    return true;
  if( p1->PredictionTargets() > p2->PredictionTargets() )
    return false;
  if( p1->m_eBeam < p2->m_eBeam )
    return true;
  if( p1->m_eBeam > p2->m_eBeam )
    return false;
  if( p1->m_i < p2->m_i )
    return true;
  return false;
}

// Sort by Q and X first, then file order

bool Cassandra::DataNodeSortQXFile( const DataNode * p1, const DataNode * p2 )
{
  if( p1->m_cQ < p2->m_cQ )
    return true;
  if( p1->m_cQ > p2->m_cQ )
    return false;
  if( p1->m_x < p2->m_x )
    return true;
  if( p1->m_x > p2->m_x )
    return false;
  return DataNodeSortFile( p1, p2 );
}

// Sort each file by W, then file order

bool Cassandra::DataNodeSortW( const DataNode * p1, const DataNode * p2 )
{
  if( p1->Type() < p2->Type() )
    return true;
  if( p1->Type() > p2->Type() )
    return false;
  if( p1->PredictionTargets() < p2->PredictionTargets() )
    return true;
  if( p1->PredictionTargets() > p2->PredictionTargets() )
    return false;
  if( p1->m_cW < p2->m_cW )
    return true;
  if( p1->m_cW > p2->m_cW )
    return false;
  return DataNodeSortFile( p1, p2 );
}

// Sort by Q2 only
bool Cassandra::DataNodeSortQ2Only( const DataNode * p1, const DataNode * p2 )
{
  return p1->m_Q2 < p2->m_Q2;
}

// Sort by X only
bool Cassandra::DataNodeSortXOnly( const DataNode * p1, const DataNode * p2 )
{
  return p1->m_x < p2->m_x;
}

// Sort by DataSet type, then Q2, then X
bool Cassandra::DataNodeSortTypeQ2X( const DataNode * p1, const DataNode * p2 )
{
  if( p1->Type() < p2->Type() )
    return true;
  if( p1->Type() > p2->Type() )
    return false;
  if( p1->PredictionTargets() < p2->PredictionTargets() )
    return true;
  if( p1->PredictionTargets() > p2->PredictionTargets() )
    return false;
  if( p1->m_Q2 < p2->m_Q2 )
    return true;
  if( p1->m_Q2 > p2->m_Q2 )
    return false;
  if( p1->m_x < p2->m_x )
    return true;
  if( p1->m_x > p2->m_x )
    return false;
  return true;
}

// Sort by DataSet type, then X, then Q2
bool Cassandra::DataNodeSortTypeXQ2( const DataNode * p1, const DataNode * p2 )
{
  if( p1->Type() < p2->Type() )
    return true;
  if( p1->Type() > p2->Type() )
    return false;
  if( p1->PredictionTargets() < p2->PredictionTargets() )
    return true;
  if( p1->PredictionTargets() > p2->PredictionTargets() )
    return false;
  if( p1->m_x < p2->m_x )
    return true;
  if( p1->m_x > p2->m_x )
    return false;
  if( p1->m_Q2 < p2->m_Q2 )
    return true;
  if( p1->m_Q2 > p2->m_Q2 )
    return false;
  return true;
}

void Cassandra::StartStep::FromTo( double dFrom, double dTo, int iStep )
{
  Start = dFrom;
  if( dFrom == dTo || iStep <= 1 )
    {
      Num = 1;
      Step = 0.; // Not that this will matter
    }
  else
    {
      Num = abs(iStep);
      Step = (dTo - dFrom) / ( Num - 1 );
    }
}

Cassandra::ModelParams::ModelParams(int iNumParams)
  : m_vSS(iNumParams)
{
  if( iNumParams >= 0 )
    {
      m_vSS[0].Start = 0.05;
      m_vSS[0].Step = 0.1;
      if( iNumParams >= 1 )
	{
	  m_vSS[1].Start = 1.;
	  m_vSS[1].Step = 1.;
	  if( iNumParams >= 2 )
	    {
	      m_vSS[2].Start = -1.8;
	      m_vSS[2].Step = 0.4;
	      if( iNumParams >= 3 )
		{
		  m_vSS[3].Start = -1.8;
		  m_vSS[3].Step = 0.4;
		}
	    }
	}
    }
}

Cassandra::ModelParams::ModelParams()
  : ModelParams(4)
{
}

void Cassandra::DataNode::SerialiseFieldName(FieldType Type, int iIndex, std::ostream& os, int iLen) const
{
  bool bAddSuffix = false;
  int iSuffix;
  const std::string& s = FieldName(Type, iIndex, bAddSuffix, iSuffix);
  if( !bAddSuffix )
    os << setw(iLen) << s;
  else
    {
      std::string sNum = std::to_string( iSuffix );
      int iPad = iLen - s.length() - sNum.length();
      if( iPad > 0 )
	os << setw(iPad + s.length()) << s << sNum;
      else
	os << s << sNum;
    }
  if( Type == FieldType::Data )
    os << "  errStat   errSys";
}

// Serialise a DataNode header to an output stream

void Cassandra::DataNode::SerialiseHeader( std::ostream & os, int iMaxTheory, int iMinTheory ) const
{
  if( iMaxTheory == std::numeric_limits<int>::max() )
      os << "Nucleon eBeam i x Q**2 Q W**2 W y";
  else
    os << "  eBeam       i        x     Q**2        Q     W**2        W        y";
  for( size_t i = 0 ; i < m_Data.size() ; i ++ )
    SerialiseFieldName(Cassandra::DataNode::FieldType::Data, i, os, 9);
  for( int i = iMinTheory ; i < iMaxTheory && i < ( int ) m_Theory.size() ; i ++ )
    SerialiseFieldName(Cassandra::DataNode::FieldType::Theory, i, os, 10);
}

// Serialise a DataNode to an output stream
void Cassandra::DataNode::Serialise( std::ostream & os, int iMaxTheory, int iMinTheory ) const
{
  const char sp = ' ';
  bool bHuman;
  if( iMaxTheory == std::numeric_limits<int>::max() )
    bHuman = false;
  else
    bHuman = true;
  if( bHuman )
    {
      // Original input
      os << fixed
	// Beam energy and index
	 << setw(7) << setprecision (0) << m_eBeam
	 << setw(8) << /*setprecision (0) <<*/ m_i
	// x and Q^2
	 << setw(9) << setprecision(5) << m_x
	 << setw(9) << setprecision(2) << m_Q2
	// Additional fields
	 << setw(9) << setprecision(2) << m_cQ
	 << setw(9) << /*setprecision(2) <<*/ m_cW2
	 << setw(9) << /*setprecision(2) <<*/ m_cW
	 << setw(9) << setprecision(5) << m_cy;
      for( Datum d : m_Data )
	os << d;
      // Theoretical predictions
      os << fixed << setprecision(5);
      for( int i = iMinTheory ; i < iMaxTheory && i < ( int ) m_Theory.size() ; i ++ )
	os << setw(10) << m_Theory[i];
    }
  else
    {
      os << Name() << scientific << setprecision(std::numeric_limits<double>::digits10)
	 << sp << m_eBeam << sp << m_i << sp << m_x << sp << m_Q2 << sp << m_cQ
	 << sp << m_cW2 << sp << m_cW << sp << m_cy;

      for( Datum d : m_Data )
	os << sp << d.m_Value << sp << d.m_errStat << sp << d.m_errSys;

      os << scientific << setprecision(std::numeric_limits<double>::digits10);
      for( int i = iMinTheory ; i < iMaxTheory && i < ( int ) m_Theory.size() ; i ++ )
	os << sp << m_Theory[i];
    }
}

const std::string& Cassandra::DataNodeSingle::FieldName( Cassandra::DataNode::FieldType t, int iIndex,
							 bool &bAddSuffix, int &iSuffix ) const
{
  static const std::string s("F2_Data");
  if( t == Cassandra::DataNode::FieldType::Data )
    {
      bAddSuffix = false;
      return s;
    }
  static const std::string sT("F2T");
  iSuffix = iIndex;
  bAddSuffix = true;
  return sT;
}

void Cassandra::NodeProton::MakePrediction(Target t, int iReplica)
{
  switch(t)
    {
    case Proton:
      m_Theory[iReplica] = APFEL::F2total(m_x);
      break;
    default:
      break;
    }
}

void Cassandra::NodeDeuteron::MakePrediction(Target t, int iReplica)
{
  switch(t)
    {
    case Deuteron:
      m_Theory[iReplica] = APFEL::F2total(m_x);
      break;
      /*case Proton:
      m_Theory[1] = APFEL::F2total(m_x);
      break;
    case Neutron:
      m_Theory[2] = APFEL::F2total(m_x);
      break;
    case None:
      m_Theory[3] = (m_Theory[1] + m_Theory[2]) / 2.0;
      m_Theory[4] = m_Theory[3] - m_Theory[0];
      break;*/
    default:
      break;
    }
}

double Cassandra::BCDMSNodeProton::Covariance(const DataNode &nr, bool bDiagonal) const
{
  const std::type_info & tiNR = typeid( nr );
  bool bSame;
  if( tiNR == typeid( BCDMSNodeProton ) )
    bSame = true;
  else if( tiNR == typeid( BCDMSNodeDeuteron ) )
    bSame = false;
  else
    return 0;

  // Construct total relative (multiplicative) normalisation error
  // by adding all the separate sources of error in quadrature
  // i.e. sigma_N2 is the square of the total relative normalisation error
  // See equation 2, pg 5 of http://www.inspirehep.net/record/585628
  // and addition in quadrature immediately below
  // I.e. (sigma_N)^2 = sigma_Na^2 + sigma_Nt^2 + sigma_Nb^2

  double sigma_N2 = 0.03 * 0.03; // sigma_Na, absolute normalisation uncertainty
  if( !bSame )
    {
      // sigma_Nt, relative normalisation for different targets
      sigma_N2 += 0.02 * 0.02;
    }
  else if(m_eBeam != nr.m_eBeam)
    {
      // sigma_Nb, relative normalisation, same target, different beam energies
      // 1.5% between deuterons if one has beam energy 280 GeV
      // 1% for all others
      sigma_N2 += 0.01 * 0.01;
    }

  // Now that we have sigma_N squared, build covariance
  // Assuming that errSys is total error due to:
  //     f_b incoming muon (beam) energy
  //     f_s outgoing muon energy (spectrometer magnetic field)
  //     f_r spectrometer resolution
  // The paper describes these as fully correlated and gives rationale at top pg 5
  double covar = ScaledErrSys() * nr.ScaledErrSys() + sigma_N2 * m_Theory[0] * nr.m_Theory[0];

  // Add in diagonal components
  if( bDiagonal )
    covar += ScaledErrStat() * nr.ScaledErrStat();

  // That's our covariance built
  return covar;
}

double Cassandra::BCDMSNodeDeuteron::Covariance(const DataNode &nr, bool bDiagonal) const
{
  const std::type_info & tiNR = typeid( nr );
  bool bSame;
  if( tiNR == typeid( BCDMSNodeDeuteron ) )
    bSame = true;
  else if( tiNR == typeid( BCDMSNodeProton ) )
    bSame = false;
  else
    return 0;

  // Construct total relative (multiplicative) normalisation error
  // by adding all the separate sources of error in quadrature
  // i.e. sigma_N2 is the square of the total relative normalisation error
  // See equation 2, pg 5 of http://www.inspirehep.net/record/585628
  // and addition in quadrature immediately below
  // I.e. (sigma_N)^2 = sigma_Na^2 + sigma_Nt^2 + sigma_Nb^2

  double sigma_N2 = 0.03 * 0.03; // sigma_Na, absolute normalisation uncertainty
  if( !bSame )
    {
      // sigma_Nt, relative normalisation for different targets
      sigma_N2 += 0.02 * 0.02;
    }
  else if(m_eBeam != nr.m_eBeam)
    {
      // sigma_Nb, relative normalisation, same target, different beam energies
      // 1.5% between deuterons if one has beam energy 280 GeV
      // 1% for all others
      if( m_eBeam == 280. || nr.m_eBeam == 280. )
	sigma_N2 += 0.015 * 0.015;
      else
	sigma_N2 += 0.01 * 0.01;
    }

  // Now that we have sigma_N squared, build covariance
  // Assuming that errSys is total error due to:
  //     f_b incoming muon (beam) energy
  //     f_s outgoing muon energy (spectrometer magnetic field)
  //     f_r spectrometer resolution
  // The paper describes these as fully correlated and gives rationale at top pg 5
  double covar = ScaledErrSys() * nr.ScaledErrSys() + sigma_N2 * m_Theory[0] * nr.m_Theory[0];

  // Add in diagonal components
  if( bDiagonal )
    covar += ScaledErrStat() * nr.ScaledErrStat();

  // That's our covariance built
  return covar;
}

// Covariance data for JLAB comes from Physical Review C 80 035207 (2009)
// Total normalisation uncertainty 1.75%, last paragraph left hand column, pg 035207-14
// which gives the total value for all the errors in table III

double Cassandra::JLABNodeProton::Covariance(const DataNode &nr, bool bDiagonal) const
{
  const std::type_info & tiNR = typeid( nr );
  if( tiNR != typeid( JLABNodeProton ) && tiNR != typeid( JLABNodeDeuteron ) )
    return 0;

  // Construct total normalisation error
  double sigma_N2 = 0.0175 * 0.0175; // total normalisation error

  // Now that we have sigma_N squared, build covariance
  double covar = ScaledErrSys() * nr.ScaledErrSys() + sigma_N2 * m_Theory[0] * nr.m_Theory[0];

  // Add in diagonal components
  if( bDiagonal )
    covar += ScaledErrStat() * nr.ScaledErrStat();

  // That's our covariance built
  return covar;
}

double Cassandra::JLABNodeDeuteron::Covariance(const DataNode &nr, bool bDiagonal) const
{
  const std::type_info & tiNR = typeid( nr );
  if( tiNR != typeid( JLABNodeProton ) && tiNR != typeid( JLABNodeDeuteron ) )
    return 0;

  // Construct total normalisation error
  double sigma_N2 = 0.0175 * 0.0175; // total normalisation error

  // Now that we have sigma_N squared, build covariance
  double covar = ScaledErrSys() * nr.ScaledErrSys() + sigma_N2 * m_Theory[0] * nr.m_Theory[0];

  // Add in diagonal components
  if( bDiagonal )
    covar += ScaledErrStat() * nr.ScaledErrStat();

  // That's our covariance built
  return covar;
}

double Cassandra::SLACNodeProton::Covariance(const DataNode &nr, bool bDiagonal) const
{
  const std::type_info & tiNR = typeid( nr );
  bool bSame;
  if( tiNR == typeid( SLACNodeProton ) )
      bSame = true;
  else if( tiNR == typeid( SLACNodeDeuteron ) )
      bSame = false;
  else
      return 0;

  // Reference for covariance matrix construction:
  //   L. W. Whitlow, E. M. Riordan, S. Dasu, S. Rock, and A. Bodek,
  //   Precise measurements of the proton and deuteron structure functions
  //   from a global analysis of the SLAC deep inelastic electron scattering cross-sections
  //   Phys. Lett. B282 (1992) 475–482.
  // cited as reference 32 in NNPDF3.1 paper, https://arxiv.org/pdf/1706.00428
  // pg 476, right column, 2nd paragraph
  //   overall normalisation
  //     2.1% for proton
  //     1.7% for deuteron
  //  1.1% relative normalisation
  // Following the lead of BCDMS code, we assume the systematic errors
  // are uncorrelated, so we add them on the diagonal

  // absolute normalisation
  double rel_error = 0.021;
  if( bSame )
    rel_error *= 0.021;
  else
    rel_error *= 0.017;

  // relative (between targets) normalisation
  if( !bSame )
    rel_error += 0.011 * 0.011;

  // Now that we have sigma_N squared, build covariance
  double covar = rel_error * m_Theory[0] * nr.m_Theory[0];

  // Add in diagonal components
  if( bDiagonal )
    covar += ScaledErrSys() * nr.ScaledErrSys() + ScaledErrStat() * nr.ScaledErrStat();

  // That's our covariance built
  return covar;
}

double Cassandra::SLACNodeDeuteron::Covariance(const DataNode &nr, bool bDiagonal) const
{
  const std::type_info & tiNR = typeid( nr );
  bool bSame;
  if( tiNR == typeid( SLACNodeDeuteron ) )
      bSame = true;
  else if( tiNR == typeid( SLACNodeProton ) )
      bSame = false;
  else
      return 0;

  // Reference for covariance matrix construction:
  //   L. W. Whitlow, E. M. Riordan, S. Dasu, S. Rock, and A. Bodek,
  //   Precise measurements of the proton and deuteron structure functions
  //   from a global analysis of the SLAC deep inelastic electron scattering cross-sections
  //   Phys. Lett. B282 (1992) 475–482.
  // cited as reference 32 in NNPDF3.1 paper, https://arxiv.org/pdf/1706.00428
  // pg 476, right column, 2nd paragraph
  //   overall normalisation
  //     2.1% for proton
  //     1.7% for deuteron
  //  1.1% relative normalisation
  // Following the lead of BCDMS code, we assume the systematic errors
  // are uncorrelated, so we add them on the diagonal

  // absolute normalisation
  double rel_error = 0.017;
  if( bSame )
    rel_error *= 0.017;
  else
    rel_error *= 0.021;

  // relative (between targets) normalisation
  if( !bSame )
    rel_error += 0.011 * 0.011;

  // Now that we have sigma_N squared, build covariance
  double covar = rel_error * m_Theory[0] * nr.m_Theory[0];

  // Add in diagonal components
  if( bDiagonal )
    covar += ScaledErrSys() * nr.ScaledErrSys() + ScaledErrStat() * nr.ScaledErrStat();

  // That's our covariance built
  return covar;
}

double Cassandra::BoNuSNode::Covariance(const DataNode &nr, bool bDiagonal) const
{
  const std::type_info & tiNR = typeid( nr );
  if( tiNR != typeid( BoNuSNode ) )
    return 0;

  // Construct total normalisation error
  double sigma_N2 = 0.03 * 0.03; // sigma_Na, global normalisation error

  // Now that we have sigma_N squared, build covariance
  double covar = ScaledErrSys() * nr.ScaledErrSys() + sigma_N2 * m_Theory[0] * nr.m_Theory[0];

  // Add in diagonal components
  if( bDiagonal )
    covar += ScaledErrStat() * nr.ScaledErrStat();

  // That's our covariance built
  return covar;
}

/*Cassandra::BoNuSNode::BoNuSNode()
{
  m_ratio = 0.0;
  m_errStat = 0.0;
  m_errSys = 0.0;
  m_cF2n = 0.0;
  m_cF2d = 0.0;
  m_cRatio = 0.0;
  }*/

/*Cassandra::BoNuSNode::BoNuSNode(const BoNuSNode &r)
{
  m_ratio = r.m_ratio;
  m_errStat = r.m_errStat;
  m_errSys = r.m_errSys;
  m_cF2n = r.m_cF2n;
  m_cF2d = r.m_cF2d;
  m_cRatio = r.m_cRatio;
  }*/

// Serialise a BoNuSNode to an output stream

/*void BoNuSNode::Serialise( std::ostream & os ) const
{
  // Serialise base class
  DataNode::Serialise(os);

  os << setw(6) << setprecision (0) << m_eBeam
     << setw(8) <<  m_i
    // Additional fields
     << setprecision(5) << setw(9) << m_cQ << setw(8) << m_cy;
}*/

const std::string& Cassandra::BoNuSNode::FieldName( Cassandra::DataNode::FieldType t, int iIndex,
						   bool &/*bAddSuffix*/, int &/*iSuffix*/ ) const
{
  static const std::string s("f2N/f2D");
  if( t == Cassandra::DataNode::FieldType::Theory )
    {
      static const std::string st0("Tf2N/f2D");
      static const std::string st1("TF2N");
      static const std::string st2("TF2D");
      switch( iIndex )
	{
	case 0:
	  return st0;
	case 1:
	  return st1;
	return st2;
	}
    }
  return s;
}

void Cassandra::BoNuSNode::MakePrediction(Target /*t*/, int /*iReplica*/)
{
}

// Reconstipate a DataNode from a stringstream

Cassandra::DataNode * Cassandra::DataSet::NewDataNode( std::stringstream & ss, int iLowestValidReplica,
						       int iNumReplicas ) const
{
  string sNodeType;
  ss >> sNodeType;
  Cassandra::DataNode::Target tNodeType;
  if( sNodeType == pszProton )
    tNodeType = Cassandra::DataNode::Target::Proton;
  else if( sNodeType == pszDeuteron )
    tNodeType = Cassandra::DataNode::Target::Deuteron;
  else
    {
      LOG( Always, "Cassandra file contains unknown node type " << sNodeType << endl );
      return nullptr;
    }

  int i;
  int iIn;
  double eBeam, dIn[9];
  ss >> eBeam >> iIn;
  for( i = 0 ; i < 9 ; i++ )
    ss >> dIn[i];
  Cassandra::DataNode * p = NewDataNode( tNodeType, eBeam, iIn, dIn[0], dIn[1],
					 Datum( dIn[6], dIn[7], dIn[8] ), iNumReplicas );

  // I'm checking the calculated values and logging if they disagree with file
  bool bOK = true;
  if( abs( p->m_cQ - dIn[2] ) > 0.001 )
    {
      bOK = false;
      LOG( Always, "Cassandra file error: Q (file) = " << dIn[2]
	   << ", m_cQ = " << p->m_cQ << endl );
    }
  if( abs( p->m_cW2 - dIn[3] ) > 0.001 )
    {
      bOK = false;
      LOG( Always, "Cassandra file error: W**2 (file) = " << dIn[3]
	   << ", m_cW2 = " << p->m_cW2 << endl );
    }
  if( abs( p->m_cW - dIn[4] ) > 0.001 )
    {
      bOK = false;
      LOG( Always, "Cassandra file error: W (file) = " << dIn[4]
	   << ", m_cW = " << p->m_cW << endl );
    }
  if( abs( p->m_cy - dIn[5] ) > 0.001 )
    {
      bOK = false;
      LOG( Always, "Cassandra file error: y (file) = " << dIn[5]
	   << ", m_cy = " << p->m_cy << endl );
    }

  while( iLowestValidReplica++ < iNumReplicas )
    {
      double d;
      ss >> d;
      p->m_Theory[iLowestValidReplica - 1] = d;
    }

  // If we had any errors, then delete this node and return error
  // ... comment this out if you don't want these to be errors (eg if changing calculations)
  if( !bOK )
    {
      delete p;
      p = nullptr;
    }
  return p;
}

// Reconstipate DataSet from stream (file)
// If load successful, iLowestValidReplica contains the lowest valid replica in file

bool Cassandra::DataSet::Load(std::istream &is, int &iLowestValidReplica,
			      const char * pszPDFSet, const char * pszMassScheme,
			      int iNumReplicas, int iCount, double W2, double Q2, bool bStrict,
			      bool bValidateParameters)
{
  std::string line;
  bool bRet = true;
  int iPhase = 0;
  m_count = 0;

  int iCheckCount = 0; // How many of the parameters were validated
  bool bGotPDFSet = false;
  bool bGotMassScheme = false;
  bool bGotNumReplicas = false;
  bool bGotCount = false;
  bool bGotW2 = false;
  bool bGotQ2 = false;
  bool bGotStrict = false;

  // Unless specified in file, all replica predictions are valid
  int iLowestValidReplicaFile = 0;

  // Read the next line from the file
  while( bRet && getline(is,line) )
    {
      // Ignore blank lines
      if( line.find_first_not_of( " \t") != string::npos )
	{
	  std::stringstream ss(line);

	  // Try to identify the file type from the header
	  if( iPhase <= 1 )
	    {
	      LOG( Rarely, line << endl ); // Print the header of each file until identified

	      string Tokens[4];
	      for( string &s : Tokens )
		ss >> s;

	      if( iPhase == 0 )
		{
		  bRet = false;
		  if( Tokens[0] != "Cassandra" || Tokens[1] != "Prediction" )
		    { LOG( Always, "Not a Cassandra Prediction file" << endl ); }
		  else if( Tokens[3] != "F2" )
		    { LOG( Always, "Cassandra file contains predictions for " << Tokens[3]
			   << " not F2" << endl ); }
		  else if( Tokens[2] != Name() )
		    { LOG( Always, "Cassandra file contains predictions for " << Tokens[2]
			   << " not " << Name() << endl ); }
		  else
		    {
		      LOG( Always, "Cassandra file contains predictions for " << Tokens[2] << endl );
		      bRet = true;
		      iPhase++;
		    }
		}
	      else if( Tokens[1] == ":" )
		{
		  if( Tokens[0] == "PDFSet" )
		    {
		      if( !bValidateParameters || Tokens[2] == pszPDFSet )
			{
			  if( !bGotPDFSet )
			    {
			      bGotPDFSet = true;
			      iCheckCount++;
			    }
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file contains predictions for " << Tokens[2]
			       << " not " << pszPDFSet << endl );
			}
		    }
		  else if( Tokens[0] == "MassScheme" )
		    {
		      if( !bValidateParameters || Tokens[2] == pszMassScheme )
			{
			  if( !bGotMassScheme )
			    {
			      bGotMassScheme = true;
			      iCheckCount++;
			    }
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file contains predictions for " << Tokens[2]
			       << " not " << pszMassScheme << endl );
			}
		    }
		  else if( Tokens[0] == "PDFReplicas" )
		    {
		      if( Tokens[2][0] >= '0' && Tokens[2][0] <= '9'
			  && ( m_iNumReplicas = atoi(Tokens[2].c_str()) ) >= iNumReplicas )
			{
			  if( !bGotNumReplicas )
			    {
			      bGotNumReplicas = true;
			      iCheckCount++;
			    }
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file contains predictions for " << m_iNumReplicas
			       << " replicas, but we need " << iNumReplicas << endl );
			}
		    }
		  else if( Tokens[0] == "Count" )
		    {
		      int iFileCount = -99;
		      if( !bValidateParameters || ( Tokens[2][0] >= '0' && Tokens[2][0] <= '9'
			  && ( iFileCount = atoi(Tokens[2].c_str()) ) == iCount ) )
			{
			  if( !bGotCount )
			    {
			      bGotCount = true;
			      iCheckCount++;
			    }
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file contains predictions for " << iFileCount
			       << " data points, but we need " << iCount << endl );
			}
		    }
		  else if( Tokens[0] == "W**2" )
		    {
		      double FileW2 = -9.9;
		      char c = Tokens[2][0];
		      if( !bValidateParameters || ( ( c == '+' || c == '-' || ( c >= '0' && c <= '9' ) )
			  && abs( ( FileW2 = atof(Tokens[2].c_str()) ) - W2 ) < 0.001 ) )
			{
			  if( !bGotW2 )
			    {
			      bGotW2 = true;
			      iCheckCount++;
			    }
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file contains predictions for W**2 cut at "
			       << scientific << setprecision(std::numeric_limits<double>::digits10)
			       << FileW2 << ", but we need " << W2 << endl );
			}
		    }
		  else if( Tokens[0] == "Q**2" )
		    {
		      double FileQ2=-8.8;
		      char c = Tokens[2][0];
		      if( !bValidateParameters || ( ( c == '+' || c == '-' || ( c >= '0' && c <= '9' ) )
			  && abs( ( FileQ2 = atof(Tokens[2].c_str()) ) - Q2 ) < 0.001 ) )
			{
			  if( !bGotQ2 )
			    {
			      bGotQ2 = true;
			      iCheckCount++;
			    }
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file contains predictions for Q**2 cut at "
			       << scientific << setprecision(std::numeric_limits<double>::digits10)
			       << FileQ2 << ", but we need " << Q2 << endl );
			}
		    }
		  else if( Tokens[0] == "Cut" )
		    {
		      //double FileQ2;
		      //char c;
		      if( !bValidateParameters || ( Tokens[2][0] == '>' &&
			  ( ( bStrict && Tokens[2][1] != '=' ) || ( !bStrict && Tokens[2][1] == '=' ) ) ) )
			{
			  if( !bGotStrict )
			    {
			      bGotStrict = true;
			      iCheckCount++;
			    }
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file contains predictions for " << Tokens[2][0]
			       << " cut, but should be" << Cassandra_StrictString( bStrict ) << endl );
			}
		    }
		  else if( Tokens[0] == "LowestValidReplica" )
		    {
		      if( Tokens[2][0] >= '0' && Tokens[2][0] <= '9'
			  && ( iLowestValidReplicaFile = atoi(Tokens[2].c_str()) ) < m_iNumReplicas )
			{
			}
		      else
			{
			  bRet = false;
			  LOG( Always, "Cassandra file LowestValidReplica=" << Tokens[2][0]
			       << ", but should be < m_iNumReplicas, i.e. < " << m_iNumReplicas << endl );
			}
		    }
		}
	      else if( Tokens[0] == "Nucleon" && Tokens[1] == "eBeam"
		       && Tokens[2] == "i" && Tokens[3] == "x" )
		{
		  // Make sure we've extracted all the info we need from the header
		  if( iCheckCount != 7 )
		    {
		      bRet = false;
		      LOG( Always, "Old Cassandra file (missing prediction parameters)" << endl );
		    }
		  else
		    {
		      // Got what we need, so parse individual records
		      iPhase++;
		      LOG( Always, scientific << setprecision(std::numeric_limits<double>::digits10) );
		    }
		}
	    }
	  else
	    {
	      // Reconstipate each node
	      DataNode * p = NewDataNode( ss, iLowestValidReplicaFile, m_iNumReplicas );
	      if( p == nullptr )
		bRet = false;
	      else
		{
		  m_list.push_back(p);
		  m_count++;
		}
	    }
	}
    }
  // Ensure we loaded the right amount of data
  if( bRet && bValidateParameters && ( int ) m_list.size() != iCount )
    {
      bRet = false;
      LOG( Always, "Cassandra file header declares " << iCount
	   << " data points, but actually contains " << m_list.size() << endl );
    }
  // If we fail, clear the list
  if( !bRet )
    EmptyList();
  else
    // Only update lowest valid replica on success
    iLowestValidReplica = iLowestValidReplicaFile;
  return bRet;
}

void Cassandra::DataSet::EmptyList(void)
{
  for( list<DataNode *>::iterator i = m_list.begin() ; i != m_list.end() ; )
    {
      DataNode * p = *i;
      i = m_list.erase(i);
      delete p;
    }
  for( list<DataNode *>::iterator i = m_l.begin() ; i != m_l.end() ; )
    {
      DataNode * p = *i;
      i = m_l.erase(i);
      delete p;
    }
  for( list<DataNode *>::iterator i = m_lParked.begin() ; i != m_lParked.end() ; )
    {
      DataNode * p = *i;
      i = m_lParked.erase(i);
      delete p;
    }
  ReleaseCutInternalOnly(false);
}

void Cassandra::DataSet::NewFile()
{
  m_stage = 0;
  m_count = 0;
  ParseBegin();
}

void Cassandra::DataSet::DumpListUncoolHelper(size_t sample, int &j, int &LastNodeCount,
		            bool & bHeaderPrinted, DataNode * & pLastNode, DataNode * p)
{
  if( pLastNode && typeid( *p ) != typeid( *pLastNode ) )
    {
      LOG( Mostly, Name() << " " << pLastNode->Name() << " contains "
	   << LastNodeCount << " data points." << endl );
      LastNodeCount = 0;
      bHeaderPrinted = false;
    }
  if( ( j % sample ) == 0 )
    {
      IF_LOG( Often )
      {
	if( !bHeaderPrinted )
	  {
	    LOG( Often, Name() << " " << p->Name() << endl );
	    p->SerialiseHeader( LOG_STREAM, 5 );
	    LOG( Often, endl );
	    bHeaderPrinted = true;
	  }
	p->Serialise( LOG_STREAM, 5 );
	LOG( Often, endl );
      }
    }
  pLastNode = p;
  LastNodeCount++;
  j++;
}

void Cassandra::DataSet::DumpList(size_t SampleSize)
{
  if( SampleSize != 0 ) IF_LOG( Mostly )
    {
      size_t sample = CutSize();
      if( sample == 0 )
	sample = m_list.size();
      if( sample > 0 )
	{
	  LOG( Mostly, Name() << " contains " << sample << " data points" );
	  IF_LOG( Often )
	  {
	    LOG( Often, ". Printing " );
	    if( ( int ) SampleSize < 0 || sample <= 2 * SampleSize )
	      sample = 1;
	    else
	      sample = sample / SampleSize;
	    if( sample == 1 )
	      LOG( Often, "every record" )
	    else
	      LOG( Often, "1 record in " << sample );
	  }
	  LOG( Mostly, "." << endl );

	  // Print each node. Formatting makes this look like original input plus extra fields
	  int j = 0;
	  int LastNodeCount = 0;
	  bool bHeaderPrinted = false;
	  DataNode * pLastNode = nullptr;
	  if( CutSize() )
	    {
	      for( DataNode * p : m_v )
		DumpListUncoolHelper(sample, j, LastNodeCount, bHeaderPrinted, pLastNode, p);
	    }
	  else
	    {
	      for( DataNode * p : m_list )
		DumpListUncoolHelper(sample, j, LastNodeCount, bHeaderPrinted, pLastNode, p);
	    }
	  if( !pLastNode )
	    LOG( Always, "Fark!!" << endl )
	  else
	    LOG( Mostly, Name() << " " << pLastNode->Name() << " contains "
		 << LastNodeCount << " data points." << endl );
	}
    }
}

void Cassandra::DataSet::ReleaseCutInternalOnly(bool bCutAll)
{
  m_v.clear();
  m_TargetsNeeded = 0;
  m_WMin = 0.0;
  m_WMax = 0.0;
  m_QMin = 0.0;
  m_QMax = 0.0;
  m_MCovar.ResizeTo( 1, 1 );
  m_MCovarInv.ResizeTo( 1, 1 );
  m_bCovarBuilt = false;

  // Cut everything on our short-list
  if( bCutAll )
    {
      for( list<DataNode *>::iterator i = m_l.begin() ; i != m_l.end() ; )
	{
	  DataNode * p = * i;
	  i = m_l.erase(i);
	  m_list.push_back(p);
	}
    }

  // The cut list should be sorted in file order
  m_list.sort( DataNodeSortFile );
}

// Compute stats for overall list - don't change unless we park a cut

void Cassandra::DataSet::ComputeOverallStats(Cassandra::LogLevel MsgLevel, size_t SampleSize)
{
  m_TargetsNeededOverall = 0;
  if( m_list.size() == 0 )
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
      bool bFirst = true;
      for( DataNode * p : m_list )
	{
	  if( bFirst )
	    {
	      bFirst = false;
	      m_WMinOverall = p->m_cW;
	      m_WMaxOverall = p->m_cW;
	      m_QMinOverall = p->m_cQ;
	      m_QMaxOverall = p->m_cQ;
	      m_W2MinOverall = p->m_cW2;
	      m_W2MaxOverall = p->m_cW2;
	      m_Q2MinOverall = p->m_Q2;
	      m_Q2MaxOverall = p->m_Q2;
	    }
	  else
	    {
	      if( p->m_cW < m_WMinOverall )
		m_WMinOverall = p->m_cW;
	      if( p->m_cW > m_WMaxOverall )
		m_WMaxOverall = p->m_cW;
	      if( p->m_cQ < m_QMinOverall )
		m_QMinOverall = p->m_cQ;
	      if( p->m_cQ > m_QMaxOverall )
		m_QMaxOverall = p->m_cQ;
	      if( p->m_cW2 < m_W2MinOverall )
		m_W2MinOverall = p->m_cW2;
	      if( p->m_cW2 > m_W2MaxOverall )
		m_W2MaxOverall = p->m_cW2;
	      if( p->m_Q2 < m_Q2MinOverall )
		m_Q2MinOverall = p->m_Q2;
	      if( p->m_Q2 > m_Q2MaxOverall )
		m_Q2MaxOverall = p->m_Q2;
	    }
	  m_TargetsNeededOverall |= p->PredictionTargets();
	}

      IF_LOGv( MsgLevel )
      {
	DumpList( SampleSize );
	LOGv( MsgLevel, Name() << " contains " << m_list.size()
	      << " data points with statistics:" << endl
	      <<"    Item         Min         Max       Min^2       Max^2"
	      << fixed << setprecision(5) << endl
	      << "       W" << setw(12) << m_WMinOverall << setw(12) << m_WMaxOverall
	      << setw(12) << m_W2MinOverall << setw(12) << m_W2MaxOverall << endl
	      << "       Q" << setw(12) << m_QMinOverall << setw(12) << m_QMaxOverall
	      << setw(12) << m_Q2MinOverall << setw(12) << m_Q2MaxOverall << endl );
      }
    }
}

// Apply a cutoff so that everything above the cutoff is on our short-list

bool Cassandra::DataSet::ApplyCut(double W2Min, double Q2Min, bool bStrictCut,
				  size_t SampleSize, bool bShowCut, NodeSortFunc fSort )
{
  bool   bChanged = false;
  int    iCutQ = 0;
  int    iCutW = 0;
  int    iCutBoth = 0;
  const char * pszCut = Cassandra_StrictString( bStrictCut );

  // First, see whether there are any items on our short-list that should be cut
  for( list<DataNode *>::iterator i = m_l.begin() ; i != m_l.end() ; )
    {
      DataNode * p = * i;
      bool bCut = true;
      if( bStrictCut )
	{
	  if( p->m_cW2 <= W2Min )
	    {
	      if( p->m_Q2 <= Q2Min )
		iCutBoth++;
	      else
		iCutW++;
	    }
	  else if( p->m_Q2 <= Q2Min )
	    iCutQ++;
	  else
	    bCut = false;
	}
      else
	{
	  if( p->m_cW2 < W2Min )
	    {
	      if( p->m_Q2 < Q2Min )
		iCutBoth++;
	      else
		iCutW++;
	    }
	  else if( p->m_Q2 < Q2Min )
	    iCutQ++;
	  else
	    bCut = false;
	}
      if( bCut )
	{
	  // Cut this node
	  bChanged = true;
	  i = m_l.erase(i);
	  m_list.push_back(p);
	}
      else
	i++; // Nothing to do, already on the short list
    }

  // Now, see whether there are any cut items that should be short-listed
  for( list<DataNode *>::iterator i = m_list.begin() ; i != m_list.end() ; )
    {
      DataNode * p = * i;
      bool bCut = true;
      if( bStrictCut )
	{
	  if( p->m_cW2 <= W2Min )
	    {
	      if( p->m_Q2 <= Q2Min )
		iCutBoth++;
	      else
		iCutW++;
	    }
	  else if( p->m_Q2 <= Q2Min )
	    iCutQ++;
	  else
	    bCut = false;
	}
      else
	{
	  if( p->m_cW2 < W2Min )
	    {
	      if( p->m_Q2 < Q2Min )
		iCutBoth++;
	      else
		iCutW++;
	    }
	  else if( p->m_Q2 < Q2Min )
	    iCutQ++;
	  else
	    bCut = false;
	}
      if( !bCut )
	{
	  // Put this node in our short list
	  bChanged = true;
	  i = m_list.erase(i);
	  m_l.push_back(p);
	}
      else
	  i++; // Nothing to do, already cut
    }

  LOG( Mostly, Name() << fixed << setprecision(5) << " cut at W**2" << pszCut << W2Min
       << " and Q**2" << pszCut << Q2Min << " GeV**2 removes " << m_list.size()
       << " nodes of " << size() << " total." << endl
       << "W**2 cuts " << iCutW << " nodes, Q**2 cuts " << iCutQ
       << " nodes, " << iCutBoth << " nodes cut due to both. "
       << CutSize() << " nodes remain after cut." << endl );

  // If our list changed, we have some recalculations to do
  if( bChanged )
    {
      // Reset internal statistics
      ReleaseCutInternalOnly( false );
      // Show cut items if requested
      if( bShowCut )
	{
	  int iCount = 0;
	  DataNode * pLastNode = nullptr;
	  for( DataNode * p : m_list )
	    {
	      if( !pLastNode || typeid(*p) != typeid(*pLastNode) )
		{
		  if( pLastNode )
		    LOG( Always, Name() << " cut removes " << iCount << " "
			 << pLastNode->Name() << " data point(s)." << endl );
		  iCount = 0;
		  p->SerialiseHeader(LOG_STREAM, 5);
		  LOG( Always, endl );
		  pLastNode = p;
		}
	      p->Serialise(LOG_STREAM, 5);
	      LOG( Always, endl );
	      iCount++;
	    }
	  if( pLastNode )
	    LOG( Always, Name() << " cut removes " << iCount << " " << pLastNode->Name()
		 << " data point(s)." << endl );
	}
      // Calculate statistics if we have any data in our short-list
      if( CutSize() )
	{
	  m_v.resize( CutSize() );
	  m_l.sort( fSort );
	  // Save some statistics on the data points that survived the cut
	  int i = 0;
	  for( DataNode * p : m_l )
	    {
	      if( i == 0 )
		{
		  m_WMin = p->m_cW;
		  m_WMax = p->m_cW;
		  m_QMin = p->m_cQ;
		  m_QMax = p->m_cQ;
		}
	      else
		{
		  if( p->m_cW < m_WMin )
		    m_WMin = p->m_cW;
		  if( p->m_cW > m_WMax )
		    m_WMax = p->m_cW;
		  if( p->m_cQ < m_QMin )
		    m_QMin = p->m_cQ;
		  if( p->m_cQ > m_QMax )
		    m_QMax = p->m_cQ;
		}
	      m_TargetsNeeded |= p->PredictionTargets();
	      // Put data set into array, in source file order - to match covariance matrix
	      m_v[i++] = p;
	    }
	  // Say what's in the cut
	  LOG( Always, "    Item         Min         Max       Min^2       Max^2"
	       << fixed << setprecision(5) << endl
	       << "       W" << setw(12) << m_WMin << setw(12) << m_WMax
	       << setw(12) << m_WMin * m_WMin << setw(12) << m_WMax * m_WMax << endl
	       << "       Q" << setw(12) << m_QMin << setw(12) << m_QMax
	       << setw(12) << m_QMin * m_QMin << setw(12) << m_QMax * m_QMax << endl );

	  // Reserve space for the covariance matrix
	  m_MCovar.ResizeTo( CutSize(), CutSize() );
	  m_MCovarInv.ResizeTo( CutSize(), CutSize() );

	  // Sort list by energy so we can save a bit of execution time, then put data in vector
	  m_l.sort( DataNodeSortQXFile );

	  // Now show what's in the list
	  DumpList( SampleSize );
	}
    }
  return bChanged;
}

// Put everything in our short-list aside, where it won't form part of our main list
// but won't be deleted

void Cassandra::DataSet::ParkCut( size_t SampleSize )
{
  for( list<DataNode *>::iterator i = m_l.begin() ; i != m_l.end() ; )
    {
      DataNode * p = * i;
      i = m_l.erase(i);
      m_lParked.push_back(p);
    }
  ReleaseCutInternalOnly( false );
  ComputeOverallStats( Cassandra::LogLevel::Mostly, SampleSize );
}

void Cassandra::DataSet::ParkXBelow( double X )
{
  int iCutCount = 0;
  for( list<DataNode *>::iterator i = m_l.begin() ; i != m_l.end() ; )
    {
      DataNode * p = * i;
      if( p->m_x < X )
	{
	  i = m_l.erase(i);
	  m_lParked.push_back(p);
	  iCutCount++;
	}
      else
	i++;
    }
  for( list<DataNode *>::iterator i = m_list.begin() ; i != m_list.end() ; )
    {
      DataNode * p = * i;
      if( p->m_x < X )
	{
	  i = m_list.erase(i);
	  m_lParked.push_back(p);
	  iCutCount++;
	}
      else
	i++;
    }
  LOG( Often, iCutCount << " nodes parked because X<" << fixed << setprecision(2) << X << endl );
  ReleaseCutInternalOnly( false );
  ComputeOverallStats( Cassandra::LogLevel::Mostly, 0 );
}

// Make covariance matrix & invert it

void Cassandra::DataSet::ConstructCovar(void)
{
  // Construct the covariance matrix
  LOG( Often, "Constructing " << CutSize() << " by " << CutSize()
       << " covariance matrix for " << Name() << endl );
  for( unsigned int i = 0 ; i < CutSize() ; i++ )
    {
      for( unsigned int j = 0 ; j <= i ; j++ )
	{
	  bool bDiagonal = ( i == j ) ? true : false;
	  double covar = m_v[i]->Covariance( * m_v[j], bDiagonal );
	  m_MCovar[i][j] = covar;
	  m_MCovarInv[i][j] = covar;
	  if( !bDiagonal )
	    {
	      m_MCovar[j][i] = covar;
	      m_MCovarInv[j][i] = covar;
	    }
	}
    }

  LOG( Often, "Inverting " << CutSize() << " by " << CutSize() << " square matrix ... " << endl );
  LOG( Sometimes, scientific << m_MCovar << endl );
  Double_t det=-9.9e99;
  auto tStart = std::chrono::high_resolution_clock::now();
  m_MCovarInv.Invert(&det);
  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> tElapsed = tEnd - tStart;
  LOG( Often, "Matrix inversion took " << fixed << tElapsed.count() << " seconds." << endl );
  LOG( Sometimes, scientific << m_MCovarInv << endl << "determinant " << setprecision(10) << det << endl );
}

// Write the current prediction to the specified file - in case of restart

bool Cassandra::DataSet::WritePrediction(const char * pPDFSet, const char * pMassScheme,
					 const std::string sFinalFileName, int /*iSeq*/,
					 double W2Min, double Q2Min, bool bStrictCut,
					 int iLowestValidReplica) const
{
  bool bRet = false;
  if( !m_v.size() )
    {
      LOG( Always, "Bug: writing empty DataSet to output file " << sFinalFileName << endl );
      return bRet;
    }

  std::string sTempFileName( sFinalFileName );
  sTempFileName.push_back( '~' );
  std::ofstream fout(sTempFileName, ios_base::out | ios_base::trunc);
  if(!fout.is_open())
    {LOG( Always, "Unable to create output file " << sTempFileName << endl );}
  else
    {
      LOG( Rarely, "Writing prediction to " << sTempFileName << endl );

      // Header
      const char * pszCutName = Cassandra_StrictString( bStrictCut );
      while( * pszCutName == ' ' )
	pszCutName++;
      fout << "Cassandra Prediction " << Name() << " F2" << endl
	   << setw(21) << "PDFSet : " << pPDFSet << endl
	   << setw(21) << "MassScheme : " << pMassScheme << endl
	   << setw(21) << "PDFReplicas : " << fixed << setprecision(0) << m_iNumReplicas << endl
	   << setw(21) << "Count : " << CutSize() << endl
	   << scientific << setprecision(std::numeric_limits<double>::digits10)
	   << setw(21) << "W**2 : " << W2Min << endl
	   << setw(21) << "Q**2 : " << Q2Min << endl
	   << setw(21) << "Cut : " << pszCutName << endl;
      if( iLowestValidReplica != 0 )
	fout << setw(21) << "LowestValidReplica : " << setprecision(0) << iLowestValidReplica << endl;
      fout << endl;

      bool bFirst = true;
      for( DataNode * p : m_v )
	{
	  if( bFirst )
	    {
	      bFirst = false;
	      p->SerialiseHeader( fout, std::numeric_limits<int>::max(), iLowestValidReplica );
	      fout << endl;
	    }
	  p->Serialise( fout, std::numeric_limits<int>::max(), iLowestValidReplica );
	  fout << endl;
	}

      //Close the file
      fout.close();

      // Now delete the file we are replacing
      bool bFileExists;
      {
	std::ifstream is(sFinalFileName.c_str());
	if( is )
	  bFileExists = true;
	else
	  bFileExists = false;
      }
      int iError;
      if( bFileExists && ( iError = std::remove( sFinalFileName.c_str() ) ) != 0 )
	{ LOG( Always, "Error " << iError << " deleting " << sFinalFileName << endl ); }
      else if( ( iError = std::rename( sTempFileName.c_str(), sFinalFileName.c_str() ) ) != 0 )
	{ LOG( Always, "Error " << iError << " renaming " << sTempFileName
	       << "to " << sFinalFileName << endl ); }
      else
	{
	  bRet = true;
	  LOG( Sometimes, "Successfully written prediction to " << sFinalFileName << endl );
	}
    }
  return bRet;
}

bool Cassandra::BCDMSSet::MyFileType( string (&s)[7] )
{
  if( s[0] == "REACtion" && s[1] == "=" && s[2] == "muon"
      && s[4] == "-->" && s[5] == "muon" && s[6] == "X" )
    {
      if( s[3] == "P" )
	{
	  m_Target = Cassandra::DataNode::Target::Proton;
	  return true;
	}
      if( s[3] == "deuterium" )
	{
	  m_Target = Cassandra::DataNode::Target::Deuteron;
	  return true;
	}
    }
  return false;
}

void Cassandra::BCDMSSet::ParseBegin( void )
{
  m_bGoteBeam = false;
  m_eBeam = 0.0;
  //m_Target = Cassandra::DataNode::Target::None;
}

bool Cassandra::BCDMSSet::ParseLine(std::stringstream & ss, bool &bContinue)
{
  bool bPrint = false;
  switch(m_stage)
    {
      // Header
    case 0:
      {
	bPrint = true;
	string s1, s2, s3;
	ss >> s1 >> s2 >> s3;
	if( s1 != "(GeV**2)" || s2 != "R=0.0" || s3 != "stat." )
	  {
	    // See whether I can find the GeV this file applies to
	    string::size_type sPos = ss.str().find( " GeV" );
	    if( sPos != string::npos && sPos > 0 )
	      {
		string::size_type sPrevPos = ss.str().find_last_of( " \t", sPos - 1 );
		if( sPrevPos == string::npos )
		  sPrevPos = 0;
		  //{
		    //sPrevPos += 2;
		    string::size_type sLen = sPos - sPrevPos;
		    if( sLen > 0 )
		      {
			string::size_type numLen;
			double thisGeV = stod( ss.str().substr( sPrevPos, sLen ), &numLen );
			if( numLen != 0 )
			  {
			    // I've read a valid number
			    if( !m_bGoteBeam )
			      {
				m_bGoteBeam = true;
				LOG( Sometimes, "Beam energy " << thisGeV << " GeV." << endl );
			      }
			    else if( m_eBeam != thisGeV )
			      {
				LOG( Always, "Error: beam energy changed from " << m_eBeam
				     << " to " << thisGeV << " GeV." << endl );
				bContinue = false;
				bPrint = true;
			      }
			    else
			      LOG( Rarely, "Beam energy " << thisGeV << " GeV (again)." << endl );
			    m_eBeam = thisGeV;
			  }
		      }
		  //}
	      }
	  }
	else
	  {
	    // Found the beginning of the body
	    m_stage++;
	    if( !m_bGoteBeam )
	      {
		LOG( Always, "Error: beam energy unspecified." << endl );
		bContinue = false;
	      }
	    else
	      LOG( Sometimes, "Beam energy " << m_eBeam << " GeV." << endl );
	  }
      }
      break;

    case 1:  // Read data
      double d[8];
      d[7] = -9.9e99;
      ss >> d[0] >> d[1] >> d[2] >> d[3] >> d[4] >> d[5] >> d[6] >> d[7];
      if( d[7] == -9.9e99)
	{
	  LOG( Always, "Error reading record " << m_count + 1 << endl );
	  bPrint = true;
	  bContinue = false;
	}
      else
	{
	  // Add this data node to the end of my list
	  ++m_count;
	  DataNode * p;
	  if( m_Target == Cassandra::DataNode::Target::Proton )
	    {
	      p = new BCDMSNodeProton( m_eBeam, m_count, d[0], d[1], Datum(d[5], d[6], d[7]),
				       m_iNumReplicas/*, MASS_PROTON*/ );
	    }
	  else
	    {
	      p = new BCDMSNodeDeuteron( m_eBeam, m_count, d[0], d[1], Datum(d[5], d[6], d[7]),
					 m_iNumReplicas/*, MASS_DEUTERON*/ );
	    }
	  m_list.push_back(p);
	}
      break;
    }
  return bPrint;
}

bool Cassandra::BoNuSSet::MyFileType( string (&s)[7] )
{
  if( s[0] == "BoNuS" && s[1] == "F2N/F2D" && s[2] == "Measured" && s[3] == "data"
      && s[4] == "---" && s[5] == "submitted" && s[6] == "by" )
    return true;
  return false;
}

bool Cassandra::BoNuSSet::ParseLine(std::stringstream & ss, bool &bContinue)
{
  bool bPrint = false;
  switch(m_stage)
    {
      // Header
    case 0:  // Skip past the header
      {
	bPrint = true;
	string s1, s2, s3;
	ss >> s1 >> s2 >> s3;
	if( s1 == "x" && s2 == "Q^2" && s3 == "f2N/f2D" )
	  {
	    // Found the beginning of the body
	    m_stage++;
	  }
      }
      break;
    case 1:  // Read data
      {
	double   n1, n2, n3, n4, n5, n6, n7, dummy = -9.9e99;
	ss >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> dummy;
	/*ss >> n.m_x >> n.m_q2 >> n.m_ratio >> i
	  >> n.m_errStat >> n.m_errSys >> n.m_eBeam >> dummy;*/
	if(dummy!=1.0)
	  {
	    LOG( Always, "Error: dummy value != 1.0 reading record " << m_count + 1 );
	    bPrint = true;
	    bContinue = false;
	  }
	else
	  {
	    ++m_count;
	    Datum d(n3,n5,n6);
	    m_list.push_back(new BoNuSNode(n7, (int) n4, n1, n2, d));
	  }
      }
      break;
    }
  return bPrint;
}

bool Cassandra::JLABSet::MyFileType( string (&s)[7] )
{
  if( s[0] == "JLab" && s[1] == "E00106" && s[2] == "F2" )
    {
      if( s[3] == "proton" )
	{
	  m_Target = Cassandra::DataNode::Target::Proton;
	  return true;
	}
      if( s[3] == "deuterium" )
	{
	  m_Target = Cassandra::DataNode::Target::Deuteron;
	  return true;
	}
    }
  return false;
}

bool Cassandra::JLABSet::ParseLine(std::stringstream & ss, bool &bContinue)
{
  bool bPrint = false;
  switch(m_stage)
    {
      // Header
    case 0:  // Skip past the header
      {
	bPrint = true;
	string s1, s2, s3;
	ss >> s1 >> s2 >> s3;
	if( s1 == "x" && s2 == "Q^2" && s3[0] == 'F' && s3[1] == '2' )
	  {
	    // Found the beginning of the body
	    m_stage++;
	  }
      }
      break;
    case 1:  // Read data
      {
	double   n1, n2, n3, n4, n5, n6, n7 = -9.9e99, dummy = -9.9e99;
	ss >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7 >> dummy;
	/*ss >> n.m_x >> n.m_q2 >> n.m_ratio >> i
	  >> n.m_errStat >> n.m_errSys >> n.m_eBeam >> dummy;*/
	if(n7 != 1.0 || dummy != 1.0)
	  {
	    LOG( Always, "Error: dummy value != 1.0 reading record " << m_count + 1 );
	    bPrint = true;
	    bContinue = false;
	  }
	else
	  {
	    ++m_count;
	    DataNode * p;
	    if( m_Target == Cassandra::DataNode::Target::Proton )
	      p = new JLABNodeProton( 0, (int) n4, n1, n2, Datum(n3, n5, n6),
				       m_iNumReplicas/*, MASS_PROTON*/ );
	    else
	      p = new JLABNodeDeuteron( 0, (int) n4, n1, n2, Datum(n3, n5, n6),
				       m_iNumReplicas/*, MASS_DEUTERON*/ );
	    /*p->m_x = n1;
	    p->m_Q2 = n2;
	    p->m_expt = n3;
	    p->m_i = (int) n4;
	    p->m_errStat = n5;
	    p->m_errSys = n6;
	    // Save calculated variables
	    p->m_cQ = sqrt(p->m_Q2);*/
	    // Add this data node to the end of my list
	    m_list.push_back(p);
	  }
      }
      break;
    }
  return bPrint;
}

bool Cassandra::SLACSet::MyFileType( string (&s)[7] )
{
  bool bValid = false;
  if( s[0] == "SLAC" && s[1] == "F2" )
    {
      if( s[2] == "proton" )
	{
	  m_Target = Cassandra::DataNode::Target::Proton;
	  bValid = true;
	}
      if( s[2] == "deuteron" )
	{
	  m_Target = Cassandra::DataNode::Target::Deuteron;
	  bValid = true;
	}
      if( bValid )
	{
	  if( s[3] == "(reanalyzed)" )
	    {
	      m_bNNPDF = false;
	      LOG( Often, "SLAC data file provenance: (reanalyzed)" << endl );
	    }
	  else if( s[3] == "(NNPDF)" )
	    {
	      m_bNNPDF = true;
	      LOG( Often, "SLAC data file provenance: (NNPDF)" << endl );
	    }
	  else
	    {
	      //bValid = false;
	      // Actually, might as well assume it's SLAC format
	      // this will be checked as we read the header
	      m_bNNPDF = false;
	      LOG( Often, "SLAC data file provenance: \"" << s[3] << "\" (unknown)" << endl );
	    }
	}
    }
  return bValid;
}

bool Cassandra::SLACSet::ParseLine(std::stringstream & ss, bool &bContinue)
{
  bool bPrint = false;
  switch(m_stage)
    {
      // Header
    case 0:  // Skip past the header
      {
	bPrint = true;
	string s1, s2, s3;
	ss >> s1 >> s2 >> s3;
	if( s1 == "x" && s2 == "Q2" && s3 == "F2" )
	  {
	    // Found the beginning of the body
	    ss >> s1 >> s2;
	    if( !( ( m_bNNPDF && s1 == "stat" && s2 == "sys" ) || ( !m_bNNPDF && s1 == "i" && s2 == "stat" ) ) )
	      {
		// This is not a type we were expecting
		LOG( Always, "SLAC [NNPDF=" << m_bNNPDF << "] header unrecognised" << endl );
		bPrint = true;
		bContinue = false;
	      }
	    else
	      m_stage++;
	  }
      }
      break;
    case 1:  // Read data
      {
	double   n1, n2, n3, n4, n5, n6, n7 = -9.9e99;
	if( m_bNNPDF )
	  {
	    ss >> n1 >> n2 >> n3 >> n5 >> n6;
	    n4 = m_count + 1;
	    n7 = 1.;
	  }
	else
	  ss >> n1 >> n2 >> n3 >> n4 >> n5 >> n6 >> n7;
	if( ss.fail() || n7 != 1.0 )
	  {
	    LOG( Always, "Error reading record " << m_count + 1 << endl );
	    bPrint = true;
	    bContinue = false;
	  }
	else
	  {
	    ++m_count;
	    DataNode * p;
	    if( m_Target == Cassandra::DataNode::Target::Proton )
	      p = new SLACNodeProton( 0, (int) n4, n1, n2, Datum(n3, n5, n6),
				       m_iNumReplicas/*, MASS_PROTON*/ );
	    else
	      p = new SLACNodeDeuteron( 0, (int) n4, n1, n2, Datum(n3, n5, n6),
				       m_iNumReplicas/*, MASS_DEUTERON*/ );
	    /*p->m_x = n1;
	    p->m_Q2 = n2;
	    p->m_expt = n3;
	    p->m_i = (int) n4;
	    p->m_errStat = n5;
	    p->m_errSys = n6;
	    // Save calculated variables
	    p->m_cQ = sqrt(p->m_Q2);*/
	    // Add this data node to the end of my list
	    m_list.push_back(p);
	  }
      }
      break;
    }
  return bPrint;
}
