#if !defined(__Cassandra_h__)
#define __Cassandra_h__

# include <cassert>
# include <cfenv>
# include <chrono>
# include <cmath>
# include <cstdio>
# include <fstream>
# include <functional>
# include <glob.h>
# include <iomanip>
# include <iostream>
# include <list>
# include <memory>
# include <ostream>
# include <sstream>
# include <typeinfo>
# include <vector>

// CERN root libraries
#include "TBox.h"
#include "TCanvas.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TRootCanvas.h"
#include "TVectorD.h"

# include "APFEL/APFEL.h"

#define MASS_PROTON   0.93827
#define MASS_DEUTERON 1.87705
#define Q2_DEFAULT_CUT 0.99999999999999

#define IF_LOG(level) if(Cassandra::LogLevel::level <= Global_Log_Level)
#define IF_LOGv(level) if(level <= Global_Log_Level)
#define LOG_STREAM ( Cassandra::Global_Log_File \
		     ? ( ( std::ostream& ) * Cassandra::Global_Log_File ) \
		     : cout )
#define LOG(level,message) { if(Cassandra::LogLevel::level <= Global_Log_Level) \
			       { LOG_STREAM << message; } }
#define LOGv(level,message) { if(level <= Global_Log_Level) \
			       { LOG_STREAM << message; } }

namespace Cassandra
{
  // Log file
  enum class LogLevel { Always, Mostly, Often, Sometimes, Rarely };
  extern LogLevel Global_Log_Level; // Show messages more important than this
  extern std::unique_ptr<std::ofstream> Global_Log_File;

  // Common strings
  extern const char * pszProton;
  extern const char * pszDeuteron;
  extern const char * pszBCDMS;
  extern const char * pszSLAC;
  extern const char * pszJLAB;
  extern const char * pszBoNuS;

  class Matrix2D
  {
  protected:
    double * m_pData;
    int m_i1;
    int m_i2;

  public:
    void Release()
        { if(m_pData) delete [] m_pData; m_pData=nullptr; m_i1=0; m_i2=0; }
    Matrix2D() : m_i1{0}, m_i2{0} { m_pData = nullptr; }
    Matrix2D(int i2, int i1) : m_i1{i1}, m_i2{i2} { m_pData=new double[i1*i2]; }
    ~Matrix2D() { Release(); }
    void Resize(int i2, int i1)
        { Release(); m_i1=i1; m_i2=i2; m_pData=new double[i1*i2]; }

  public:
    double * operator[](int i2) { return m_pData + i2 * m_i1; }
  };

  class HT
  {
  public:
    double p[4];
    double weight;
    double chisq;
  };

  class StartStep
  {
  public:
    double Start;
    double Step;
    int    Num;
  public:
    StartStep( double dStart = 1., double dStep = 1., int iNum = 10 )
      : Start{ dStart }, Step{ dStep }, Num{ iNum } {}

  public:
    void FromTo( double dFrom, double dTo, int iStep = 10 );
  };

  class ModelParams
  {
  protected:
    std::vector<StartStep> m_vSS;

  public:
    ModelParams(int iNumParams);
    ModelParams();
    StartStep & operator[](int iIndex) { return m_vSS[iIndex]; }
  };

  class Datum
  {
    //friend std::ostream& operator<<(std::ostream& os, const Datum& d);
  public:
    double m_Value;   // Experimental measurement
    double m_errStat; // Statistical error
    double m_errSys;  // Systematic / uncorrelated error
  public:
  Datum() : m_Value{0}, m_errStat{0}, m_errSys{0} {};
  Datum(double Value, double errStat, double errSys)
    : m_Value{Value}, m_errStat{errStat}, m_errSys{errSys} {};
  Datum(const Datum &r)
    : m_Value{r.m_Value}, m_errStat{r.m_errStat}, m_errSys{r.m_errSys} {};
  public:
    bool operator==(Datum &r) const { return m_Value == r.m_Value
		  && m_errStat == r.m_errStat && m_errSys == r.m_errSys;}
    bool operator!=(Datum &r) const { return m_Value != r.m_Value
		  || m_errStat != r.m_errStat || m_errSys != r.m_errSys;}

    void Set( double Value, double errStat, double errSys )
    { m_Value = Value; m_errStat = errStat; m_errSys = errSys; }
    double ErrorSquared(void) const
    { return m_errStat * m_errStat + m_errSys * m_errSys; }
    double LowerLimit(void) const
    { return m_Value - (m_errStat + m_errSys); }
    double UpperLimit(void) const
    { return m_Value + (m_errStat + m_errSys); }
    double RelativeError(double Value)
    {
      if( Value < LowerLimit() )
	return Value - LowerLimit();
      if( Value > UpperLimit() )
	return Value - UpperLimit();
      return 0.0;
    }
  };

  // Base class for data elements read from file
  // (overload to use)

  class DataNode
  {
    //friend std::ostream& operator<<(std::ostream& os, const DataNode& dn);
    friend class NodeParser;

  public:
    enum Target { None = 0, Proton = 1, Neutron = 2, Deuteron = 4 };
 
    // Data from  input file
    const double m_eBeam;
    const int    m_i;
    const double m_x;
    const double m_Q2;

    // Calculated values
    const double m_cQ, m_cy;
    const double m_cW2, m_cW;

    // Each constructed with one member, but derived classes can overload
    std::vector<Datum>  m_Data;     // Data points with errors from file
    std::vector<double> m_Theory;   // Theoretical prediction

  public:
    virtual ~DataNode() {}
    // Default constructor
    DataNode() : m_eBeam{0.0}, m_i{0}, m_x{0.0}, m_Q2{0.0},
		 m_cQ{0.0}, m_cy{0.0}, m_cW2{0.0}, m_cW{0.0},
		 m_Data(1,Datum()), m_Theory(1,0) {}
    // Initialising constructor
    DataNode(double eBeam, int i, double x, double Q2,
	     size_t NumData, const Datum &Data, size_t NumTheory,
	     double WMassTerm = 0.)
               : m_eBeam{eBeam}, m_i{i}, m_x{x}, m_Q2{Q2},
		 m_cQ{sqrt(m_Q2)},
		 m_cy{ ( m_eBeam == 0. ) ? 0. : m_cQ / m_eBeam },
		 m_cW2{m_Q2*(1-m_x)/m_x + WMassTerm * WMassTerm},
		 m_cW{sqrt(m_cW2)},
		 m_Data(NumData,Data), m_Theory(NumTheory,0) {}
    // Copy constructor
      DataNode( const DataNode &r) : m_eBeam{r.m_eBeam}, m_i{r.m_i},
                  m_x{r.m_x}, m_Q2{r.m_Q2},
			      m_cQ{r.m_cQ}, m_cy{r.m_cy},
			      m_cW2{r.m_cW2}, m_cW{r.m_cW},
			      m_Data(r.m_Data), m_Theory(r.m_Theory) {}
    // Move constructor
    DataNode( DataNode &&r) : m_eBeam{r.m_eBeam}, m_i{r.m_i},
                  m_x{r.m_x}, m_Q2{r.m_Q2},
                  m_cQ{r.m_cQ}, m_cy{r.m_cy},
			      m_cW2{r.m_cW2}, m_cW{r.m_cW},
			      m_Data( std::move( r.m_Data ) ),
			      m_Theory( std::move( r.m_Theory ) ) {}
  public:
    // Helper functions - only really for intended for DataNode
  protected:
    enum class FieldType { Data, Theory };
    void SerialiseFieldName(FieldType Type, int iIndex,
			    std::ostream& os, int iLen ) const;

  public:
    // Derived classes may choose to overide
    virtual double ScalingFactor( int /*iDataIndex*/ = 0 ) const { return 1.; }

    // Derived classes can overload, but MUST call base first
    virtual void SerialiseHeader(std::ostream& os, int iMaxTheory,
				 int iMinTheory = 0) const;
    virtual void Serialise(std::ostream& os, int iMaxTheory,
				 int iMinTheory = 0) const;

    // Derived classes must overload
    virtual const char * Name(void) const noexcept = 0;
    virtual Target PredictionTargets( void ) const = 0;
    virtual const std::string& FieldName(FieldType iType, int iIndex,
				  bool &bAddSuffix, int &iSuffix ) const = 0;
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const = 0;
    virtual void MakePrediction(Target t, int iReplica) = 0;
    // Identifies type: SLAC; JLAB; BCDMS; etc
    // NB: Pointer equality guarantees string equality
    virtual const char * Type(void) const noexcept = 0;

    // Services provided to consumers
    //void SetW2( double W2 ) { m_cW2 = W2; m_cW = sqrt( W2 ); }
    double ScaledValue( int iDataIndex = 0 ) const
             { return m_Data[iDataIndex].m_Value * ScalingFactor(iDataIndex); }
    double ScaledErrSys( int iDataIndex = 0 ) const
             { return m_Data[iDataIndex].m_errSys * ScalingFactor(iDataIndex); }
    double ScaledErrStat( int iDataIndex = 0 ) const
             { return m_Data[iDataIndex].m_errStat * ScalingFactor(iDataIndex); }
  };

  class DataNodeSingle : public DataNode
  {
  public:
    //virtual ~DataNodeSingle() {}
    DataNodeSingle(double eBeam, int i, double x, double Q2, const Datum &Data,
		   int iNumReplicas, double WMassTerm = 0.)
      : DataNode(eBeam, i, x, Q2, 1, Data, iNumReplicas, WMassTerm) {}
    virtual const std::string& FieldName(FieldType iType, int iIndex,
					 bool &bAddSuffix, int &iSuffix ) const;
    // Name of the target: proton; deuteron; etc.
    // NB: Pointer equality guarantees string equality
    virtual const char * Name(void) const noexcept = 0;
    virtual Target PredictionTargets( void ) const = 0;
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const = 0;
    virtual void MakePrediction(Target t, int iReplica) = 0;
    // Identifies type: SLAC; JLAB; BCDMS; etc
    // NB: Pointer equality guarantees string equality
    virtual const char * Type(void) const noexcept = 0;
  };

  class NodeProton : public DataNodeSingle
  {
  public:
    //virtual ~NodeProton() {}
    NodeProton(double eBeam, int i, double x, double Q2, const Datum &Data,
	       int iNumReplicas, double WMassTerm = 0.)
      : DataNodeSingle(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual const char * Name(void) const noexcept { return pszProton; }
    virtual Target PredictionTargets( void ) const { return Target::Proton; }
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const = 0;
    virtual void MakePrediction(Target t, int iReplica);
    // Identifies type: SLAC; JLAB; BCDMS; etc
    // NB: Pointer equality guarantees string equality
    virtual const char * Type(void) const noexcept = 0;
  };

  class NodeDeuteron : public DataNodeSingle
  {
  public:
    //virtual ~NodeDeuteron() {}
    NodeDeuteron(double eBeam, int i, double x, double Q2,const Datum &Data,
		 int iNumReplicas, double WMassTerm = 0.)
      : DataNodeSingle(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual const char * Name(void) const noexcept { return pszDeuteron; }
    virtual Target PredictionTargets( void ) const { return Target::Deuteron; }
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const = 0;
    virtual void MakePrediction(Target t, int iReplica);
    // Identifies type: SLAC; JLAB; BCDMS; etc
    // NB: Pointer equality guarantees string equality
    virtual const char * Type(void) const noexcept = 0;
  };

  class BCDMSNodeProton : public NodeProton
  {
  public:
    //virtual ~BCDMSNodeProton() {}
    BCDMSNodeProton(double eBeam, int i, double x, double Q2, const Datum &Data,
		    int iNumReplicas, double WMassTerm = 0.)
    : NodeProton(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const;
    virtual const char * Type(void) const noexcept { return pszBCDMS; }
  };

  class BCDMSNodeDeuteron : public NodeDeuteron
  {
  public:
    //virtual ~BCDMSNodeDeuteron() {}
    BCDMSNodeDeuteron(double eBeam, int i, double x, double Q2,const Datum &Data,
		      int iNumReplicas, double WMassTerm = 0.)
      : NodeDeuteron(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const;
    virtual const char * Type(void) const noexcept { return pszBCDMS; }
  };

  class JLABNodeProton : public NodeProton
  {
  public:
    //virtual ~JLABNodeProton() {}
    JLABNodeProton(double eBeam, int i, double x, double Q2, const Datum &Data,
		   int iNumReplicas, double WMassTerm = 0.)
    : NodeProton(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const;
    virtual const char * Type(void) const noexcept { return pszJLAB; }
  };

  class JLABNodeDeuteron : public NodeDeuteron
  {
  public:
    //virtual ~JLABNodeDeuteron() {}
    JLABNodeDeuteron(double eBeam, int i, double x, double Q2,const Datum &Data,
		     int iNumReplicas, double WMassTerm = 0.)
      : NodeDeuteron(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const;
    virtual const char * Type(void) const noexcept { return pszJLAB; }
    // The JLAB Deuteron data are *2, so scale them down (see data file header)
    virtual double ScalingFactor( int /*iDataIndex*/ = 0 ) const { return 0.5; }
  };

  class SLACNodeProton : public NodeProton
  {
  public:
    //virtual ~SLACNodeProton() {}
    SLACNodeProton(double eBeam, int i, double x, double Q2, const Datum &Data,
		   int iNumReplicas, double WMassTerm = 0.)
    : NodeProton(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const;
    virtual const char * Type(void) const noexcept { return pszSLAC; }
  };

  class SLACNodeDeuteron : public NodeDeuteron
  {
  public:
    //virtual ~SLACNodeDeuteron() {}
    SLACNodeDeuteron(double eBeam, int i, double x, double Q2,const Datum &Data,
		     int iNumReplicas, double WMassTerm = 0.)
      : NodeDeuteron(eBeam, i, x, Q2, Data, iNumReplicas, WMassTerm) {}
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const;
    virtual const char * Type(void) const noexcept { return pszSLAC; }
  };

  class BoNuSNode : public DataNode
  {
  public:
    //double m_ratio, m_errStat, m_errSys;
    //double m_cF2n, m_cF2d, m_cRatio;

  public:
    //virtual ~BoNuSNode() {}
    BoNuSNode(double eBeam, int i, double x, double Q2, const Datum &Data)
      : DataNode(eBeam, i, x, Q2, 1, Data, 3) {}
    virtual const char * Name(void) const noexcept { return "BoNuS node"; }
    /*BoNuSNode(const BoNuSNode &r);
    virtual void Serialise( std::ostream & os ) const;
    virtual void SerialiseHeader( std::ostream & os ) const;*/
    virtual Target PredictionTargets( void ) const
    {return static_cast<Target>( (int) Target::Proton | (int)Target::Deuteron);}
    virtual const std::string& FieldName(FieldType iType, int iIndex,
				  bool &bAddSuffix, int &iSuffix ) const;
    virtual void MakePrediction(Target t, int iReplica);
    virtual double Covariance(const DataNode &nr, bool bDiagonal) const;
    virtual const char * Type(void) const noexcept { return pszBoNuS; }
  };

  // Generic parser
  // (overload to use)

  // Helper functions to sort lists of DataNode
  typedef bool ( * NodeSortFunc ) ( const DataNode * p1, const DataNode * p2 );
  // File order
  bool DataNodeSortFile( const DataNode * p1, const DataNode * p2 );
  // Sort by Q and X first, then file order
  bool DataNodeSortQXFile( const DataNode * p1, const DataNode * p2 );
  // Sort each file by W, then file order
  bool DataNodeSortW( const DataNode * p1, const DataNode * p2 );
  // Sort by Q2 only
  bool DataNodeSortQ2Only( const DataNode * p1, const DataNode * p2 );
  // Sort by X only
  bool DataNodeSortXOnly( const DataNode * p1, const DataNode * p2 );
  // Sort by DataSet type, then Q2, then X
  bool DataNodeSortTypeQ2X( const DataNode * p1, const DataNode * p2 );
  // Sort by DataSet type, then X, then Q2
  bool DataNodeSortTypeXQ2( const DataNode * p1, const DataNode * p2 );

  class DataSet
  {
    friend class DataManager;
  protected:
    ModelParams & m_Params;
    int           m_iNumReplicas;
    std::string   m_sName;

    std::list<DataNode *>   m_list; // Cut list (owns DataNode members)
    std::list<DataNode *>   m_lParked;  // Nodes we're not using just now

    // Following variables are the short list - i.e. data surviving cut
    std::list<DataNode *>   m_l;  // Short list (owns DataNode members)
    std::vector<DataNode *> m_v;  // Never  owns DataNode members
    TMatrixDSym m_MCovar, m_MCovarInv;
    bool        m_bCovarBuilt;
    // Stats for the whole dataset ... don't change unless we park
    double      m_WMinOverall, m_WMaxOverall;
    double      m_QMinOverall, m_QMaxOverall;
    double      m_W2MinOverall, m_W2MaxOverall;
    double      m_Q2MinOverall, m_Q2MaxOverall;
    int         m_TargetsNeededOverall;
    // Stats for the current cut
    double      m_WMin, m_WMax;
    double      m_QMin, m_QMax;
    int         m_TargetsNeeded;

    // Variables/functions used during parsing only
  protected:
    int m_stage;
    int m_count;

  private:
    void DumpListUncoolHelper(size_t sample, int &j,
			      int &LastNodeCount,
			      bool & bHeaderPrinted,
			      DataNode * & pLastNode,
			      DataNode * p);
    void ReleaseCutInternalOnly(bool bCutAll);
    void ComputeOverallStats(Cassandra::LogLevel MsgLevel, size_t SampleSize);

  public:
    DataSet( ModelParams & mp, int iNumReplicas, const char * pszName )
      : m_Params{ mp }, m_iNumReplicas{iNumReplicas}, m_sName{ pszName }
      { ReleaseCut(); }
    virtual ~DataSet() { EmptyList(); }

  protected:
    // Derived classes must provide implementations of these functions
    //virtual const char * Name(void) const noexcept = 0;
    virtual bool MyFileType( string (& Tokens)[7] ) = 0;
    virtual bool ParseLine(std::stringstream & ss, bool &bContinue) = 0;
    virtual DataNode * NewDataNode(DataNode::Target tTarget, double eBeam, int i,
			      double x, double Q2, const Datum &Data,
			      int iNumReplicas, double WMassTerm = 0.) const = 0;
    virtual DataSet * NewDataSet(void) const = 0;

    // Derived classes may choose to overide
    virtual void ParseBegin(void){} // Do nothing when parse new file

  public:
    // Services provided by the base class implementation
    const char * Name(void) const noexcept { return m_sName.c_str(); }
    size_t CutSize(void) const noexcept { return m_l.size(); }
    size_t size(void) const noexcept { return CutSize() + m_list.size(); }
    void   EmptyList(void);
    void   DumpList(size_t SampleSize);
    bool   ApplyCut(double W2Min, double Q2Min, bool bStrictCut,
		    size_t SampleSize, bool bShowCut,
		    NodeSortFunc fSort );
    void   ReleaseCut(void) { ReleaseCutInternalOnly(true); }
    void   ParkCut( size_t SampleSize );
    void   ParkXBelow( double X );
    void   ConstructCovar(void); // Make covariance matrix & invert it
    // Construct a new data node from a stream
    DataNode * NewDataNode(std::stringstream & ss,
			   int iLowestValidReplica,
			   int iNumReplicas) const;

    // Reconstipate from stream (file)
    bool Load(std::istream &is, int &iLowestValidReplica,
	      const char * pszPDFSet, const char * pszMassScheme,
	      int iNumReplicas, int iCount, double W2, double Q2, bool bStrict,
	      bool bValidateParameters = true);
    // Write all replicas >= iLowestValidReplica to disk
    bool WritePrediction(const char * pPDFSet, const char * pMassScheme,
			 const std::string sFinalFileName,
			 int iSeq, double W2Min, double Q2Min,
			 bool bStrictCut, int iLowestValidReplica) const;

  private:
    void NewFile(void);
  };

  // Parser for BCDMS data

  class BCDMSSet : public DataSet
  {
  protected:
    double m_eBeam;
    bool   m_bGoteBeam;
    DataNode::Target m_Target;

  protected:
    //virtual const char * Name(void) const noexcept { return pszBCDMS; }
    virtual bool MyFileType( string (& Tokens)[7] );
    virtual void ParseBegin(void);
    virtual bool ParseLine(std::stringstream & ss, bool &bContinue);
    virtual DataNode * NewDataNode(DataNode::Target tTarget, double eBeam, int i,
				   double x, double Q2, const Datum &Data,
				   int iNumReplicas, double WMassTerm = 0.) const
    { if( tTarget == DataNode::Target::Proton )
	return new BCDMSNodeProton(eBeam,i,x,Q2,Data,iNumReplicas,WMassTerm);
      return new BCDMSNodeDeuteron(eBeam,i,x,Q2,Data,iNumReplicas,WMassTerm);}
    virtual DataSet * NewDataSet(void) const
                      { return new BCDMSSet( m_Params, m_iNumReplicas ); }

  public:
  BCDMSSet( ModelParams & mp, int iNumReplicas, const char * szName = pszBCDMS )
    : DataSet( mp, iNumReplicas, szName ) {}
  };

  // Parser for BoNuS data

  class BoNuSSet : public DataSet
  {
  protected:
    //virtual const char * Name(void) const noexcept { return pszBoNuS; }
    virtual bool MyFileType( string (& Tokens)[7] );
    virtual bool ParseLine(std::stringstream & ss, bool &bContinue);
    virtual DataNode * NewDataNode(DataNode::Target /*tTarget*/, double eBeam, int i,
				   double x, double Q2, const Datum &Data,
				   int /*iNumReplicas*/, double /*WMassTerm*/ = 0.) const
    { return new BoNuSNode(eBeam,i,x,Q2,Data); }
    virtual DataSet * NewDataSet(void) const
                      { return new BoNuSSet( m_Params, m_iNumReplicas ); }

  public:
  BoNuSSet( ModelParams & mp, int iNumReplicas, const char * szName = pszBoNuS )
    : DataSet( mp, iNumReplicas, szName ) {}
  };

  class JLABSet : public DataSet
  {
  protected:
    DataNode::Target m_Target;

  protected:
    //virtual const char * Name(void) const noexcept { return pszJLAB; }
    virtual bool MyFileType( string (& Tokens)[7] );
    virtual bool ParseLine(std::stringstream & ss, bool &bContinue);
    virtual DataNode * NewDataNode(DataNode::Target tTarget, double eBeam, int i,
				   double x, double Q2, const Datum &Data,
				   int iNumReplicas, double WMassTerm = 0.) const
    { if( tTarget == DataNode::Target::Proton )
	return new JLABNodeProton(eBeam,i,x,Q2,Data,iNumReplicas,WMassTerm);
      return new JLABNodeDeuteron(eBeam,i,x,Q2,Data,iNumReplicas,WMassTerm);}
    virtual DataSet * NewDataSet(void) const
                      { return new JLABSet( m_Params, m_iNumReplicas ); }

  public:
  JLABSet( ModelParams & mp, int iNumReplicas, const char * szName = pszJLAB )
    : DataSet( mp, iNumReplicas, szName ) {}
  };

  class SLACSet : public DataSet
  {
  protected:
    DataNode::Target m_Target;
    bool m_bNNPDF; // There's a different version of data from NNPDF

  protected:
    //virtual const char * Name(void) const noexcept { return pszSLAC; }
    virtual bool MyFileType( string (& Tokens)[7] );
    virtual bool ParseLine(std::stringstream & ss, bool &bContinue);
    virtual DataNode * NewDataNode(DataNode::Target tTarget, double eBeam, int i,
				   double x, double Q2, const Datum &Data,
				   int iNumReplicas, double WMassTerm = 0.) const
    { if( tTarget == DataNode::Target::Proton )
	return new SLACNodeProton(eBeam,i,x,Q2,Data,iNumReplicas,WMassTerm);
      return new SLACNodeDeuteron(eBeam,i,x,Q2,Data,iNumReplicas,WMassTerm);}
    virtual DataSet * NewDataSet(void) const
    { return new SLACSet( m_Params, m_iNumReplicas ); }

  public:
  SLACSet( ModelParams & mp, int iNumReplicas, const char * szName = pszSLAC )
    : DataSet( mp, iNumReplicas, szName ) {}
  };

  // Serialise a TMatrixD to an output stream
  std::ostream& operator<<(std::ostream& os, const TMatrixD &m);
  // Serialise a Datum to an output stream
  std::ostream& operator<<(std::ostream& os, const Datum& d);
  // String describing strictness of cut
  const char * Cassandra_StrictString( bool bStrictCut );
}
#endif
