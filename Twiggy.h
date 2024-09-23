#if !defined(__CassMain_h__)
#define __CassMain_h__

# include "Cassandra.h"

namespace Cassandra
{
  // This is an entire set of data - doesn't need to understand details

  class Bin
  {
  public:
    double central; // i.e. central value
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ey;
  public:
    Bin( double central_value = 0 ) { central = central_value; }
  public:
    inline void Add( double x_coord, double y_coord, double error_y )
    {
      x.push_back( x_coord );
      y.push_back( y_coord );
      ey.push_back( error_y );
    }
  };

  class Model
  {
  protected:
    const int  m_iModelNum;
    const bool m_bAdditive;

  public:
    std::vector<HT> m_Model;

  public:
  Model( int iModelNum = 0, bool bAdditive = false )
    : m_iModelNum{ iModelNum }, m_bAdditive{ bAdditive } {}

    //Services
  public:
    inline HT & operator[]( int index ) { return m_Model[index]; }
    inline std::size_t size( void ) const { return m_Model.size(); }
    double RawModel( const DataNode * p, int k ) const;
    double HigherTwist( const DataNode * p, int k, int iReplica = 0 ) const;
    inline double NewTheory( const DataNode * p, int k, int iReplica = 0 ) const
        { return p->m_Theory[iReplica] + HigherTwist( p, k, iReplica ); }
  };

  struct Parton
  {
    int          m_APFELIndex;
    const char * m_szName;    // File / programmatic short name
    const char * m_szChartName; // Full-text
  };

  struct CutChi
  {
    static constexpr int NoSeq{ std::numeric_limits<int>::min() };
    int    iSeq;
    double W2;
    double Q2;
    bool   bStrict;
    double Chi;
    double ChiTwist;
    double ChiTwistReweight;
    CutChi( int Seq, double w2, double q2, bool Strict ) : iSeq{Seq},W2{w2},Q2{q2},bStrict{Strict} {}
    std::string SeqString( const char *pszPrefix ) const
    {
      std::string s;
      if( pszPrefix )
	s = pszPrefix;
      if( iSeq != NoSeq )
      {
	if( !s.empty() )
	  s.append( 1, '_' );
	s.append( std::to_string( iSeq ) );
      }
      return s;
    }
  };

  class DataManager
  {
  protected:
    static const char * m_pszChiSuffix;
    static const char m_szModelSuffix[];
    static const char m_szModelHeader[];
    static const char m_szPDFSuffix[];
    static const char m_szPDFHeader[];
    static const char m_szOptionStdDev[];

  public:
    // Unknown must be last
    enum NodeType : int {
      BCDMS, SLAC, JLAB, BoNuS, Unknown };

  protected:
    Model m_ht;	// Which model to run
    DataSet     * m_DataSets[Unknown];
    ModelParams & m_Params;
    NodeSortFunc const m_TwiggySort;
    const int m_NumOutliers;
    const char * const m_pszModelPrefix;
    const double m_XCutoff;
    const char * const m_pszPDFSet;
    const double m_dQChartScale;
    const double m_xBinSize;
    const double m_Q2BinSize;

    const char * m_pszTarget;

    // Stats for the whole dataset ... don't change unless we park
    int         m_NDat; // Total number of data points
    double      m_WMinOverall, m_WMaxOverall;
    double      m_QMinOverall, m_QMaxOverall;
    double      m_W2MinOverall, m_W2MaxOverall;
    double      m_Q2MinOverall, m_Q2MaxOverall;

  public:
    DataManager(ModelParams & mp, int iNumReplicas, NodeSortFunc TwiggySort,
		int NumOutliers, const char * pszModelPrefix, double XCutoff,
		int iModelNum, bool bAdditive, const char * pszDefaultPDFSet,
		double dQChartScale, double xBinSize, double Q2BinSize );
    ~DataManager();
    DataSet & operator[]( NodeType nt ) { return * m_DataSets[nt]; }

  protected:
    double SeriesValue( const DataNode &n, int iSeries ) const;
    double SeriesValueBinned( const DataNode &n, int iSeries ) const;
    Int_t ValidateMatrices( const TMatrixDSym &MData,
			    const TMatrixDSym &MTheory,
			    const std::vector<DataNode *> & vNodes ) const;

  public:
    size_t size(void) const noexcept;
    void PlotHT( const DataSet &ds,
		 const std::string &sFileNameSeq, const std::string &sCut,
		 double * HTMean, const TMatrixDSym & MTheoryCovar ) const;
    void CreateGraphs( const TMatrixDSym & MData, const TMatrixDSym & MTheory,
		       const std::vector<DataNode *> & vNodes, double dNeff,
		       const std::string &sNamePrefix, 
		       const std::string &sFileNameSeq,
		       const std::string &sCut,
		       const std::string &sCounts,
		       double * HTMean ) const;
    void PlotCorrelation( const TMatrixDSym & M, double dNeff,
			  const std::string & sChartTitle,
			  const std::string &sFileName,
			  const std::string &sCut,
			  const std::string &sCounts ) const;
    void PlotDiagonals( const TMatrixDSym & MData, const TMatrixDSym & MTheory,
			const std::vector<DataNode *> & vNodes, double dNeff,
			const std::string &sNamePrefix, 
			const std::string &sFileNameSeq,
			const std::string &sCut,
			const std::string &sCounts ) const;
    bool LoadHTModel( std::vector<HT> &ht, int iModelNum ) const;
    // Plot the data set and save useful info
    void PlotDataSet( DataSet &ds, const char *pszOutFilePrefix, const CutChi &cc,
		      int *piNumDataSets = nullptr, std::string *psCounts = nullptr,
		      std::string *psCuts = nullptr );
    bool MakeOnePrediction( DataSet & ds, const char * pszOutFilePrefix,
			    CutChi & cc, size_t SampleSize, bool bShowCut );
    bool LoadData( const char * fileName );
    void ApplyCutoff( double W2, double Q2, bool bStrictCut,
		      size_t SampleSize, bool bShowCut );
    void ParkXBelow( double X );
    void ComputeOverallStats( size_t SampleSize );
    bool InitPrediction( const char * pszOutFilePrefix, int iNumReplicas );
    bool MakePrediction( const char * pszOutFilePrefix,
			 bool bStrictFirstCut,
			 int SampleSize, bool bDumpOnly, bool bShowCut,
			 double W2Min, double W2Max,
			 double Q2Min, double Q2Max, int NumSteps );
  };
}
#endif
