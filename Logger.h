#if !defined(__Logger_h__)
#define __Logger_h__

# include <cassert>
# include <fstream>
# include <glob.h>
# include <iomanip>
# include <iostream>
# include <ostream>
# include <sstream>

namespace Cassandra
{
  // Simple logger
  class Logger : public std::ostream
  {
  public:
    enum class LogLevel { Always, Mostly, Often, Sometimes, Rarely };
    LogLevel m_LogLevel; // Only show messages more important than this
    //static const char endl;

  protected:
    LogLevel m_MsgLevel; // Message level used when logging
    std::ofstream * o;

  public:
    Logger(LogLevel l) : m_LogLevel{l}, o{nullptr} {}
  Logger(void) : Logger(LogLevel::Always) {std::cout << "Logger()" << std::endl;}
    Logger(LogLevel l, const char * pFileName) : m_LogLevel{l},
      o{new std::ofstream( pFileName, std::ios::out | std::ios::ate ) } {}
    Logger(const char * pFileName) : Logger(LogLevel::Always, pFileName) {}
    ~Logger() { delete o; std::cout << "~Logger()" << std::endl;}

  public:
    LogLevel operator()(void)
      { return m_LogLevel; }
    Logger& operator()(LogLevel MsgLevel)
      { m_MsgLevel = MsgLevel; return *this; }
    template <typename T>
    inline Logger& operator<<(T const & val)
      {
	if( m_MsgLevel <= m_LogLevel )
	  {
	    if( o ) (*o) << val;
	    else std::cout << val;
	  }
	return * this;
      }
    void LogToFile(const char * pFileName)
      {
	delete o;
	o = new std::ofstream( pFileName, std::ios::out | std::ios::ate );
	}
  };

  //const char Logger::endl = '\n';
}
#endif
