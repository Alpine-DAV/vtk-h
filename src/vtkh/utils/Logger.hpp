#ifndef VTK_H_LOGGER_HPP
#define VTK_H_LOGGER_HPP

#include <vtkh/filters/util.hpp>
#include <stack>

//from rover logging
namespace vtkh
{

class Logger
{
public:
  static Logger *GetInstance(const std::string& name);

  ~Logger();
  void write(const int level, const std::string &message, const char *file, int line);
  std::ofstream & GetStream() { return Stream; }

protected:
  Logger(const std::string& name);
  Logger(Logger const &);
  std::ofstream Stream;

  static std::map<std::string, Logger*> Loggers;
};

class DataLogger
{
public:
  ~DataLogger();
  static DataLogger *GetInstance();
  void OpenLogEntry(const std::string &entryName);
  void CloseLogEntry(const double &entryTime);

  template<typename T>
  void AddLogData(const std::string key, const T &value)
  {
    this->Stream<<key<<" "<<value<<"\n";
  }

  std::stringstream& GetStream() { return Stream; }
  void WriteLog();
protected:
  DataLogger();
  DataLogger(DataLogger const &);
  std::stringstream Stream;
  static class DataLogger* Instance;
  std::stack<std::string> Entries;
};

#ifdef TRACE_DEBUG
#define DBG(msg) vtkh::Logger::GetInstance("out")->GetStream()<<msg
#define WDBG(msg) vtkh::Logger::GetInstance("wout")->GetStream()<<msg
#else
#define DBG(msg)
#define WDBG(msg)
#endif

#ifdef ROVER_ENABLE_LOGGING
#define ROVER_INFO(msg) rover::Logger::GetInstance()->GetStream() <<"<Info>\n" \
  <<"  message: "<< msg <<"\n  file: " <<__FILE__<<"\n  line:  "<<__LINE__<<std::endl;
#define ROVER_WARN(msg) rover::Logger::GetInstance()->GetStream() <<"<Warn>\n" \
  <<"  message: "<< msg <<"\n  file: " <<__FILE__<<"\n  line:  "<<__LINE__<<std::endl;
#define ROVER_ERROR(msg) rover::Logger::GetInstance()->GetStream() <<"<Error>\n" \
  <<"  message: "<< msg <<"\n  file: " <<__FILE__<<"\n  line:  "<<__LINE__<<std::endl;

#define ROVER_DATA_OPEN(name) rover::DataLogger::GetInstance()->OpenLogEntry(name);
#define ROVER_DATA_CLOSE(time) rover::DataLogger::GetInstance()->CloseLogEntry(time);
#define ROVER_DATA_ADD(key,value) rover::DataLogger::GetInstance()->AddLogData(key, value);

#else
#define ROVER_INFO(msg)
#define ROVER_WARN(msg)
#define ROVER_ERROR(msg)

#define ROVER_DATA_OPEN(name)
#define ROVER_DATA_CLOSE(name)
#define ROVER_DATA_ADD(key,value)
#endif
}; // namespace vtkh

#endif //VTK_H_LOGGER_HPP
