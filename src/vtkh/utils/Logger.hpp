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
  void Write(const int level, const std::string &message, const char *file, int line);
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

#ifdef ENABLE_LOGGING
#define VTKH_INFO(msg) vtkh::Logger::GetInstance("info")->GetStream()<<msg<<std::endl;
#define VTKH_WARN(msg) vtkh::Logger::GetInstance("warning")->GetStream()<<msg<<std::endl;
#define VTKH_ERROR(msg) vtkh::Logger::GetInstance("error")->GetStream()<<msg<<std::endl;
#define VTKH_DATA_ADD(key,value) vtkh::DataLogger::GetInstance()->AddLogData(key, value);

#else
#define VTKH_INFO(msg)
#define VTKH_WARN(msg)
#define VTKH_ERROR(msg)
#define VTKH_DATA_ADD(key,value)
#endif

}; // namespace vtkh

#endif //VTK_H_LOGGER_HPP
