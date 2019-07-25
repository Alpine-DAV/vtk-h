#include <vtkh/vtkh.hpp>
#include <vtkh/Logger.hpp>

namespace vtkh
{
std::map<std::string, Logger*> Logger::Loggers;

Logger::Logger(const std::string& name)
{
  std::stringstream logName;
  logName<<name;
#ifdef VTKH_PARALLEL
  logName<<"."<<vtkh::GetMPIRank();
#endif
  logName<<".log";

  Stream.open(logName.str().c_str(), std::ofstream::out);
  if(!Stream.is_open())
    std::cout<<"Warning: could not open the vtkh log file\n";
}

Logger::~Logger()
{
  if (Stream.is_open())
    Stream.close();
}

Logger* Logger::GetInstance(const std::string& name)
{
  if (Loggers.find(name) == Loggers.end())
    Loggers[name] = new Logger(name);

  return Loggers[name];
}

void
Logger::Write(const int level, const std::string &message, const char *file, int line)
{
  if(level == 0)
    Stream<<"<Info> \n";
  else if (level == 1)
    Stream<<"<Warning> \n";
  else if (level == 2)
    Stream<<"<Error> \n";
  Stream<<"  message: "<<message<<" \n  file: "<<file<<" \n  line: "<<line<<"\n";
}

// ---------------------------------------------------------------------------------------

DataLogger* DataLogger::Instance  = NULL;

DataLogger::DataLogger()
{
}

DataLogger::~DataLogger()
{
  Stream.str("");
}

DataLogger*
DataLogger::GetInstance()
{
  if(DataLogger::Instance == NULL)
  {
    DataLogger::Instance =  new DataLogger();
  }
  return DataLogger::Instance;
}

void
DataLogger::WriteLog()
{
  std::stringstream log_name;
  std::ofstream stream;
  log_name<<"rover_data";
#ifdef ROVER_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  log_name<<"_"<<rank;
#endif
  log_name<<".log";
  stream.open(log_name.str().c_str(), std::ofstream::out);
  if(!stream.is_open())
  {
    std::cerr<<"Warning: could not open the rover data log file\n";
    return;
  }
  stream<<Stream.str();
  stream.close();
}

void
DataLogger::OpenLogEntry(const std::string &entryName)
{
    Stream<<entryName<<" "<<"<\n";
    Entries.push(entryName);
}
void
DataLogger::CloseLogEntry(const double &entryTime)
{
  this->Stream<<"total_time "<<entryTime<<"\n";
  this->Stream<<this->Entries.top()<<" >\n";
  Entries.pop();
}
};
