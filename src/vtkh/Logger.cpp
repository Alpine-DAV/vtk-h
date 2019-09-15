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

DataLogger DataLogger::Instance;

DataLogger::DataLogger()
  : AtBlockStart(true)
{
  Blocks.push(Block(0));
}

DataLogger::~DataLogger()
{
  //WriteLog();
  //std::cout<<"Instance "<<DataLogger::Instance<<"\n";
  std::cout<<Stream.str();
  Stream.str("");
}

DataLogger*
DataLogger::GetInstance()
{
  return &DataLogger::Instance;
}

DataLogger::Block&
DataLogger::CurrentBlock()
{
  return Blocks.top();
}

void
DataLogger::WriteIndent()
{
  int indent = this->CurrentBlock().Indent;
  bool listStart = this->CurrentBlock().AtListItemStart;

  if (listStart)
  {
    --indent;
  }

  for (int i = 0; i < indent; ++i)
  {
    Stream << "  ";
  }

  if (listStart)
  {
    Stream << "- ";
    CurrentBlock().AtListItemStart = false;
  }
}

void
DataLogger::WriteLog()
{
  std::stringstream log_name;
  std::ofstream stream;
  log_name<<"vtkh_data";
#ifdef ROVER_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  log_name<<"_"<<rank;
#endif
  log_name<<".log";
  stream.open(log_name.str().c_str(), std::ofstream::out);
  if(!stream.is_open())
  {
    std::cerr<<"Warning: could not open the vtkh data log file\n";
    return;
  }
  stream<<Stream.str();
  stream.close();
}

void
DataLogger::OpenLogEntry(const std::string &entryName)
{
    WriteIndent();
    Stream<<entryName<<":"<<"\n";
    int indent = this->CurrentBlock().Indent;
    Blocks.push(Block(indent+1));

    Timer timer;
    Timers.push(timer);
    AtBlockStart = true;

}
void
DataLogger::CloseLogEntry()
{
  WriteIndent();
  this->Stream<<"time : "<<Timers.top().elapsed()<<"\n";
  Timers.pop();
  Blocks.pop();
  AtBlockStart = false;
}
};
