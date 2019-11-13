#ifndef VTK_H_PARTICLE_ADVECTION_TASK_HPP
#define VTK_H_PARTICLE_ADVECTION_TASK_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/StatisticsDB.hpp>
#include <vtkh/utils/ThreadSafeContainer.hpp>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkh/filters/communication/BoundsMap.hpp>

#ifdef VTKH_ENABLE_LOGGING
#define DBG(msg) vtkh::Logger::GetInstance("out")->GetStream()<<msg
#define WDBG(msg) vtkh::Logger::GetInstance("wout")->GetStream()<<msg
#else
#define DBG(msg)
#define WDBG(msg)
#endif

namespace vtkh
{
template <typename ResultT>
class ParticleAdvectionTask
{
public:
    ParticleAdvectionTask(MPI_Comm comm, const vtkh::BoundsMap &bmap, ParticleAdvection *pa) :
        numWorkerThreads(-1),
        done(false),
        begin(false),
        communicator(comm, bmap),
        boundsMap(bmap),
        filter(pa),
        sleepUS(100)
    {
        m_Rank = vtkh::GetMPIRank();
        m_NumRanks = vtkh::GetMPISize();
        communicator.RegisterMessages(2, std::min(64, m_NumRanks-1), 128, std::min(64, m_NumRanks-1));
        ADD_TIMER("worker_sleep");
        ADD_COUNTER("worker_naps");
    }
    ~ParticleAdvectionTask()
    {
    }

    void Init(const std::vector<Particle> &particles, int N, int _sleepUS)
    {
        numWorkerThreads = 1;
        TotalNumParticles = N;
        sleepUS = _sleepUS;
        active.Assign(particles);
        inactive.Clear();
        terminated.Clear();
    }

    bool CheckDone()
    {
        bool val;
        stateLock.Lock();
        val = done;
        stateLock.Unlock();
        return val;
    }
    void SetDone()
    {
        stateLock.Lock();
        done = true;
        stateLock.Unlock();
    }

    bool GetBegin()
    {
        bool val;
        stateLock.Lock();
        val = begin;
        stateLock.Unlock();
        return val;
    }

    void SetBegin()
    {
        stateLock.Lock();
        begin = true;
        stateLock.Unlock();
    }

    void Go()
    {
        DBG("Go_bm: "<<boundsMap<<std::endl);
        DBG("actives= "<<active<<std::endl);

#ifdef VTKH_USE_OPENMP
        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            #pragma omp parallel num_threads(1)
            {
                this->Manage();
            }

            #pragma omp section
            #pragma omp parallel num_threads(numWorkerThreads)
            #pragma omp master
            {
                this->Work();
            }
        }
#else
        workerThreads.push_back(std::thread(ParticleAdvectionTask::Worker, this));
        this->Manage();
        for (auto &t : workerThreads)
            t.join();
#endif
    }

#ifndef VTKH_USE_OPENMP
    static void Worker(ParticleAdvectionTask *t)
    {
      t->Work();
    }
#endif

    void Work()
    {
      std::vector<ResultT> traces;

        while (!CheckDone())
        {
            std::vector<Particle> particles;
            if (active.Get(particles))
            {
                std::vector<Particle> I, T, A;

                DataBlockIntegrator *blk = filter->GetBlock(particles[0].blockIds[0]);

                TIMER_START("advect");
                WDBG("WORKER: Integrate "<<particles<<" --> "<<std::endl);
                int n = filter->InternalIntegrate<ResultT>(*blk, particles, I, T, A, traces);
                TIMER_STOP("advect");
                COUNTER_INC("advectSteps", n);
                WDBG("TIA: "<<T<<" "<<I<<" "<<A<<std::endl<<std::endl);

                worker_terminated.Insert(T);
                worker_active.Insert(A);
                worker_inactive.Insert(I);
            }
            else
            {
                TIMER_START("worker_sleep");
                usleep(sleepUS);
                TIMER_STOP("worker_sleep");
                COUNTER_INC("worker_naps", 1);
            }
        }
        WDBG("WORKER is DONE"<<std::endl);
        results.Insert(traces);
    }

    void Manage()
    {
        DBG("manage_bm: "<<boundsMap<<std::endl);

        int N = 0;

        DBG("Begin TIA: "<<terminated<<" "<<inactive<<" "<<active<<std::endl);
        MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());

        while (true)
        {
            DBG("MANAGE TIA: "<<terminated<<" "<<worker_inactive<<" "<<active<<std::endl<<std::endl);
            std::vector<Particle> out, in, term;
            worker_inactive.Get(out);
            worker_terminated.Get(term);

            int numTermMessages;
            communicator.Exchange(out, in, term, numTermMessages);
            int numTerm = term.size() + numTermMessages;

            if (!in.empty())
                active.Insert(in);
            if (!term.empty())
                terminated.Insert(term);

            N += numTerm;
            if (N > TotalNumParticles)
                throw "Particle count error";
            if (N == TotalNumParticles)
                break;

            if (active.Empty())
            {
                TIMER_START("sleep");
                usleep(sleepUS);
                TIMER_STOP("sleep");
                COUNTER_INC("naps", 1);
                communicator.CheckPendingSendRequests();
            }
        }
        DBG("TIA: "<<terminated<<" "<<inactive<<" "<<active<<" WI= "<<worker_inactive<<std::endl);
        DBG("RESULTS= "<<results.Size()<<std::endl);
        DBG("DONE_"<<m_Rank<<" "<<terminated<<" "<<active<<" "<<inactive<<std::endl);
        SetDone();
    }

    int m_Rank, m_NumRanks;
    int TotalNumParticles;

#ifndef VTKH_USE_OPENMP
    std::vector<std::thread> workerThreads;
#endif

    using ParticleList = vtkh::ThreadSafeContainer<Particle, std::vector>;
    using ResultsVec = vtkh::ThreadSafeContainer<ResultT, std::vector>;

    ParticleMessenger communicator;
    ParticleList active, inactive, terminated;
    ParticleList worker_active, worker_inactive, worker_terminated;
    ResultsVec results;

    int numWorkerThreads;
    int sleepUS;

    bool done, begin;
    vtkh::Mutex stateLock;
    BoundsMap boundsMap;
    ParticleAdvection *filter;
};
}

#endif //VTK_H_PARTICLE_ADVECTION_TASK_HPP
