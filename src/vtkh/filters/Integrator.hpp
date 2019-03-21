#ifndef VTK_H_INTEGRATOR_HPP
#define VTK_H_INTEGRATOR_HPP

#include "adapter.h"

#include <list>
#include <vector>
#include <deque>
#include <vector>
#include <string>

#include <vtkm/cont/DataSet.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/Particles.h>

#include <vtkh/filters/Particle.hpp>

using namespace std;

class Integrator
{
    typedef WORKER_DEVICE_ADAPTER DeviceAdapter;
    typedef vtkm::Float64 FieldType;
    typedef vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, 3>> FieldHandle;
    typedef FieldHandle::template ExecutionTypes<DeviceAdapter>::PortalConst FieldPortalConstType;

    using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldHandle>;
    using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;

public:
    Integrator(vtkm::cont::DataSet &ds, const string &fieldName, FieldType _stepSize) : stepSize(_stepSize)
    {
        vecField = ds.GetField(fieldName).GetData().Cast<FieldHandle>();

        gridEval = GridEvalType(ds.GetCoordinateSystem(), ds.GetCellSet(), vecField);
        rk4 = RK4Type(gridEval, stepSize);
    }

    ~Integrator()
    {
        //coords.ReleaseResourcesExecution();
        //cells.ReleaseResourcesExecution();
        //vecField.ReleaseResourcesExecution();
        //cout<<"Toss an integrator: refCnt= "<<vecField.Internals.use_count()<<endl;
    }

    /*
    int Go(bool recordPath,
           std::vector<Particle> &particles,
           const vtkm::Id &maxSteps,
           std::list<Particle> &I,
           std::list<Particle> &T,
           std::list<Particle> &A,
           std::vector<vtkm::worklet::StreamlineResult<FieldType>> *particleTraces=NULL)
    {
        if (recordPath)
            return Trace(particles, maxSteps, I, T, A, particleTraces);
        else
            return Advect(particles, maxSteps, I, T, A, particleTraces);
    }
    */

    int Advect(std::vector<Particle> &particles,
               const vtkm::Id &maxSteps,
               std::list<Particle> &I,
               std::list<Particle> &T,
               std::list<Particle> &A,
               std::vector<vtkm::worklet::ParticleAdvectionResult> *particleTraces=NULL)
    {
        size_t nSeeds = particles.size();
        vtkm::cont::ArrayHandle<vtkm::Vec<FieldType,3>> seedArray;
        vtkm::cont::ArrayHandle<vtkm::Id> stepsTakenArray;

        int steps0 = SeedPrep(particles, seedArray, stepsTakenArray);

        vtkm::worklet::ParticleAdvection particleAdvection;
        vtkm::worklet::ParticleAdvectionResult result;

        result = particleAdvection.Run(rk4, seedArray, stepsTakenArray, maxSteps);
        auto posPortal = result.positions.GetPortalConstControl();
        auto statusPortal = result.status.GetPortalConstControl();
        auto stepsPortal = result.stepsTaken.GetPortalConstControl();

        //Update particle data.
        //Need a functor to do this...
        int steps1 = 0;
        for (int i = 0; i < nSeeds; i++)
        {
            particles[i].coords = posPortal.Get(i);
            particles[i].nSteps = stepsPortal.Get(i);
            UpdateStatus(particles[i], statusPortal.Get(i), maxSteps, I,T,A);
            steps1 += stepsPortal.Get(i);
        }

        /*
        {
            for (int i = 0; i < nSeeds; i++)
            {
                vtkm::Vec<FieldType,4> p(posPortal.Get(i)[0],
                                     posPortal.Get(i)[1],
                                     posPortal.Get(i)[2],
                                     particles[i].id);
//                                     particles[i].blockId);
                particleTraces->push_back(p);
            }
        }
        */

        if (particleTraces)
          (*particleTraces).push_back(result);

        int totalSteps = steps1-steps0;
        return totalSteps;
    }

    int Trace(std::vector<Particle> &particles,
              const vtkm::Id &maxSteps,
              std::list<Particle> &I,
              std::list<Particle> &T,
              std::list<Particle> &A,
              std::vector<vtkm::worklet::StreamlineResult> *particleTraces=NULL)
    {
        size_t nSeeds = particles.size();
        vtkm::cont::ArrayHandle<vtkm::Vec<FieldType,3>> seedArray;
        vtkm::cont::ArrayHandle<vtkm::Id> stepsTakenArray;

        int steps0 = SeedPrep(particles, seedArray, stepsTakenArray);

        vtkm::worklet::Streamline streamline;
        vtkm::worklet::StreamlineResult result;
        result = streamline.Run(rk4, seedArray, stepsTakenArray, maxSteps);

        auto posPortal = result.positions.GetPortalConstControl();
        auto statusPortal = result.status.GetPortalConstControl();
        auto stepsPortal = result.stepsTaken.GetPortalConstControl();

        //Update particle data.
        int steps1 = 0;
        for (int i = 0; i < nSeeds; i++)
        {
            vtkm::cont::ArrayHandle<vtkm::Id> ids;
            result.polyLines.GetIndices(i, ids);
            auto idPortal = ids.GetPortalConstControl();
            vtkm::Id nPts = idPortal.GetNumberOfValues();

            particles[i].coords = posPortal.Get(idPortal.Get(nPts-1));
            particles[i].nSteps = stepsPortal.Get(i);
            steps1 += stepsPortal.Get(i);
            UpdateStatus(particles[i], statusPortal.Get(i), maxSteps, I,T,A);

            /*
            {
                for (int j = 0; j < nPts; j++)
                {
                    vtkm::Vec<FieldType,4> p(posPortal.Get(idPortal.Get(j))[0],
                                         posPortal.Get(idPortal.Get(j))[1],
                                         posPortal.Get(idPortal.Get(j))[2],
                                         particles[i].id);
//                                         particles[i].blockId);

                    particleTraces->push_back(p);
                }
            }
            */
        }
        if (particleTraces)
          (*particleTraces).push_back(result);

        int totalSteps = steps1-steps0;
        return totalSteps;
    }

private:

    void UpdateStatus(Particle &p,
                      const vtkm::Id &status,
                      const vtkm::Id &maxSteps,
                      std::list<Particle> &I,
                      std::list<Particle> &T,
                      std::list<Particle> &A)
    {
        /*
        p.status = Particle::TERMINATE;
        T.push_back(p);
        return;
        */

        if (p.nSteps >= maxSteps || status == vtkm::worklet::particleadvection::ParticleStatus::TERMINATED)
        {
            p.status = Particle::TERMINATE;
            T.push_back(p);
        }
        else if (status == vtkm::worklet::particleadvection::ParticleStatus::STATUS_OK)
        {
            p.status = Particle::ACTIVE;
            A.push_back(p);
        }
        else
        {
            p.status = Particle::OUTOFBOUNDS;
            I.push_back(p);
        }
    }

    int SeedPrep(const std::vector<Particle> &particles,
                 vtkm::cont::ArrayHandle<vtkm::Vec<FieldType,3>> &seedArray,
                 vtkm::cont::ArrayHandle<vtkm::Id> &stepArray)
    {
        int stepsTaken = 0;
        size_t nSeeds = particles.size();
        seedArray.Allocate(nSeeds);
        stepArray.Allocate(nSeeds);
        auto seedPortal = seedArray.GetPortalControl();
        auto stepPortal = stepArray.GetPortalControl();
        for (int i = 0; i < nSeeds; i++)
        {
            seedPortal.Set(i, particles[i].coords);
            stepPortal.Set(i, particles[i].nSteps);
            stepsTaken += particles[i].nSteps;
        }
        return stepsTaken;
    }

    FieldType stepSize;
    GridEvalType gridEval;
    RK4Type rk4;
    FieldHandle vecField;
};

#endif
