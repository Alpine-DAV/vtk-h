#ifndef VTK_H_INTEGRATOR_HPP
#define VTK_H_INTEGRATOR_HPP

//#include "adapter.h"

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

#include <vtkh/vtkh_exports.h>
#include <vtkh/filters/Particle.hpp>
#include <vtkh/utils/ThreadSafeContainer.hpp>


class VTKH_API Integrator
{
    typedef vtkm::Float64 FieldType;
    typedef vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, 3>> FieldHandle;

    using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldHandle>;
    using RK4Type = vtkm::worklet::particleadvection::RK4Integrator<GridEvalType>;

public:
    Integrator(vtkm::cont::DataSet *ds, const std::string &fieldName, FieldType _stepSize, int _batchSize, int _rank)
      : stepSize(_stepSize), batchSize(_batchSize), rank(_rank)
    {
        vecField = ds->GetField(fieldName).GetData().Cast<FieldHandle>();
        gridEval = GridEvalType(ds->GetCoordinateSystem(), ds->GetCellSet(), vecField);
        rk4 = RK4Type(gridEval, stepSize);
    }

    int Advect(std::vector<vtkh::Particle> &particles,
               vtkm::Id MaxSteps,
               std::vector<vtkh::Particle> &I,
               std::vector<vtkh::Particle> &T,
               std::vector<vtkh::Particle> &A,
               std::vector<vtkm::worklet::ParticleAdvectionResult> *particleTraces,
               vtkh::ThreadSafeContainer<vtkh::Particle, std::vector> &workerInactive,
               vtkh::StatisticsDB& statsDB);

    int Advect(std::vector<vtkh::Particle> &particles,
               vtkm::Id maxSteps,
               std::vector<vtkh::Particle> &I,
               std::vector<vtkh::Particle> &T,
               std::vector<vtkh::Particle> &A,
               std::vector<vtkm::worklet::ParticleAdvectionResult> *particleTraces=NULL)
    {
        size_t nSeeds = particles.size();
        vtkm::cont::ArrayHandle<vtkm::Particle> seedArray;

        int steps0 = SeedPrep(particles, seedArray);

        vtkm::worklet::ParticleAdvection particleAdvection;
        vtkm::worklet::ParticleAdvectionResult result;

        result = particleAdvection.Run(rk4, seedArray, maxSteps);
        auto parPortal = result.Particles.GetPortalConstControl();

        //Update particle data.
        //Need a functor to do this...
        int steps1 = 0;
        for (int i = 0; i < nSeeds; i++)
        {
            particles[i].p = parPortal.Get(i);
            UpdateParticle(particles[i], I,T,A);
            steps1 += particles[i].p.NumSteps;
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

    int Trace(std::vector<vtkh::Particle> &particles,
              vtkm::Id maxSteps,
              std::vector<vtkh::Particle> &I,
              std::vector<vtkh::Particle> &T,
              std::vector<vtkh::Particle> &A,
              std::vector<vtkm::worklet::StreamlineResult> *particleTraces=NULL)
    {
        size_t nSeeds = particles.size();
        vtkm::cont::ArrayHandle<vtkm::Particle> seedArray;

        int steps0 = SeedPrep(particles, seedArray);

        vtkm::worklet::Streamline streamline;
        vtkm::worklet::StreamlineResult result;
        result = streamline.Run(rk4, seedArray, maxSteps);
        auto parPortal = result.Particles.GetPortalConstControl();

        //Update particle data.
        int steps1 = 0;
        for (int i = 0; i < nSeeds; i++)
        {
            /*
            vtkm::cont::ArrayHandle<vtkm::Id> ids;
            result.PolyLines.GetIndices(i, ids);
            auto idPortal = ids.GetPortalConstControl();
            vtkm::Id nPts = idPortal.GetNumberOfValues();
            */

            particles[i].p = parPortal.Get(i);
            steps1 += particles[i].p.NumSteps;

            UpdateParticle(particles[i], I,T,A);

#if 0
            {
                for (int j = 0; j < nPts; j++)
                {
                    vtkm::Vec3f p(parPortal.Get(idPortal.Get(j))[0],
                                  parPortal.Get(idPortal.Get(j))[1],
                                  parPortal.Get(idPortal.Get(j))[2]);
//                                         particles[i].id);
//                                         particles[i].blockId);

                    particleTraces->push_back(p);
                }
            }
#endif

        }
        if (particleTraces)
          (*particleTraces).push_back(result);

        int totalSteps = steps1-steps0;
        return totalSteps;
    }

private:

    void UpdateParticle(vtkh::Particle &p,
                        std::vector<vtkh::Particle> &I,
                        std::vector<vtkh::Particle> &T,
                        std::vector<vtkh::Particle> &A)
    {
      if (p.p.Status.CheckTerminate())
        T.push_back(p);
      else if (p.p.Status.CheckSpatialBounds())
        I.push_back(p);
      else if (p.p.Status.CheckOk())
        A.push_back(p);
      else
        T.push_back(p);
    }

    int SeedPrep(const std::vector<vtkh::Particle> &particles,
                 vtkm::cont::ArrayHandle<vtkm::Particle> &seedArray)
    {
        int stepsTaken = 0;
        size_t nSeeds = particles.size();
        seedArray.Allocate(nSeeds);
        auto seedPortal = seedArray.GetPortalControl();
        for (int i = 0; i < nSeeds; i++)
        {
            seedPortal.Set(i, particles[i].p);
            stepsTaken += particles[i].p.NumSteps;
        }
        return stepsTaken;
    }

    int rank;
    int batchSize;
    FieldType stepSize;
    GridEvalType gridEval;
    RK4Type rk4;
    FieldHandle vecField;
};

#endif
