#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
//#include "g4root.hh"
#include "G4AnalysisManager.hh"


class HistoManager
{
        public:
                HistoManager();
                ~HistoManager();

        public:
                void Book();
                void Save();

                G4String fFileName;

                void FillNtuple(G4double zTurningPoint, G4double kineticEnergy, G4double ParticleMass);

	private:
                G4bool fFactoryOn;
};

#endif

