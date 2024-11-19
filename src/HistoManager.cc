// the HistoManager is used to store the information from the Run, Event, Step, Track
// in the form of ROOT ntuples. 
// As an example, here I'm saving the kinetic energy, particle mass and z position from the step. 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//#include "g4root.hh"
#include "G4AnalysisManager.hh"


HistoManager::HistoManager():fFactoryOn(false)
{}

HistoManager::~HistoManager() {}
/*
void HistoManager::Book()
{
        G4AnalysisManager* ana = G4AnalysisManager::Instance();
	ana->SetFileName("Cosmic_MMs");
        ana->SetVerboseLevel(1);
        //ana->SetNtupleMerging(true);
        //ana->SetActivation(true);
        //ana->SetNtupleDirectoryName("Cosmics");

        G4bool fileOpen = ana->OpenFile();
        if (!fileOpen)
        {
                G4cerr << "\n ---> HistoManager::Book: cannot open "
                        <<ana->GetFileName() << G4endl;
                return;
        }


	ana->CreateNtuple("monopole", "MAGCOS Monopoles");
        ana->CreateNtupleDColumn(0, "zTurningPoint");
        ana->CreateNtupleDColumn(0, "kineticEnergy");
        ana->CreateNtupleDColumn(0, "ParticleMass");
        ana->FinishNtuple(0);

        fFactoryOn = true;
        G4cout << " --->> Ntuples created "
                << ana->GetFileName() << "." << ana->GetFileType() << G4endl;

}


void HistoManager::Save()
{
        if (! fFactoryOn) return;
        auto ana = G4AnalysisManager::Instance();
        ana->Write();
        ana->CloseFile();

        G4cout << "\n --> Ntuple saved " << G4endl;
        delete G4AnalysisManager::Instance();
        fFactoryOn = false;
}



void HistoManager::FillNtuple(G4double zTurningPoint, G4double kineticEnergy, G4double ParticleMass)
{
        auto ana = G4AnalysisManager::Instance();
        ana->FillNtupleDColumn(0, 0, zTurningPoint);
        ana->FillNtupleDColumn(0, 1, kineticEnergy);
        ana->FillNtupleDColumn(0, 2, ParticleMass);
        ana->AddNtupleRow(0);

}

*/
