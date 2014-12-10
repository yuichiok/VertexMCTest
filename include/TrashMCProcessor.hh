#ifndef TrashMCProcessor_h
#define TrashMCProcessor_h 1
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include "marlin/VerbosityLevels.h"
#include "marlin/Processor.h"
#include "lcio.h"

// ----- include for verbosity dependend logging ---------
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <map>

#include "ConstantStorage.hh"
#include "MathOperator.hh"
#include "MCOperator.hh"
#include "VertexMCOperator.hh"
using namespace lcio ;
using namespace marlin ;


namespace TTbarAnalysis 
{
	
	class TrashMCProcessor : public Processor 
	{
	  
	 public:
	  
	  virtual Processor*  newProcessor() { return new TrashMCProcessor ; }
	  
	  
	  TrashMCProcessor() ;
	  
	  /** Called at the begin of the job before anything is read.
	   * Use to initialize the processor, e.g. book histograms.
	   */
	  virtual void init() ;
	  
	  /** Called for every run.
	   */
	  virtual void processRunHeader( LCRunHeader* run ) ;
	  
	  /** Called for every event - the working horse.
	   */
	  virtual void processEvent( LCEvent * evt ) ; 
	  
	  
	  virtual void check( LCEvent * evt ) ; 
	  
	  /** Called after data processing for clean up.
	   */
	  virtual void end() ;
	 /////////////////CUSTOM/////////////////////////// 
	  void PrintParticle(MCParticle * particle);
	  void PrintChain(std::vector< MCParticle * > * chain);

	  void WriteVertexCollection(LCEvent * evt, std::vector< Vertex * > * bvertexes, std::vector< Vertex * > * bbarvertexes);
	  void Write(MCOperator & opera, DecayChain * chain, int & number);
	  void ClearVariables(); 
	 protected:
	
	  /** Input collection name.
	   */
	  std::string _colName ;
	  std::string _outputcolName;
	  std::vector<MESONS> _pdgs;
	  
	  TFile * _hfile;
	  TTree * _hTree;
	  TTree * _hVertexTree;
	  std::string _hfilename ;

	  int _charmedPDG[2];
	  int _strangePDG[2];
	  int _bottomPDG[2];
	  float _strangeCharge[2];
	  int _numberOfB0;
	  float _firstVertexDistance[2];
	  float _secondVertexDistance[2];
	  float _bdistance;
	  float _bbardistance;
	  static const int MAXV = 15;
	  int _tag;
	  int _numberOfVertexes;
	  float _distanceFromIP[MAXV];
	  float _coordinates[MAXV][3];
	  int _PDG[MAXV];
	  int _generation[MAXV];
	  int _numberOfParticles[MAXV];
	  float _energyOfParticles[MAXV][MAXV];
	  float _momentumOfParticles[MAXV][MAXV];
	  float _massOfParticles[MAXV][MAXV];


	  int _nRun ;
	  int _nEvt ;
	} ;
} /* TTbarAnalysis */
#endif



