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
	  void Write(MCOperator & opera,DecayChain * chain, std::vector< Vertex * > * bvertexes);
	  void AddProngs( VertexMCOperator & vertexOperator, MCOperator & opera, DecayChain * chain, std::vector< Vertex * > * vertices);
	  void Write(std::vector< MCParticle * > particle , int v);
	  double getMissingPt(std::vector< MCParticle * > & bdaugthers, std::vector< MCParticle * > & cdaughters, Vertex * vertex);
	  void WriteQuarksCollection(LCEvent * evt, std::vector< MCParticle * > & quarks);
	  void ExtractStarParticles(LCEvent * evt, MCOperator opera,  DecayChain * bChainRaw, DecayChain * bChain, int v);
	  void ClearVariables(); 
	 protected:
	
	  /** Input collection name.
	   */
	  std::string _colName ;
	  std::string _outputcolName;
	  std::string _outputquarkcolName;
	  std::vector<MESONS> _pdgs;
	  std::string _outputBStarName; 
	  std::string _outputBbarStarName; 
	  int _tagParameter;
	  float _aParameter;
	  float _bParameter;
	  int _writeBonlyParameter;
	  TFile * _hfile;
	  TTree * _hTree;
	  TTree * _hVertexTree;
	  TTree * _hTrackTree;
	  TTree * _hBStarTree;
	  std::string _hfilename ;
	  
	  int _tag;
	  int _numberOfB0;
	  float _firstVertexDistance[2];
	  float _secondVertexDistance[2];
	  int _totalBcharge;
	  int _ccharge;
	  int _cbarcharge;
	  int _bcharge;
	  int _bbarcharge;
	  float _baccuracy;
	  float _bbaraccuracy;
	  float _bIPdistance;
	  float _bbarIPdistance;
	  float _btracks;
	  float _bbartracks;
	  float _ctracks;
	  float _cbartracks;
	  float _bdistance;
	  float _bbardistance;
	  float _bmomentum;
	  float _bbarmomentum;
	  float _cmomentum;
	  float _cbarmomentum;
	  float _caccuracy;
	  float _cbaraccuracy;
	  int _bnumber;
	  int _bbarnumber;
	  int _cnumber;
	  int _cbarnumber;
	  int _btotalnumber;
	  int _bbartotalnumber;
	  int _bnumber_f;
	  int _bbarnumber_f;
	  int _cnumber_f;
	  int _cbarnumber_f;
	  float _bptmiss;
	  float _bbarptmiss;

	
	  static const int MAXV = 15;
	  int _numberOfVertexes;
	  float _distanceFromIP[MAXV];
	  float _coordinates[MAXV][3];
	  int _PDG[MAXV];
	  int _generation[MAXV];
	  int _charge[MAXV];
	  int _numberOfParticles[MAXV];
	  float _energyOfParticles[MAXV][MAXV];
	  float _momentumOfParticles[MAXV][MAXV];
	  float _massOfParticles[MAXV][MAXV];
	
	  float _bptrack[MAXV];
	  float _betatrack[MAXV];
	  float _boffsettrack[MAXV];
	  float _bbarptrack[MAXV];
	  float _bbaretatrack[MAXV];
	  float _bbaroffsettrack[MAXV];
	  float _cptrack[MAXV];
	  float _cetatrack[MAXV];
	  float _coffsettrack[MAXV];
	  float _cbarptrack[MAXV];
	  float _cbaretatrack[MAXV];
	  float _cbaroffsettrack[MAXV];

	  int _bstarnumber;
	  float _bstarmomentum[MAXV];
	  float _bstaroffset[MAXV];
	  double ip[3];
	  int _nRun ;
	  int _nEvt ;
	} ;
} /* TTbarAnalysis */
#endif



