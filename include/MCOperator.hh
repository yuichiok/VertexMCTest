#include <stdlib.h>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <EVENT/MCParticle.h>
#include <EVENT/LCCollection.h>
#include "DecayChain.hh"
#include "ConstantStorage.hh"
#include "MathOperator.hh"
#ifndef _MCOperator_hh
#define _MCOperator_hh
namespace TTbarAnalysis 
{
	class MCOperator 
	{
		public:
		//
		//	Constants
		//
		
		//
		//	Constructors
		//
			MCOperator (EVENT::LCCollection * col);
			virtual ~MCOperator () {};
		//
		//	Methods
		//
			bool CheckProcessForPair(int pdg);
			DecayChain * Construct(std::string name, int pdg, std::vector<MESONS> typeOfProducts);
			bool CheckParticle(EVENT::MCParticle * particle, MESONS type);
			std::vector< EVENT::MCParticle * > SelectDaughtersOfType(EVENT::MCParticle * parent, MESONS type);
			MESONS GetParticleType(EVENT::MCParticle * particle); 
			EVENT::MCParticle * FindExceptionalChild(EVENT::MCParticle * parent, MESONS parentType);
			EVENT::MCParticle * FindYoungestChild(EVENT::MCParticle * parent, MESONS type);
			bool CheckForServiceParticle(EVENT::MCParticle * daughter, int pdgOfParent);
			EVENT::MCParticle * GetConsistentDaughter(EVENT::MCParticle * parent, EVENT::MCParticle * service, MESONS type);
			std::vector< EVENT::MCParticle * > ScanForVertexParticles(const float * vertex, double precision);
			std::vector< EVENT::MCParticle * > SelectStableCloseDaughters(EVENT::MCParticle * parent);
			std::vector< EVENT::MCParticle * > ScanForVertexParticles(const double * vertex, double precision);
		private:
		//
		//	Data
		//
			EVENT::LCCollection * myCollection;
			double myPrimaryVertex[3];
			float myAngleCut;
			int myCurrentParentPDG;
		//
		//	Private methods
		//
	};
} /* TTbarAnalysis */
#endif
