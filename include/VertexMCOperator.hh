#include <EVENT/MCParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <IMPL/VertexImpl.h>
#include "DecayChain.hh"
#include "ConstantStorage.hh"
#include "MathOperator.hh"
#ifndef _VertexMCOperator_hh
#define _VertexMCOperator_hh
namespace TTbarAnalysis
{
	class VertexMCOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			VertexMCOperator ();
			virtual ~VertexMCOperator () {};
		//
		//	Methods
		//
			std::vector< EVENT::Vertex * > * Construct(DecayChain * chain);
			void AddProngs(EVENT::Vertex * vertex, std::vector< EVENT::MCParticle * > & particles);
		private:
		//
		//	Data
		//
			/* data */
		//
		//	Private methods
		//
			EVENT::Vertex * construct(EVENT::MCParticle * particle, const double * ip, int pdg, int number);
			void addParticle(EVENT::Vertex * vertex, EVENT::MCParticle * particle);
			EVENT::ReconstructedParticle * translate(EVENT::MCParticle * particle);
	};
}
#endif
