#include "VertexMCOperator.hh"
using EVENT::Vertex;
using std::vector;
using std::string;
using EVENT::MCParticle;
using IMPL::VertexImpl;
using IMPL::ReconstructedParticleImpl;
namespace TTbarAnalysis 
{
	VertexMCOperator:: VertexMCOperator ()
	{
	}
	vector< Vertex * > * VertexMCOperator::Construct(DecayChain * chain)
	{
		if (!chain || chain->GetSize() == 0) 
		{
			return NULL;
		}
		
		vector< Vertex * > * result = new vector< Vertex * >();
		const double * ip = chain->Get(0)->getVertex();
		for (int i = 1; i < chain->GetSize(); i++) // <<==============================
		{
			result->push_back(construct(chain->Get(i), ip, chain->GetParentPDG(), i+1));
		}
		for (int i = 0; i < result->size(); i++) 
		{
			addParticle(result->at(i), chain->Get(i));
		}
		return result;
	}

	Vertex * VertexMCOperator::construct(EVENT::MCParticle * particle, const double * ip, int pdg, int number)
	{
		VertexImpl * result = new VertexImpl();
		
		const double * initial;
			initial = particle->getVertex();
		float distance = MathOperator::getDistance(ip, initial);
		result->setPrimary(false);
		result->setAlgorithmType("VertexMCOperator");
		
		result->setPosition(initial[0], initial[1], initial[2]);
		result->addParameter (distance);
		result->addParameter (pdg);
		result->addParameter (number);

		return result;
	}
	void VertexMCOperator::addParticle(Vertex * vertex, MCParticle * particle)
	{
		if (!vertex || !particle) 
		{
			std::cout << "ERRORMC: argument is null!\n";
			return;
		}
		ReconstructedParticleImpl * reco = new ReconstructedParticleImpl();
		reco->setType(particle->getPDG());
		reco->setMass(particle->getMass());
		reco->setCharge(particle->getCharge());
		reco->setMomentum(particle->getMomentum());
		reco->setEnergy(particle->getEnergy());
		VertexImpl * ivertex = static_cast<VertexImpl*>(vertex);
		ivertex->setAssociatedParticle(reco);
		if (vertex->getAssociatedParticle()) 
		{
			std::cout << "Particle attached!\n";
		}
	}
} /* TTbarAnalysis */
