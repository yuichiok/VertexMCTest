#include "VertexMCOperator.hh"
using EVENT::Vertex;
using std::vector;
using std::string;
using EVENT::MCParticle;
using EVENT::ReconstructedParticle;
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
	void VertexMCOperator::AddProngs(Vertex * vertex, vector< MCParticle * > & particles)
	{
		if (!vertex || particles.size() == 0) 
		{
			std::cout << "ERRORMC: argument is null!\n";
			return;
		}
		ReconstructedParticle * reco = vertex->getAssociatedParticle();
		for (unsigned int i = 0; i < particles.size(); i++) 
		{
			ReconstructedParticle * prong = translate(particles[i]);
			reco->addParticle(prong);
		}
		//std::cout << "Added " << reco->getParticles().size() << " particles!\n";

	}
	void VertexMCOperator::addParticle(Vertex * vertex, MCParticle * particle)
	{
		if (!vertex || !particle) 
		{
			std::cout << "ERRORMC: argument is null!\n";
			return;
		}
		ReconstructedParticle * reco = translate(particle);
		VertexImpl * ivertex = static_cast<VertexImpl*>(vertex);
		ivertex->setAssociatedParticle(reco);
	}
	ReconstructedParticle * VertexMCOperator::translate(EVENT::MCParticle * particle)
	{
		ReconstructedParticleImpl * reco = new ReconstructedParticleImpl();
		reco->setType(particle->getPDG());
		reco->setMass(particle->getMass());
		reco->setCharge(particle->getCharge());
		reco->setMomentum(particle->getMomentum());
		reco->setEnergy(particle->getEnergy());
		return reco;
	}
} /* TTbarAnalysis */
