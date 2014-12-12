#include "MCOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
namespace TTbarAnalysis
{
	MCOperator:: MCOperator (LCCollection * col)
	{
		myCollection = col;
		myPrimaryVertex[0] = 0.0;
		myPrimaryVertex[1] = 0.0;
		myPrimaryVertex[2] = 0.0;
		myAngleCut = 0.7854; // \pi/4
	}
	DecayChain * MCOperator::Construct(string name, int pdg, vector<MESONS> typeOfProducts)
	{
		DecayChain * result = new DecayChain(name, pdg);
		int number = myCollection->getNumberOfElements();
		bool finished = false;
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) ) ;
			if (particle->getPDG() == pdg) 
			{
				if (CheckForUnification(particle, pdg)) 
				{
					std::cout << "WARNING: found a quark-antiquark merging!\n";
					continue;
				}
				for (int j = 0; j < typeOfProducts.size(); j++) 
				{
					MESONS currentPDG = typeOfProducts[j];
					MCParticle * found = FindYoungestChild(particle, currentPDG);
					if (found) 
					{
						if (CheckForColorString(found, pdg)) 
						{
							MCParticle * service = found->getParents()[0];
							MCParticle * consistent = GetConsistentDaughter(particle, service, currentPDG);
							if (consistent) 
							{
								result->Add(consistent);
								found = consistent;
								finished = true;
							}
						}
						else 
						{
							if (found->getParents()[0]->getPDG() == 92 && found->getCharge() * particle->getCharge() < -0.001) 
							{
								std::cout << "FATAL: Charge is wrong!\n";
							}
							result->Add(found);
							finished = true;
						}
					}
					else 
					{
						break;
					}
					particle =(finished)? found : NULL;
				}
			}
			if (finished)
			{
				break;
			}
		}
		return result;
	}
	bool MCOperator::CheckForUnification(MCParticle * particle, int pdg)
	{
		vector< MCParticle * > daughters = particle->getDaughters();
		if (daughters.size() == 0) 
		{
			return false;
		}
		if (daughters.size() == 1 && daughters[0]->getPDG() == pdg) 
		{
			return CheckForUnification(daughters[0], pdg);
		}
		if (daughters.size() == 1 && (daughters[0]->getPDG() == 92 || daughters[0]->getPDG() == 94)) 
		{
			vector< MCParticle * > grandDaughters = daughters[0]->getDaughters();
			bool foundParticle = false;
			bool foundAntiparticle = false;
			for (int i = 0; i < grandDaughters.size(); i++) 
			{
				if (grandDaughters[i]->getPDG() == pdg) 
				{
					foundParticle = true;
				}
				if (grandDaughters[i]->getPDG() == -pdg) 
				{
					foundAntiparticle = true;
				}
			}
			return foundParticle && foundAntiparticle;
		}
		return false;

	}
	bool MCOperator::CheckParticle(MCParticle * particle, MESONS type)
	{
		bool result = false;
		vector<int> pdg =  ConstantStorage::GET_PDG(type);
		for (int i = 0; i < pdg.size(); i++) 
		{
			if (abs(particle->getPDG()) == pdg[i]) 
			{
				result = true;
				break;
			}
		}
		return result;
	}
	bool MCOperator::CheckForColorString(EVENT::MCParticle * daughter, int pdgOfParent)
	{
		vector< MCParticle * > parents = daughter->getParents();
		if (parents.size() == 1) 
		{
			MCParticle * parent = parents[0];
			if (parent->getPDG() == 92) 
			{
				vector< MCParticle * > grandParents = parent->getParents();
				int count = 0;
				for (int i = 0; i < grandParents.size(); i++) 
				{
					if (abs(grandParents[i]->getPDG()) == abs(pdgOfParent)) 
					{
						count++;
					}
					if (count > 1) 
					{
						std::cout << "United color string found for PDG " << pdgOfParent <<"!\n";
						return true;
					}
				}
			}
		}
		return false;
	}
	MCParticle * MCOperator::GetConsistentDaughter(MCParticle * parent, MCParticle * service, MESONS type)
	{
		vector< MCParticle * > daughters = SelectDaughtersOfType(service, type); //service->getDaughters();
		std::cout << "Parent PDG: " << parent->getPDG() << ";\n";
		int size = daughters.size();
		float angle[size];
		if (size == 1) 
		{
			 //std::cout << "We have only 1 meson of type " << type << '\n';
		}
		for (int i = 0; i < size; i++) 
		{
			MCParticle * daughter = daughters[i];
			std::cout << "Daughter PDG: " << daughter->getPDG() << ";\n";
			if (daughter->getEnergy() > parent->getEnergy())// || 
			{
				continue;
			}
			if (daughter->getCharge()*parent->getCharge() > 0.0)
			{
				std::cout << "Same sign of charge!\n";
				return daughter;
			}
			angle[i] = MathOperator::getAngle(daughter->getMomentum(), parent->getMomentum());
		}
		float minAngle = myAngleCut;
		int winner = -1;
		std::cout << "Checking angles...\n";
		for (int i = 0; i < size; i++) 
		{
			std::cout << "Angle " << i << ": " << angle[i] << '\n';
			if (angle[i] < minAngle && 
			   daughters[i]->getEnergy() < parent->getEnergy() &&
			   daughters[i]->getCharge()*parent->getCharge() > -0.0001) 
			{
				minAngle = angle[i];
				winner = i;
			}
		}
		if (winner > -1) 
		{
			std::cout << "Choosing angle " << winner <<  "\n";
			if (daughters[winner]->getCharge()*parent->getCharge() < -0.0001) 
			{
				std::cout << "FATAL: Charge is wrong!\n";
				
			}
			return daughters[winner];
		}
		std::cout << "Angle is wrong!\n";
		return NULL;
	}
	MESONS MCOperator::GetParticleType(MCParticle * particle)
	{
		MESONS result = EXCEPTIONAL_MESONS;
		int count = _max_MESONS;
		for (int i = 1; i < count; i++) 
		{
			MESONS candidate = static_cast<MESONS>(i);
			vector<int> pdg =  ConstantStorage::GET_PDG(candidate);
			for (int j = 0; j < pdg.size(); j++) 
			{
				if (abs(particle->getPDG()) == pdg[j]) 
				{
					result = candidate;
					break;
				}
			}
		}
		return result;

	}
	MCParticle * MCOperator::FindExceptionalChild(MCParticle * parent, MESONS parentType)
	{
		if (!parent) 
		{
			//std::cout << "ERROR: Input particle is NULL!\n";
			return NULL;
		}
		vector< MCParticle * > daughters = parent->getDaughters();
		if (daughters.size() < 1) 
		{
			//std::cout<<"ERROR: Found 0 daughters!\n";
			return NULL;
		}
		if (parentType == EXCEPTIONAL_MESONS) 
		{
			parentType = GetParticleType(parent);
			if (parentType == EXCEPTIONAL_MESONS) 
			{
				std::cout << "ERROR: Probably, parent is not meson, to be continue....\n";
				return NULL;
			}
		}
		bool found = false;
		for (int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (CheckParticle(daughter,parentType)) 
			{
				//std::cout << "Found PDG to exclude: " << daughter->getPDG() << "; Going recursion\n";
				found = true;
			}
		}
		if (!found) 
		{
			//std::cout << "Found state without initial mesons!\n";
			return daughters[0];
		}
		//std::cout<<"Not found in I gen. Checking next gen.\n";
		for (int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (daughter->getEnergy() < 1.0 && daughter->getMass() < 300.0) 
			{
				//std::cout<<"Daughter does not pass the selection cuts\n";
				continue;
			}
			return FindExceptionalChild(daughter, parentType);
		}
		return NULL;
	}
	MCParticle * MCOperator::FindYoungestChild(MCParticle * parent, MESONS type)
	{
		if (!parent) 
		{
			//std::cout << "ERROR: Input particle is NULL!\n";
			return NULL;
		}
		if (type == EXCEPTIONAL_MESONS) 
		{
			return FindExceptionalChild(parent, EXCEPTIONAL_MESONS);
		}
		vector< MCParticle * > daughters = parent->getDaughters();
		if (daughters.size() < 1) 
		{
			//std::cout<<"ERROR: Found 0 daughters!\n";
			return NULL;
		}
		//std::cout<<"Checking " << daughters.size() <<" I gen of daughters. \n";
		for (int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (CheckParticle(daughter,type)) 
			{
				//std::cout << "Found PDG: " << daughter->getPDG() << ";\n";
				return daughter;
			}
		}
		//std::cout<<"Not found in I gen. Checking next gen.\n";
		for (int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters.at(i);
			if (daughter->getEnergy() < 1.0 && daughter->getMass() < 300.0) 
			{
				//std::cout<<"Daughter does not pass the selection cuts\n";
				continue;
			}
			return FindYoungestChild(daughter, type);
		}
		return NULL;
	}
	vector< MCParticle * > MCOperator::SelectDaughtersOfType(MCParticle * parent, MESONS type)
	{
		vector< MCParticle * > daughters = parent->getDaughters();
		vector< MCParticle * > result;
		for (int i = 0; i < daughters.size(); i++) 
		{
			if (CheckParticle(daughters[i], type)) 
			{
				result.push_back(daughters[i]);
			}
		}
		return result;
	}
	vector< MCParticle * > MCOperator::ScanForVertexParticles(const double * vertex, double precision)
	{
		int number = myCollection->getNumberOfElements();
		vector< MCParticle * > result;
		if (MathOperator::approximatelyEqual(myPrimaryVertex, vertex, precision)) 
		{
			 std::cout<<"ERROR: Vertex too close to primary!\n";
			 return result;
		}
		for (int i = 0; i < number; i++)
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) ) ;
			const double * another = particle->getVertex();
			if (MathOperator::approximatelyEqual(another, vertex, precision)) 
			{
				//std::cout<<"Found particle " << particle->getPDG() << " from this vertex! Checking its visability...\n";
				if (particle->isDecayedInCalorimeter() && abs(particle->getCharge())>0.00001) 
				{
					//std::cout<<"Particle reached calorimeter!\n";
					result.push_back(particle);
				}
				else 
				{
					//std::cout<<"Not stable, checking daughters:\n";
					vector< MCParticle * > daughters = SelectStableCloseDaughters(particle);
					for (int j = 0; j < daughters.size(); j++) 
					{
						result.push_back(daughters[j]);
					}
				}
			}
		}
		return result;
	}
	
	vector< MCParticle * > MCOperator::ScanForVertexParticles(const float * vertex, double precision)
	{
		int number = myCollection->getNumberOfElements();
		vector< MCParticle * > result;
		if (MathOperator::approximatelyEqual(myPrimaryVertex, vertex, precision)) 
		{
			 std::cout<<"ERROR: Vertex too close to primary!\n";
			 return result;
		}
		for (int i = 0; i < number; i++)
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) ) ;
			const double * another = particle->getVertex();
			if (MathOperator::approximatelyEqual(another, vertex, precision)) 
			{
				//std::cout<<"Found particle " << particle->getPDG() << " from this vertex! Checking its visability...\n";
				if ((particle->isDecayedInCalorimeter() || particle->isDecayedInTracker()) && abs(particle->getCharge())>0.00001) 
				{
					//std::cout<<"Particle reached calorimeter!\n";
					result.push_back(particle);
				}
				else 
				{
					//std::cout<<"Not stable, checking daughters:\n";
					vector< MCParticle * > daughters = SelectStableCloseDaughters(particle);
					for (int j = 0; j < daughters.size(); j++) 
					{
						result.push_back(daughters[j]);
					}
				}
			}
		}
		return result;
	}
	bool MCOperator::CheckProcessForPair(int pdg)
	{
		pdg = abs(pdg);
		if (pdg < 1) 
		{
			return false;
		}
		int countParticle = 0;
		int countAntiparticle = 0;
		int number = myCollection->getNumberOfElements();
		for (int i = 0; i < number; i++) 
		{
			MCParticle * particle = dynamic_cast<MCParticle*>( myCollection->getElementAt(i) );
			if (particle->getPDG() == pdg) 
			{
				countParticle++;
			}
			if (particle->getPDG() == 0-pdg) 
			{
				countAntiparticle++;
			}
		}
		std::cout<<"INFO: " << countParticle << " t-quarks and " << countAntiparticle << " tbar-quarks\n";
		return countParticle && countAntiparticle;
	}
	vector< MCParticle * > MCOperator::SelectStableCloseDaughters(MCParticle * parent)
	{
		vector< MCParticle * > result;
		if (!parent || 
		    parent->getPDG() == 22 ||
		    parent->getPDG() == 111 ||
		    CheckParticle(parent, CHARMED_MESONS)) 
		{
			return result;
		}
		vector< MCParticle * > daughters = parent->getDaughters();
		if (daughters.size()<1) 
		{
			return result;
		}
		if ( parent && MathOperator::getDistance(daughters[0]->getVertex(), parent->getVertex()) > 50.0) // Long distance!!!
		{
			//std::cout<<"WARNING: Long-lived noninteracting particle!\n";
			return result;
		}
		bool finished = false;
		for (int i = 0; i < daughters.size(); i++) 
		{
			MCParticle * daughter = daughters[i];
			//std::cout<<"\tFound daughter " << daughter->getPDG() << "\n";
			if ((daughter->isDecayedInCalorimeter() || daughter->isDecayedInTracker()) && 
			    abs(daughter->getCharge())>0.00001)
			{
				//std::cout<<"\tDaughter reached calorimeter!\n";
				finished = true;
				result.push_back(daughter);
			}
			else 
			{
				vector< MCParticle * > grannies = SelectStableCloseDaughters(daughter);
				for (int j = 0; j < grannies.size(); j++) 
				{
					result.push_back(grannies[j]);
				}
			}
		}
		return result;
	}
} /*  */
