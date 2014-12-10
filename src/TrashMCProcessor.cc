#include "TrashMCProcessor.hh"
using std::string;
using std::vector;
using std::map;
namespace TTbarAnalysis
{
	TrashMCProcessor aTrashMCProcessor ;
	TrashMCProcessor::TrashMCProcessor() : Processor("TrashMCProcessor") 
	{

	    _description = "TrashMCProcessor does whatever it does ..." ;


	    registerInputCollection( LCIO::MCPARTICLE,
        	    "CollectionName" , 
	            "Name of the MCParticle collection"  ,
        	    _colName ,
           	 std::string("MCParticleSkimmed")
	    );
	    registerInputCollection( LCIO::VERTEX,
        	    "OutputCollectionName" , 
	            "Name of the Vertex collection"  ,
        	    _outputcolName ,
           	 std::string("MCVertex")
	    );
		_pdgs.push_back(BOTTOM_MESONS);
		_pdgs.push_back(CHARMED_MESONS);
		_pdgs.push_back(EXCEPTIONAL_MESONS);
	}	



	void TrashMCProcessor::init() 
	{ 
		streamlog_out(DEBUG) << "   init called  " << std::endl;
		printParameters() ;
		_nRun = 0 ;
		_nEvt = 0 ;
		_hfilename = "TrashMCTest.root";
		_hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;
		_hTree = new TTree( "Stats", "My test tree!" );
		_hTree->Branch("bdistance", &_bdistance, "bdistance/F");
		_hTree->Branch("bbardistance", &_bbardistance, "bbardistance/F");
		_hTree->Branch("firstVertexDistance", _firstVertexDistance, "firstVertexDistance[numberOfB0]/F");
		_hTree->Branch("secondVertexDistance", _secondVertexDistance, "secondVertexDistance[numberOfB0]/F");
		_hVertexTree = new TTree( "Vertexes", "My test tree!" );
		_hVertexTree->Branch("tag", &_tag, "tag/I");
		_hVertexTree->Branch("numberOfVertexes", &_numberOfVertexes, "numberOfVertexes/I");
		_hVertexTree->Branch("distance", _distanceFromIP, "distance[numberOfVertexes]/F");
		_hVertexTree->Branch("coordinates", _coordinates, "coordinates[numberOfVertexes][3]/F");
		_hVertexTree->Branch("PDG", _PDG, "PDG[numberOfVertexes]/I");
		_hVertexTree->Branch("generation", _generation, "generation[numberOfVertexes]/I");
		_hVertexTree->Branch("numberOfParticles", _numberOfParticles, "numberOfParticles[numberOfVertexes]/I");
		_hVertexTree->Branch("energyOfParticles", _energyOfParticles, "energyOfParticles[numberOfVertexes][15]/F");
		_hVertexTree->Branch("momentumOfParticles", _momentumOfParticles, "momentumOfParticles[numberOfVertexes][15]/F");
		_hVertexTree->Branch("massOfParticles", _massOfParticles, "massOfParticles[numberOfVertexes][15]/F");

		
	}


	void TrashMCProcessor::processRunHeader( LCRunHeader* run) 
	{ 
		_nRun++ ;
	} 
	
	void TrashMCProcessor::PrintParticle(MCParticle * particle)
	{
		if (!particle) 
		{
			return;
		}
		std::cout << std::fixed << std::setw( 6 ) << std::setprecision( 3 ) << std::setfill( ' ' );
		std::cout<<"|"<<particle->getPDG() <<"\t\t|"<<particle->getMass()<<"\t\t|"<<particle->getCharge()  <<"\t\t|"<<particle->getEnergy()<<"\t\t|"<<particle->getVertex()[0]<<"\t\t|"<<particle->getVertex()[1]<<"\t\t|"<<particle->getVertex()[2] <<"\t\t|\n";
	
	}
	void TrashMCProcessor::Write(MCOperator & opera, DecayChain * chain, int & number)
	{
		if (!chain || chain->GetSize() < 2) 
		{
			return;
		}
		for (int i = 1; i < chain->GetSize(); i++) 
		{
			_numberOfVertexes = number++;
			_PDG[_numberOfVertexes] = chain->GetParentPDG();
			_generation[_numberOfVertexes] = i+2;
			MCParticle * parent = chain->Get(i);
			if (!parent) 
			{
				break;
			}
			vector< MCParticle * > daughters = opera.ScanForVertexParticles(parent->getVertex(), 1e-3);
			if (daughters.size() < 1) 
			{
				std::cout<<"ERROR: NUMBER IS 0!\n";
			}
			_numberOfParticles[_numberOfVertexes] = daughters.size();
			for (int j = 0; j < daughters.size(); j++) 
			{
				MCParticle * daughter = daughters[j];
				if (abs(daughter->getCharge()) > 0.0 && daughter->isDecayedInCalorimeter()) 
				{
					_energyOfParticles[_numberOfVertexes][j] = daughter->getEnergy();
					_momentumOfParticles[_numberOfVertexes][j] = MathOperator::getModule( daughter->getMomentum());
					_massOfParticles[_numberOfVertexes][j] = daughter->getMass();
				}
			}
		}
	}
	void TrashMCProcessor::PrintChain(vector< MCParticle * > * chain)
	{
		if (!chain) 
		{
			return;
		}
		for (int i = 0; i < chain->size(); i++) 
		{
			PrintParticle(chain->at(i));
		}
	}

	void TrashMCProcessor::processEvent( LCEvent * evt ) 
	{ 
		try
		{
			LCCollection* col = evt->getCollection( _colName );
			std::cout<< "***********TrashMCProcessor*"<<_nEvt<<"***************\n";
			MCOperator opera(col);
			VertexMCOperator vertexOperator;
	 		_tag = opera.CheckProcessForPair(6);
			DecayChain * bChain = opera.Construct(string("b-quark decay chain"), 5, _pdgs);
			DecayChain * bbarChain = opera.Construct(string("bbar-quark decay chain"), -5, _pdgs);
			std::cout<<"\t|PDG\t\t|Mass\t\t|Charge\t\t|Energy\t\t|Vtx X\t\t|Vtx Y\t\t|Vtx Z\t\t|\n";
			for (int i = 0; i < bChain->GetSize(); i++) 
			{
				PrintParticle(bChain->Get(i));
			}
			for (int i = 0; i < bbarChain->GetSize(); i++) 
			{
				PrintParticle(bbarChain->Get(i));
			}
			vector< Vertex * > * bvetrexes = vertexOperator.Construct(bChain);
			vector< Vertex * > * bbarvetrexes = vertexOperator.Construct(bbarChain);
			if (bvetrexes) 
			{
				for (int i = 0; i < bvetrexes->size(); i++)
				{
				        std::cout<<"Vertex b-quark"<< i <<": " << bvetrexes->at(i)->getParameters()[0] << '\n';
					if (i == 0 && bvetrexes->at(i)) 
					{
						_firstVertexDistance[0] = bvetrexes->at(i)->getParameters()[0];
					}
					if (i == 1 && bvetrexes->at(i))
					{
					        _secondVertexDistance[0] = bvetrexes->at(i)->getParameters()[0];
						_bdistance = MathOperator::getDistance(bvetrexes->at(1)->getPosition(), bvetrexes->at(0)->getPosition());
					}

				}

			}
			if (bbarvetrexes) 
			{
				for (int i = 0; i < bbarvetrexes->size(); i++)
				{
				        std::cout<<"Vertex bbar-quark"<< i <<": " << bbarvetrexes->at(i)->getParameters()[0] << '\n';
					if (i == 0 && bbarvetrexes->at(i)) 
					{
						_firstVertexDistance[1] = bbarvetrexes->at(i)->getParameters()[0];

					}
					if (i == 1 && bbarvetrexes->at(i)) 
					{
						_secondVertexDistance[1] = bbarvetrexes->at(i)->getParameters()[0];
						_bbardistance = MathOperator::getDistance(bbarvetrexes->at(1)->getPosition(), bbarvetrexes->at(0)->getPosition());
					}

				}

			}
			if (bvetrexes && bbarvetrexes && bvetrexes->size() == 2 && bbarvetrexes->size() == 2 ) 
			{
				std::cout << "Everything went fine!\n";
			}
			WriteVertexCollection(evt, bvetrexes, bbarvetrexes);
			int number = 0;
			Write(opera,bChain,number);
			Write(opera,bbarChain,number);
			_numberOfVertexes = number;

			//std::cout<< "There was " << _numberOfB0 << " B-mesons.\n";
			_numberOfB0 = 2;
			_hTree->Fill();
			_hVertexTree->Fill();
			ClearVariables();
	
		}
		catch( DataNotAvailableException &e)
		{
			streamlog_out(DEBUG) << "No collection!" << std::endl ;
		}
		_nEvt ++ ;
	}

	void TrashMCProcessor::WriteVertexCollection(LCEvent * evt, vector< Vertex * > * bvertexes, vector< Vertex * > * bbarvertexes)
	{
		IMPL::LCCollectionVec * mc = new IMPL::LCCollectionVec ( EVENT::LCIO::VERTEX ) ;
		if (bvertexes) 
		{
			for (int i = 0; i < bvertexes->size(); i++) 
			{
				
				mc->addElement(bvertexes->at(i));
			}
		}
		if (bbarvertexes) 
		{
			for (int i = 0; i < bbarvertexes->size(); i++)
			{
			        mc->addElement(bbarvertexes->at(i));
			}
		}
		evt->addCollection( mc , _outputcolName ) ;
	}

	void TrashMCProcessor::ClearVariables()
	{
		_tag = false;
		_bdistance = -1.0;
		_bbardistance = -1.0;
		for (int i = 0; i < 2; i++) 
		{
			_firstVertexDistance[i] = 0.0;
			_secondVertexDistance[i] = 0.0;
			_numberOfB0 = 0;
		}
		for (int i = 0; i < MAXV; i++) 
		{
			_PDG[i] = 0;
			_generation[i] = 0;
			for (int j = 0; j < MAXV; j++) 
			{
				_energyOfParticles[i][j] = -1.0;
				_momentumOfParticles[i][j] = -1.0;
				_massOfParticles[i][j] = -1.0;
			}
		}
	}
	
	
	void TrashMCProcessor::check( LCEvent * evt ) 
	{ 
		// nothing to check here - could be used to fill checkplots in reconstruction processor
		
	}
	
	
	void TrashMCProcessor::end()
	{ 
		_hfile->cd();
		_hfile->Write();
		_hfile->Close();
	
	    //   std::cout << "TrashMCProcessor::end()  " << name() 
	    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    // 	    << std::endl ;
	
	}
} /* TTbarAnalysis */
