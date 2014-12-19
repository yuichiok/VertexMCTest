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
	    _tagParameter = 6;
	    registerProcessorParameter("tagPDG" , 
	            "PDG of desired particle"  ,
        	    _tagParameter,
           	 _tagParameter
	    );
	    _aParameter = 0.005;
	    registerProcessorParameter("a" , 
	            "a parameter of accuracy in mm"  ,
        	    _aParameter,
           	 _aParameter
	    );
	    _bParameter = 0.01;
	    registerProcessorParameter("b" , 
	            "b parameter of accuracy in mm"  ,
        	    _bParameter,
           	 _bParameter
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
		_hTree->Branch("tag", &_tag, "tag/I");
		_hTree->Branch("baccuracy", &_baccuracy, "baccuracy/F");
		_hTree->Branch("bbaraccuracy", &_bbaraccuracy, "bbaraccuracy/F");
		_hTree->Branch("bIPdistance", &_bIPdistance, "bIPdistance/F");
		_hTree->Branch("bbarIPdistance", &_bbarIPdistance, "bbarIPdistance/F");
		_hTree->Branch("bdistance", &_bdistance, "bdistance/F");
		_hTree->Branch("bbardistance", &_bbardistance, "bbardistance/F");
		_hTree->Branch("bmomentum", &_bmomentum, "bmomentum/F");
		_hTree->Branch("bbarmomentum", &_bbarmomentum, "bbarmomentum/F");
		_hTree->Branch("cmomentum", &_cmomentum, "cmomentum/F");
		_hTree->Branch("cbarmomentum", &_cbarmomentum, "cbarmomentum/F");
		_hTree->Branch("caccuracy", &_caccuracy, "caccuracy/F");
		_hTree->Branch("cbaraccuracy", &_cbaraccuracy, "cbaraccuracy/F");
		_hTree->Branch("bnumber", &_bnumber, "bnumber/I");
		_hTree->Branch("bbarnumber", &_bbarnumber, "bbarnumber/I");
		_hTree->Branch("cnumber", &_cnumber, "cnumber/I");
		_hTree->Branch("cbarnumber", &_cbarnumber, "cbarnumber/I");
		//_hTree->Branch("firstVertexDistance", _firstVertexDistance, "firstVertexDistance[numberOfB0]/F");
		//_hTree->Branch("secondVertexDistance", _secondVertexDistance, "secondVertexDistance[numberOfB0]/F");
		_hVertexTree = new TTree( "Vertexes", "My test tree!" );
		_hVertexTree->Branch("numberOfVertexes", &_numberOfVertexes, "numberOfVertexes/I");
		_hVertexTree->Branch("distance", _distanceFromIP, "distance[numberOfVertexes]/F");
		_hVertexTree->Branch("coordinates", _coordinates, "coordinates[numberOfVertexes][3]/F");
		_hVertexTree->Branch("PDG", _PDG, "PDG[numberOfVertexes]/I");
		_hVertexTree->Branch("charge", _charge, "charge[numberOfVertexes]/I");
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
			_generation[_numberOfVertexes] = i+1;
			MCParticle * parent = chain->Get(i);
			if (!parent) 
			{
				break;
			}
			/*vector< MCParticle * > daughters = opera.ScanForVertexParticles(parent->getVertex(), 1e-3);
			if (daughters.size() < 1) 
			{
				std::cout<<"ERROR: NUMBER IS 0!\n";
			}
			_numberOfParticles[_numberOfVertexes] = daughters.size();
			int charge = 0;
			for (int j = 0; j < daughters.size(); j++) 
			{
				MCParticle * daughter = daughters[j];
				charge += daughter->getCharge();
				_energyOfParticles[_numberOfVertexes][j] = daughter->getEnergy();
				_momentumOfParticles[_numberOfVertexes][j] = MathOperator::getModule( daughter->getMomentum());
				_massOfParticles[_numberOfVertexes][j] = daughter->getMass();
			}
			_charge[_numberOfVertexes] = chain->Get(i-1)->getCharge();*/
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
	 		_tag = opera.CheckProcessForPair(_tagParameter);
			_nEvt ++ ;
			if (!_tag) 
			{
				return;
			}
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
			vector< Vertex * > * bverticies = vertexOperator.Construct(bChain);
			vector< Vertex * > * bbarverticies = vertexOperator.Construct(bbarChain);
			if (bverticies && bverticies->size() > 1) 
			{
			        std::cout<<"Vertex b-quark: " << bverticies->at(0)->getParameters()[0] << '\n';
				vector< MCParticle * > daughters = opera.SelectStableCloseDaughters(bChain->Get(0)); //opera.ScanForVertexParticles(bverticies->at(0)->getPosition(), 1e-3);
				bool compatible = opera.CheckCompatibility(daughters, bChain->Get(0),(int)bChain->Get(1)->getCharge());
				if (!compatible) 
				{
					for (unsigned int i = 0; i < daughters.size(); i++) 
					{
						PrintParticle(daughters[i]);
					}
				}
				_bnumber = daughters.size();
				std::cout<<"Vertex c-quark: " << bverticies->at(1)->getParameters()[0] << '\n';
				vector< MCParticle * > cdaughters =opera.SelectStableCloseDaughters(bChain->Get(1), false); //opera.ScanForVertexParticles(bverticies->at(1)->getPosition(), 1e-3);
				compatible = opera.CheckCompatibility(cdaughters, bChain->Get(1),0);
				if (!compatible) 
				{
					for (unsigned int i = 0; i < cdaughters.size(); i++) 
					{
						PrintParticle(cdaughters[i]);
					}
				}
				_cnumber = cdaughters.size();

				_bIPdistance = bverticies->at(0)->getParameters()[0];
			        //_secondVertexDistance[0] = bverticies->at(i)->getParameters()[0];
				_bdistance = MathOperator::getDistance(bverticies->at(1)->getPosition(), bverticies->at(0)->getPosition());
			}
			if (bbarverticies && bbarverticies->size() > 1) 
			{
			        std::cout<<"Vertex bbar-quark"<< 0 <<": " << bbarverticies->at(0)->getParameters()[0] << '\n';
				vector< MCParticle * > daughters =opera.SelectStableCloseDaughters(bbarChain->Get(0)); // opera.ScanForVertexParticles(bbarverticies->at(0)->getPosition(), 1e-3);
				_bbarnumber = daughters.size();
				std::cout<<"Vertex cbar-quark: " << bbarverticies->at(1)->getParameters()[0] << '\n';
				vector< MCParticle * > cdaughters = opera.SelectStableCloseDaughters(bbarChain->Get(1), false); //opera.ScanForVertexParticles(bbarverticies->at(1)->getPosition(), 1e-3);
				_cbarnumber = cdaughters.size();
				_bbarIPdistance = bbarverticies->at(0)->getParameters()[0];
				//_secondVertexDistance[1] = bbarverticies->at(i)->getParameters()[0];
				_bbardistance = MathOperator::getDistance(bbarverticies->at(1)->getPosition(), bbarverticies->at(0)->getPosition());
			}
			if (bChain && bChain->GetSize() > 1) 
			{
				_bmomentum = MathOperator::getModule(bChain->Get(0)->getMomentum());
				_cmomentum = MathOperator::getModule(bChain->Get(1)->getMomentum());
				_baccuracy = opera.GetAccuracy(bChain->Get(0), _aParameter, _bParameter); 
				_caccuracy = opera.GetAccuracy(bChain->Get(1), _aParameter, _bParameter);
			}
			if (bbarChain && bbarChain->GetSize() > 1) 
			{
				_bbarmomentum = MathOperator::getModule(bbarChain->Get(0)->getMomentum());
				_bbaraccuracy = opera.GetAccuracy(bbarChain->Get(0), _aParameter, _bParameter); 
				_cbarmomentum = MathOperator::getModule(bbarChain->Get(1)->getMomentum());
				_cbaraccuracy = opera.GetAccuracy(bbarChain->Get(1), _aParameter, _bParameter); 
			}
			WriteVertexCollection(evt, bverticies, bbarverticies);
			int number = 0;
			Write(opera,bbarChain,number);
			Write(opera,bChain,number);
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
		_bmomentum = -1.0;
		_bbarmomentum = -1.0;
		_bbaraccuracy = -1.0;
		_baccuracy = -1.0;
		_cmomentum = -1.0;
		_cbarmomentum = -1.0;
		_cbaraccuracy = -1.0;
		_caccuracy = -1.0;
		_bnumber = -1;
		_bbarnumber = -1;
		_cnumber = -1;
		_cbarnumber = -1;
		for (int i = 0; i < 2; i++) 
		{
			_firstVertexDistance[i] = 0.0;
			_secondVertexDistance[i] = 0.0;
			_numberOfB0 = 0;
		}
		for (int i = 0; i < MAXV; i++) 
		{
			_charge[i] = 0;
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
