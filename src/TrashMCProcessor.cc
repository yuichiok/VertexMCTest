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
	    registerOutputCollection( LCIO::VERTEX,
        	    "OutputCollectionName" , 
	            "Name of the Vertex collection"  ,
        	    _outputcolName ,
           	 std::string("MCVertex")
	    );
	    registerOutputCollection( LCIO::MCPARTICLE,
        	    "OutputBStarName" , 
	            "Name of the Vertex collection"  ,
        	    _outputBStarName ,
           	 std::string("BStar")
	    );
	    registerOutputCollection( LCIO::MCPARTICLE,
        	    "OutputBStarName" , 
	            "Name of the Vertex collection"  ,
        	    _outputBStarName ,
           	 std::string("BStar")
	    );
	    registerOutputCollection(LCIO::MCPARTICLE,
	    		"QuarkCollectionName" , 
	            "Name of the b-quark collection"  ,
        	    _outputquarkcolName,
           	 std::string("MCbquarks")
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
	    _writeBonlyParameter = 1;
	    registerProcessorParameter("writeBonly" , 
	            "b parameter"  ,
        	    _writeBonlyParameter,
           	 _writeBonlyParameter
	    );
		_pdgs.push_back(BOTTOM_HADRONS);
		_pdgs.push_back(CHARMED_MESONS);
		_pdgs.push_back(EXCEPTIONAL_MESONS);
                ip[0] = 0.0;
                ip[1] = 0.0;
                ip[2] = 0.0;
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
		_hTree->Branch("totalBcharge", &_totalBcharge, "totalBcharge/I");
		_hTree->Branch("ccharge", &_ccharge, "ccharge/I");
		_hTree->Branch("cbarcharge", &_cbarcharge, "cbarcharge/I");
		_hTree->Branch("bcharge", &_bcharge, "bcharge/I");
		_hTree->Branch("bbarcharge", &_bbarcharge, "bbarcharge/I");
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
		_hTree->Branch("btotalnumber", &_btotalnumber, "btotalnumber/I");
		_hTree->Branch("bbartotalnumber", &_bbartotalnumber, "bbartotalnumber/I");
		_hTree->Branch("bptmiss", &_bptmiss, "bptmiss/F");
		_hTree->Branch("bbarptmiss", &_bbarptmiss, "bbarptmiss/F");
		//_hTree->Branch("bnumber_f", &_bnumber_f, "bnumber_f/I");
		//_hTree->Branch("bbarnumber_f", &_bbarnumber_f, "bbarnumber_f/I");
		//_hTree->Branch("cnumber_f", &_cnumber_f, "cnumber_f/I");
		//_hTree->Branch("cbarnumber_f", &_cbarnumber_f, "cbarnumber_f/I");
		//_hTree->Branch("firstVertexDistance", _firstVertexDistance, "firstVertexDistance[numberOfB0]/F");
		//_hTree->Branch("secondVertexDistance", _secondVertexDistance, "secondVertexDistance[numberOfB0]/F");
		_hVertexTree = new TTree( "Vertices", "My test tree!" );
		_hVertexTree->Branch("numberOfVertices", &_numberOfVertexes, "numberOfVertexes/I");
		_hVertexTree->Branch("distance", _distanceFromIP, "distance[numberOfVertexes]/F");
		_hVertexTree->Branch("coordinates", _coordinates, "coordinates[numberOfVertexes][3]/F");
		_hVertexTree->Branch("PDG", _PDG, "PDG[numberOfVertexes]/I");
		_hVertexTree->Branch("charge", _charge, "charge[numberOfVertexes]/I");
		_hVertexTree->Branch("generation", _generation, "generation[numberOfVertexes]/I");
		_hVertexTree->Branch("numberOfParticles", _numberOfParticles, "numberOfParticles[numberOfVertexes]/I");
		_hVertexTree->Branch("energyOfParticles", _energyOfParticles, "energyOfParticles[numberOfVertexes][15]/F");
		_hVertexTree->Branch("momentumOfParticles", _momentumOfParticles, "momentumOfParticles[numberOfVertexes][15]/F");
		_hVertexTree->Branch("massOfParticles", _massOfParticles, "massOfParticles[numberOfVertexes][15]/F");
		//******************************************************************************************************//
		_hTrackTree = new TTree( "Tracks", "My test tree!" );
		_hTrackTree->Branch("bnumber", &_bnumber, "bnumber/I");
		_hTrackTree->Branch("bbarnumber", &_bbarnumber, "bbarnumber/I");
		_hTrackTree->Branch("cnumber", &_cnumber, "cnumber/I");
		_hTrackTree->Branch("cbarnumber", &_cbarnumber, "cbarnumber/I");

		_hTrackTree->Branch("bptrack", _bptrack, "bptrack[bnumber]/F");
		_hTrackTree->Branch("cptrack", _cptrack, "cptrack[cnumber]/F");
		_hTrackTree->Branch("bbarptrack", _bbarptrack, "bbarptrack[bbarnumber]/F");
		_hTrackTree->Branch("cbarptrack", _cbarptrack, "cbarptrack[cbarnumber]/F");

		_hTrackTree->Branch("betatrack", _betatrack, "betatrack[bnumber]/F");
		_hTrackTree->Branch("cetatrack", _cetatrack, "cetatrack[cnumber]/F");
		_hTrackTree->Branch("bbaretatrack", _bbaretatrack, "bbaretatrack[bbarnumber]/F");
		_hTrackTree->Branch("cbaretatrack", _cbaretatrack, "cbaretatrack[cbarnumber]/F");
		
		_hTrackTree->Branch("boffsettrack", _boffsettrack, "boffsettrack[bnumber]/F");
		_hTrackTree->Branch("coffsettrack", _coffsettrack, "coffsettrack[cnumber]/F");
		_hTrackTree->Branch("bbaroffsettrack", _bbaroffsettrack, "bbaroffsettrack[bbarnumber]/F");
		_hTrackTree->Branch("cbaroffsettrack", _cbaroffsettrack, "cbaroffsetttrack[cbarnumber]/F");

		_hBStarTree = new TTree( "BStar", "My test tree!" );
		_hBStarTree->Branch("bstarnumber",  &_bstarnumber, "bstarnumber/I");
		_hBStarTree->Branch("bstarmomentum",  _bstarmomentum, "bstarmomentum[bstarnumber]/F");
		_hBStarTree->Branch("bstaroffset",  _bstaroffset, "bstaroffset[bstarnumber]/F");
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
			_charge[_numberOfVertexes] = chain->Get(i-1)->getCharge();
			MCParticle * parent = chain->Get(i);
			if (!parent) 
			{
				break;
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
	void TrashMCProcessor::WriteQuarksCollection(LCEvent * evt, std::vector< MCParticle * > & quarks)
	{
		std::cout<< "HERE!!!!\n";
		IMPL::LCCollectionVec * mc = new IMPL::LCCollectionVec ( LCIO::MCPARTICLE ) ;
		if (quarks.size() == 2) 
		{
			//std::cout << "PDG: " << quarks[0]->getPDG() << '\n';
			//std::cout << "PDG: " << quarks[1]->getPDG() << '\n';
			mc->addElement(quarks[0]);
			mc->addElement(quarks[1]);
		}
		evt->addCollection( mc , _outputquarkcolName ) ;
	}

	void TrashMCProcessor::AddProngs( VertexMCOperator & vertexOperator, MCOperator & opera, DecayChain * chain, vector< Vertex * > * verticies)
	{
		if (!verticies || !chain) 
		{
			return;
		}
		int vsize = verticies->size();
		if (vsize > 1) 
		{
			vector< MCParticle * > bdaughters = opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG());
			vertexOperator.AddProngs(verticies->at(0), bdaughters);
			vector< MCParticle * > cdaughters = opera.SelectStableCloseDaughters(chain->Get(1));
			vertexOperator.AddProngs(verticies->at(1), cdaughters);
		}
		else 
		{
			vector< MCParticle * > bdaughters = opera.SelectStableCloseDaughters(chain->Get(0));
			vertexOperator.AddProngs(verticies->at(0), bdaughters);
		}
	}

	void TrashMCProcessor::ExtractStarParticles(LCEvent * evt, MCOperator opera,  DecayChain * bChainRaw, DecayChain * bChain, int v)
	{
		if (!bChain) 
		{
			return;
		}
		vector< MCParticle * > daughters = opera.SelectStableCloseDaughters(bChainRaw->Get(0), bChain->Get(0)->getPDG());
		if (daughters.size() < 1) 
		{
			//return;
		}
		IMPL::LCCollectionVec * mc = new IMPL::LCCollectionVec ( LCIO::MCPARTICLE ) ;
		for (unsigned int i = 0; i < daughters.size(); i++) 
		{
			mc->addElement(new MCParticleImpl((const IMPL::MCParticleImpl&)(*daughters[i])));
			PrintParticle(daughters[i]);
		}
		switch(v)
		{
		        case 1:
			evt->addCollection( mc , _outputBStarName ) ;
			break;		
		        case -1:
			evt->addCollection( mc , _outputBbarStarName ) ;
			break;		
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
			vector< MCParticle * > bquarks = opera.GetPairParticles(_pdgs[0]);
			_nEvt ++ ;
			
			std::cout<<"\t|PDG\t\t|Mass\t\t|Charge\t\t|Energy\t\t|Vtx X\t\t|Vtx Y\t\t|Vtx Z\t\t|\n";
			DecayChain * bChainRaw = opera.Construct(string("b-quark decay chain"), 5, _pdgs);
			DecayChain * bChain = opera.RefineDecayChain(bChainRaw, _pdgs);
			IMPL::LCCollectionVec * mc = new IMPL::LCCollectionVec ( LCIO::MCPARTICLE ) ;
			if (bChain) 
			{
				vector< MCParticle * > daughters = opera.SelectStableCloseDaughters(bChainRaw->Get(0), bChain->Get(0)->getPDG());	
				std::cout<<"Additional B particles: \n";
				for (int i = 0; i < daughters.size(); i++) 
				{
					vector< float > direction = MathOperator::getDirection(daughters[i]->getMomentum());
					_bstaroffset[_bstarnumber] = MathOperator::getDistanceTo(ip, direction, bChain->Get(1)->getVertex());
					_bstarmomentum[_bstarnumber++] = MathOperator::getModule(daughters[i]->getMomentum());
					mc->addElement(new MCParticleImpl((const IMPL::MCParticleImpl&)(*daughters[i])));
					PrintParticle(daughters[i]);
				}
			}
			DecayChain * bbarChainRaw = opera.Construct(string("bbar-quark decay chain"), -5, _pdgs);
			DecayChain * bbarChain = opera.RefineDecayChain(bbarChainRaw, _pdgs);
			if (bbarChain) 
			{
				vector< MCParticle * > daughters = opera.SelectStableCloseDaughters(bbarChainRaw->Get(0), bbarChain->Get(0)->getPDG());	
				std::cout<<"Additional Bbar particles: \n";
				for (int i = 0; i < daughters.size(); i++) 
				{
					vector< float > direction = MathOperator::getDirection(daughters[i]->getMomentum());
					_bstaroffset[_bstarnumber] = MathOperator::getDistanceTo(ip, direction, bbarChain->Get(1)->getVertex());
					_bstarmomentum[_bstarnumber++] = MathOperator::getModule(daughters[i]->getMomentum());
					mc->addElement(new MCParticleImpl((const IMPL::MCParticleImpl&)(*daughters[i])));
					PrintParticle(daughters[i]);
				}
			}
			evt->addCollection( mc , _outputBStarName ) ;

			vector< Vertex * > * bverticies = vertexOperator.Construct(bChain);
			vector< Vertex * > * bbarverticies = vertexOperator.Construct(bbarChain);
			
			AddProngs(vertexOperator, opera, bChain, bverticies);
			AddProngs(vertexOperator, opera, bbarChain, bbarverticies);

			Write(opera, bChain,bverticies);
			Write(opera, bbarChain,bbarverticies);
			
			std::cout<<"Total b number: " << _btotalnumber << '\n';
			std::cout<<"Total bbar number: " << _bbartotalnumber << '\n';
			
			WriteQuarksCollection(evt, bquarks);
			WriteVertexCollection(evt, bverticies, bbarverticies);
			
			int number = 0;
			Write(opera,bbarChain,number);
			Write(opera,bChain,number);
			_numberOfVertexes = number;

			//std::cout<< "There was " << _numberOfB0 << " B-mesons.\n";
			_numberOfB0 = 2;
			_hTree->Fill();
			_hBStarTree->Fill();
			_bnumber = (_bnumber < 0)? 0: _bnumber;
			_bbarnumber = (_bbarnumber < 0)? 0: _bbarnumber;
			_cnumber = (_cnumber < 0)? 0: _cnumber;
			_cbarnumber = (_cbarnumber < 0)? 0: _cbarnumber;
			_hTrackTree->Fill();
			_hVertexTree->Fill();
			ClearVariables();
	
		}
		catch( DataNotAvailableException &e)
		{
			streamlog_out(DEBUG) << "No collection!" << std::endl ;
		}
	}
	void TrashMCProcessor::Write(MCOperator & opera, DecayChain * chain, vector< Vertex * > * verticies)
	{
		if (!chain || !chain->Get(0)) 
		{
			return;
		}
		for (int i = 0; i < chain->GetSize(); i++) 
		{
			PrintParticle(chain->Get(i));
		}
		if (chain->GetParentPDG() > 0 && chain->GetSize() > 1) 
		{
			vector< MCParticle * > daughters = opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG()); //opera.ScanForVertexParticles(bverticies->at(0)->getPosition(), 1e-3);
			_bcharge = (int)chain->Get(0)->getCharge();
			_bnumber = daughters.size();
			_bnumber_f = opera.CheckDaughterVisibility(daughters).size();
			_bIPdistance = verticies->at(0)->getParameters()[0];
			Write(daughters, 1);
			std::cout<<"Vertex b-quark: " << verticies->at(0)->getParameters()[0]<< " n-tracks: " << _bnumber << '\n';
			vector< MCParticle * > cdaughters =opera.SelectStableCloseDaughters(chain->Get(1)); //opera.ScanForVertexParticles(bverticies->at(1)->getPosition(), 1e-3);
			_cnumber = cdaughters.size();
			_cnumber_f = opera.CheckDaughterVisibility(cdaughters).size();
			Write(cdaughters, 2);
			std::cout<<"Vertex c-quark: " << verticies->at(1)->getParameters()[0] <<" n-tracks: " << _cnumber <<  '\n';
			_bdistance = MathOperator::getDistance(verticies->at(1)->getPosition(), verticies->at(0)->getPosition());
			_btotalnumber = _cnumber + _bnumber;
			std::cout<<"Checking b-quark meson...\n";
			bool compatible = opera.CheckCompatibility(daughters, chain->Get(0), chain->Get(1)->getCharge());
				
			std::cout<<"Checking c-quark meson...\n";
			compatible = opera.CheckCompatibility(cdaughters, chain->Get(1));
			_bptmiss = getMissingPt(daughters, cdaughters, verticies->at(0));
			std::cout<<"Missing pt for b-quark hadron: " << _bptmiss << "\n";	
			_ccharge = (int)chain->Get(1)->getCharge();
			_bmomentum = MathOperator::getModule(chain->Get(0)->getMomentum());
			_baccuracy = opera.GetAccuracy(chain->Get(0), _aParameter, _bParameter); 
			_cmomentum = MathOperator::getModule(chain->Get(1)->getMomentum());
			_caccuracy = opera.GetAccuracy(chain->Get(1), _aParameter, _bParameter);
		}
		if (chain->GetParentPDG() < 0 && chain->GetSize() > 1) 
		{
			vector< MCParticle * > daughters =opera.SelectStableCloseDaughters(chain->Get(0), chain->Get(1)->getPDG()); // opera.ScanForVertexParticles(bbarverticies->at(0)->getPosition(), 1e-3);
			Write(daughters, -1);
			_bbarcharge = (int) chain->Get(0)->getCharge();
			_bbarnumber = daughters.size();
			_bbarIPdistance = verticies->at(0)->getParameters()[0];
			_bbarnumber_f = opera.CheckDaughterVisibility(daughters).size();
		        std::cout<<"Vertex bbar-quark"<< 0 <<": " << verticies->at(0)->getParameters()[0] <<" n-tracks: " << _bbarnumber <<  '\n';
			vector< MCParticle * > cdaughters= opera.SelectStableCloseDaughters(chain->Get(1)); //opera.ScanForVertexParticles(bbarverticies->at(1)->getPosition(), 1e-3);
			_cbarnumber = cdaughters.size();
			std::cout<<"Vertex cbar-quark: " << verticies->at(1)->getParameters()[0] <<" n-tracks: " <<_cbarnumber << '\n';
			_cbarnumber_f = opera.CheckDaughterVisibility(cdaughters).size();
			_bbardistance = MathOperator::getDistance(verticies->at(1)->getPosition(), verticies->at(0)->getPosition());
			_bbartotalnumber = _cbarnumber + _bbarnumber;
			_cbarcharge = (int) chain->Get(1)->getCharge();
			Write(cdaughters, -2);
			_bbarptmiss = getMissingPt(daughters, cdaughters, verticies->at(0));
			std::cout<<"Missing pt for bbar-quark hadron: " << _bbarptmiss << "\n";	
			_bbarmomentum = MathOperator::getModule(chain->Get(0)->getMomentum());
			_bbaraccuracy = opera.GetAccuracy(chain->Get(0), _aParameter, _bParameter); 
			_cbarmomentum = MathOperator::getModule(chain->Get(1)->getMomentum());
			_cbaraccuracy = opera.GetAccuracy(chain->Get(1), _aParameter, _bParameter); 

		}
	}
	double TrashMCProcessor::getMissingPt(vector< MCParticle * > & bdaugthers, vector< MCParticle * > & cdaughters, Vertex * vertex)
	{
		const float * position = vertex->getPosition();
		for (int i = 0; i < 3; i++) 
		{
			std::cout << i << ": " << position[i];
		}
		std::cout << '\n';
		vector< const double * > vectors;
		for (int i = 0; i < bdaugthers.size(); i++) 
		{
			vectors.push_back(bdaugthers[i]->getMomentum());
		}
		for (int i = 0; i < cdaughters.size(); i++) 
		{
			vectors.push_back(cdaughters[i]->getMomentum());
		}
		double missing =  MathOperator::getMissingPt(vectors, position);
		return missing;
	}
	void TrashMCProcessor::Write(vector< MCParticle * > daughters, int v)
	{
		float * offset = NULL;
		float * pt = NULL;
		float * p = NULL;
		float * eta = NULL;
		switch(v)
		{
			case 1:
			offset = _boffsettrack;
			p = _bptrack;
			eta = _betatrack;
			break;
			case 2:
			offset = _coffsettrack;
			p = _cptrack;
			eta = _cetatrack;
			break;
			case -2:
			offset = _cbaroffsettrack;
			p = _cbarptrack;
			eta = _cbaretatrack;
			break;
			case -1:
			offset = _bbaroffsettrack;
			p = _bbarptrack;
			eta = _bbaretatrack;
			break;

		}
		int size = daughters.size();
		for (int i = 0; i < size; i++) 
		{
			MCParticle * daughter = daughters[i];
			vector< float > direction = MathOperator::getDirection(daughter->getMomentum());
			offset[i] = MathOperator::getDistanceTo(ip, direction, daughter->getVertex());
			//float pt[i] = MathOperator::getPt(daughter->getMomentum());
			p[i] = MathOperator::getModule(daughter->getMomentum());
			eta[i] = MathOperator::getAngles(direction)[1];
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
		_bstarnumber = 0;
	  	_bptmiss = -1.0;
	  	_bbarptmiss = -1.0;
		_totalBcharge = -3.0;
		_ccharge = -3.0;
		_cbarcharge = -3.0;
		_bcharge = -3.0;
		_bbarcharge = -3.0;
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
		_bbartotalnumber = -1;
		_btotalnumber = -1;
		_bnumber_f = -1;
		_bbarnumber_f = -1;
		_cnumber_f = -1;
		_cbarnumber_f = -1;
		for (int i = 0; i < 2; i++) 
		{
			_firstVertexDistance[i] = 0.0;
			_secondVertexDistance[i] = 0.0;
			_numberOfB0 = 0;
		}
		for (int i = 0; i < MAXV; i++) 
		{
			_bstaroffset[i]  = -1.0;
			_bstarmomentum[i] = -1.0;
			_bptrack[i] = -1.0;
			_betatrack[i] = -1.0;
			_boffsettrack[i] = -1.0;
			_bbarptrack[i] = -1.0;
			_bbaretatrack[i] = -1.0;
			_bbaroffsettrack[i] = -1.0;
			_cptrack[i] = -1.0;
			_cetatrack[i] = -1.0;
			_coffsettrack[i] = -1.0;
			_cbarptrack[i] = -1.0;
			_cbaretatrack[i] = -1.0;
			_cbaroffsettrack[i]  = -1.0;

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
