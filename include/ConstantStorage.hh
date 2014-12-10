#include <iostream>
#include <vector>
#ifndef _TTbarAnalysis_h
#define _TTbarAnalysis_h 1
namespace TTbarAnalysis 
{
	enum MESONS
	{
		EXCEPTIONAL_MESONS,
		BOTTOM_MESONS,
		CHARMED_MESONS,
		STRANGE_MESONS,
		_max_MESONS = STRANGE_MESONS
	};
	class ConstantStorage 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			ConstantStorage ( ){}
			virtual ~ConstantStorage (){}
		//
		//	Methods
		//
		static const std::vector<int> & CHARMED_MESONS_PGD();
		static const std::vector<int> & STRANGE_MESONS_PGD();
		static const std::vector<int> & BOTTOM_MESONS_PGD();
		
		static const std::vector<int> & GET_PDG(MESONS type);
		//static constructor pattern
		static class _init
		{
			public:
			_init()
			{
				///////////////CHARMED/////////////////////
				myCharmedMesonsPDGs.push_back(433);
				myCharmedMesonsPDGs.push_back(413);
				myCharmedMesonsPDGs.push_back(431);
				myCharmedMesonsPDGs.push_back(411);
				myCharmedMesonsPDGs.push_back(421);
				myCharmedMesonsPDGs.push_back(10411);
				myCharmedMesonsPDGs.push_back(10421);
				myCharmedMesonsPDGs.push_back(423);
				myCharmedMesonsPDGs.push_back(10413);
				myCharmedMesonsPDGs.push_back(10423);
				myCharmedMesonsPDGs.push_back(20413);
				myCharmedMesonsPDGs.push_back(20423);
				myCharmedMesonsPDGs.push_back(415);
				myCharmedMesonsPDGs.push_back(425);
				myCharmedMesonsPDGs.push_back(10431);
				myCharmedMesonsPDGs.push_back(10433);
				myCharmedMesonsPDGs.push_back(20433);
				myCharmedMesonsPDGs.push_back(435);
				///////////////BOTTOM//////////////////////
				myBottomMesonsPDGs.push_back(511);
				myBottomMesonsPDGs.push_back(521);
				myBottomMesonsPDGs.push_back(513);
				myBottomMesonsPDGs.push_back(523);
				myBottomMesonsPDGs.push_back(515);
				myBottomMesonsPDGs.push_back(525);
				myBottomMesonsPDGs.push_back(531);
				myBottomMesonsPDGs.push_back(533);
				myBottomMesonsPDGs.push_back(535);
				myBottomMesonsPDGs.push_back(541);
				myBottomMesonsPDGs.push_back(543);
				myBottomMesonsPDGs.push_back(545);
				myBottomMesonsPDGs.push_back(10511);
				myBottomMesonsPDGs.push_back(10521);
				myBottomMesonsPDGs.push_back(10513);
				myBottomMesonsPDGs.push_back(20513);
				myBottomMesonsPDGs.push_back(20523);
				myBottomMesonsPDGs.push_back(10531);
				myBottomMesonsPDGs.push_back(10533);
				myBottomMesonsPDGs.push_back(20533);
				myBottomMesonsPDGs.push_back(10541);
				myBottomMesonsPDGs.push_back(10543);
				myBottomMesonsPDGs.push_back(20543);
				///////////////STRANGE//////////////////
				myStrangeMesonsPDGs.push_back(130);
				myStrangeMesonsPDGs.push_back(310);
				myStrangeMesonsPDGs.push_back(311);
				myStrangeMesonsPDGs.push_back(321);
				myStrangeMesonsPDGs.push_back(10311);
				myStrangeMesonsPDGs.push_back(10321);
				myStrangeMesonsPDGs.push_back(100311);
				myStrangeMesonsPDGs.push_back(100321);
				myStrangeMesonsPDGs.push_back(200311);
				myStrangeMesonsPDGs.push_back(200321);
				myStrangeMesonsPDGs.push_back(9000311);
				myStrangeMesonsPDGs.push_back(9000321);
				myStrangeMesonsPDGs.push_back(313);
				myStrangeMesonsPDGs.push_back(323);
				myStrangeMesonsPDGs.push_back(10313);
				myStrangeMesonsPDGs.push_back(10323);
				myStrangeMesonsPDGs.push_back(20313);
				myStrangeMesonsPDGs.push_back(20323);
				myStrangeMesonsPDGs.push_back(100313);
				myStrangeMesonsPDGs.push_back(10323);
				myStrangeMesonsPDGs.push_back(9000313);
				myStrangeMesonsPDGs.push_back(9000323);
				myStrangeMesonsPDGs.push_back(30313);
				myStrangeMesonsPDGs.push_back(30323);
				myStrangeMesonsPDGs.push_back(315);
				myStrangeMesonsPDGs.push_back(325);
				myStrangeMesonsPDGs.push_back(9000315);
				myStrangeMesonsPDGs.push_back(9000325);
				myStrangeMesonsPDGs.push_back(10315);
				myStrangeMesonsPDGs.push_back(10325);
				myStrangeMesonsPDGs.push_back(20315);
				myStrangeMesonsPDGs.push_back(20325);
				myStrangeMesonsPDGs.push_back(100315);
				myStrangeMesonsPDGs.push_back(100325);
				myStrangeMesonsPDGs.push_back(9010315);
				myStrangeMesonsPDGs.push_back(9010325);
				myStrangeMesonsPDGs.push_back(317);
				myStrangeMesonsPDGs.push_back(327);
				myStrangeMesonsPDGs.push_back(9010317);
				myStrangeMesonsPDGs.push_back(9010327);
				myStrangeMesonsPDGs.push_back(319);
				myStrangeMesonsPDGs.push_back(329);
				myStrangeMesonsPDGs.push_back(9000319);
				myStrangeMesonsPDGs.push_back(9000329);
				
			}
		} _initializer;
		
		private:
		//
		//	Data
		//
			static std::vector<int> myCharmedMesonsPDGs; 
			static std::vector<int> myBottomMesonsPDGs; 
			static std::vector<int> myStrangeMesonsPDGs; 
		//
		//	Private methods
		//


	};
} /* TTbarAnalysis */
#endif
