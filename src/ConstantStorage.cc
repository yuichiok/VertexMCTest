#include "ConstantStorage.hh"
using std::vector;
namespace TTbarAnalysis 
{
	vector<int> ConstantStorage::myCharmedMesonsPDGs;
	vector<int> ConstantStorage::myBottomMesonsPDGs;
	vector<int> ConstantStorage::myStrangeMesonsPDGs;

	ConstantStorage::_init ConstantStorage::_initializer;

	const std::vector<int> & ConstantStorage::CHARMED_MESONS_PGD()
	{
		return myCharmedMesonsPDGs;
	}
	const std::vector<int> & ConstantStorage::STRANGE_MESONS_PGD()
	{
		return myStrangeMesonsPDGs;
	}
	const std::vector<int> & ConstantStorage::BOTTOM_MESONS_PGD()
	{
		return myBottomMesonsPDGs;
	}
	const std::vector<int> & ConstantStorage::GET_PDG(MESONS type)
	{
		switch(type)
		{
			case BOTTOM_MESONS:
			return BOTTOM_MESONS_PGD();
			case CHARMED_MESONS:
			return CHARMED_MESONS_PGD();
			case STRANGE_MESONS:
			return STRANGE_MESONS_PGD();
			default:
			return vector<int>();
		}
	}

} /* TTbarAnalysis */
