#include <array>
#include <map>
#include <set>
#include <tuple>
#define NO_SEC_RETURN_TRUE if(this->irreducible_sector_.empty()) return true;
#define NO_QUADS_RETURN_TRUE if(this->irreducible_quads_.empty()) return true;
#include "../global/Array_Operator.h"
namespace RI
{
	using namespace Array_Operator;
	template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
	class Symmetry_Filter
	{
		using TC = std::array<Tcell, Ndim>;
		using TAC = std::pair<TA, TC>;

		using TIJ = std::pair<TA, TA>;
		using TIJR = std::pair<TIJ, TC>;
		using Tsec = std::map<TIJ, std::set<TC>>;
		using Tquad_IJR = std::tuple<TIJR, TC>;
		using Tquad = std::map<TIJR, std::set<Tquad_IJR>>;
	public:
		Symmetry_Filter(const TC& period_in, const Tsec& irsec, const Tquad& irquad = {})
			:period(period_in), irreducible_sector_(irsec), irreducible_quads_(irquad) {}
		bool in_irreducible_sector(const TA& Aa, const TAC& Ab) const
		{
			NO_SEC_RETURN_TRUE;
			const TIJ& ap = { Aa, Ab.first };
			if (irreducible_sector_.find(ap) != irreducible_sector_.end())
				if (irreducible_sector_.at(ap).find(Ab.second % this->period) != irreducible_sector_.at(ap).end())
					return true;
			return false;
		}
		bool in_irreducible_sector(const TAC& Aa, const TAC& Ab) const
		{
			NO_SEC_RETURN_TRUE;
			const TC dR = (Ab.second - Aa.second) % this->period;
			const std::pair<TA, TA> ap = { Aa.first, Ab.first };
			if (irreducible_sector_.find(ap) != irreducible_sector_.end())
				if (irreducible_sector_.at(ap).find(dR) != irreducible_sector_.at(ap).end())
					return true;
			return false;
		}
		bool is_I_in_irreducible_sector(const TA& Aa) const
		{
			NO_SEC_RETURN_TRUE;
			for (const auto& apRs : irreducible_sector_)
				if (apRs.first.first == Aa)return true;
			return false;
		}
		bool is_J_in_irreducible_sector(const TA& Ab) const
		{
			NO_SEC_RETURN_TRUE;
			for (const auto& apRs : irreducible_sector_)
				if (apRs.first.second == Ab)return true;
			return false;
		}
		bool in_irreducible_quad(const TIJR& IJRij, const TIJR& KLRkl, const TC& Rik) const
		{
			NO_SEC_RETURN_TRUE;
			NO_QUADS_RETURN_TRUE;
			auto find_IJRij = irreducible_quads_.find(IJRij);
			if (find_IJRij == irreducible_quads_.end()) return false;
			return find_IJRij->find({ KLRkl, Rik }) != find_IJRij->end();
		}
		bool is_IKL_in_irreducible_quad(const TA& I, const TIJR& KLRkl, const TC& Rik)
		{
			NO_SEC_RETURN_TRUE;
			NO_QUADS_RETURN_TRUE;
			for (auto& apRs : irreducible_sector_)
				if (apRs.first.first == I)
					for (auto& Rij : apRs.second)
						if (in_irreducible_quad({ apRs.first, Rij }, KLRkl, Rik))
							return true;
			return false;
		}
		bool is_JKL_in_irreducible_quad(const TA& J, const TIJR& KLRkl, const TC& Rjk)
		{
			NO_SEC_RETURN_TRUE;
			NO_QUADS_RETURN_TRUE;
			for (auto& apRs : irreducible_sector_)
				if (apRs.first.second == J)
					for (auto& Rij : apRs.second)
						if (in_irreducible_quad({ apRs.first, Rij }, KLRkl, (Rij + Rjk) % period))
							return true;
			return false;
		}
		TIJR get_IJR(const TA& I, const TAC& J) const
		{
			return { {I,J.first}, J.second % this->period };
		}
		TIJR get_IJR(const TAC& I, const TAC& J) const
		{
			return { {I.first,J.first}, (J.second - I.second) % this->period };
		}

		bool is_IKL_in_irreducible_quad(const TA& I, const TAC& K, const TAC& L)
		{
			return is_IKL_in_irreducible_quad(I, get_IJR(K, L), K.second % this->period);
		}
		bool is_IKL_in_irreducible_quad(const TAC& I, const TA& K, const TAC& L)
		{
			return is_IKL_in_irreducible_quad(I, get_IJR(K, L), (-I.second) % this->period);
		}
		bool is_JKL_in_irreducible_quad(const TA& J, const TAC& K, const TAC& L)
		{
			return is_JKL_in_irreducible_quad(J, get_IJR(K, L), K.second % this->period);
		}
		bool is_JKL_in_irreducible_quad(const TA& J, const TA& K, const TAC& L)
		{
			return is_JKL_in_irreducible_quad(J, get_IJR(K, L), (-J.second) % this->period);
		}
	private:
		const Tsec& irreducible_sector_;
		const Tquad& irreducible_quads_;
		const TC& period;
	};

}