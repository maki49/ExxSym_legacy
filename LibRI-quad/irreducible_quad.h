#include "irreducible_sector.h"

namespace ModuleSymmetry
{
    // using Tquad_IJR = std::tuple<TapR, TC, TC>;    // {KLR, R_IK, R_JL}
    using Tquad_IJR = std::pair<TapR, TC>;    // {KLR, R_IK}
    using Tquad = std::pair<TapR, Tquad_IJR>;   // {irIJR, {KLR, R_IK, R_JL}}
    class Irreducible_Quad : public Irreducible_Sector
    {
    public:
        void find_irsector_invariant_operations(const Symmetry& symm);
        void find_irreducible_quads(const Symmetry& symm, const Atom* atoms, const Statistics& st,
            const std::vector<TC>& Rs, const TC& period, const Lattice& lat);
        void print_quads_stars()const;
        void print_irreducible_quads()const;
    private:
        std::map<TapR, std::set<int>> irsector_invariant_ops_;
        std::map<TapR, std::vector<std::map<int, Tquad_IJR>>> quads_stars_;   // irIJR, {isym, {KLR, R_IK, R_JL}}
        std::map<TapR, std::set<Tquad_IJR>> irreducible_quads_;
    };
}
