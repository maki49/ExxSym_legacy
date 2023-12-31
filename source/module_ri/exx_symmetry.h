#pragma once
#include <vector>
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_base/vector3.h"
#include "module_psi/psi.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_cell/module_symmetry/symmetry.h"
namespace ExxSym
{
    /// @brief Rearrange the $\nu$ index of a matrix $S_{\mu,\nu}(k)$ to its rotated matrix $S_{\mu,\nu}(gk)$
    /// @param nbasis [in] global number of basis
    /// @param p2d [in]2d parallel info
    /// @param ucell [in] unitcell
    /// @param col_inside [in] whether the matrix is column-major (major means memory continuity)
    /// @param iat_rotiat [in]  atom index map of the symmetry operation g: g(iat0)=iat1
    /// @param sloc_ikibz [in] input matrix S(k) (local)
    /// @return S(gk) (local)
    std::vector<std::complex<double>> rearrange_col(
        const int& nbasis,
        const Parallel_2D& p2d,
        const UnitCell& ucell,
        const bool col_inside,
        const std::vector<int>& iat_rotiat,
        const std::vector<std::complex<double>>& sloc_in);


    /// @brief restore c_k from c_gk: $c_k=\tilde{S}^{-1}(k)S(gk)c_{gk}$ for all the ibz-kpoints gk
    /// @param nkstot_full  [in] number of kpoints in all the kstars
    /// @param psi_ikibz [in] c_gk: wavefunction of current ibz-kpoint
    /// @param sloc_ikibz [in] local overlap matrices of current ibz-kpoint: S(gk)
    /// @param sloc_ik [in] $\nu$-rearranged local overlap matrices of each k in current kstar : S(k)
    /// @param nbasis [in] global number of basis
    /// @param nbands [in] global number of bands
    /// @param pv [in] parallel orbitals (for both matrix and wavefunction)
    psi::Psi<std::complex<double>, psi::DEVICE_CPU> restore_psik(
        const int& nkstot_full,
        const psi::Psi<std::complex<double>, psi::DEVICE_CPU>& psi_ibz,
        const std::vector<std::vector<std::complex<double>>>& sloc_ibz,
        const std::vector<std::vector<std::vector<std::complex<double>>>>& sloc_full,
        const int& nbasis,
        const int& nbands,
        const Parallel_Orbitals& pv);
    ///@param invSkrot_Sgk [in] $\tilde{S}^{-1}(k)S(gk)$ for all the ibz-kpoints gk
    psi::Psi<std::complex<double>, psi::DEVICE_CPU> restore_psik(
        const int& nkstot_full,
        const psi::Psi<std::complex<double>, psi::DEVICE_CPU>& psi_ibz,
        const std::vector<std::vector<std::vector<std::complex<double>>>>& invSkrot_Sgk,
        const int& nbasis,
        const int& nbands,
        const Parallel_Orbitals& pv);

#ifdef __MPI
    /// @brief calculate $\tilde{S}^{-1}(k)S(gk)$ for one ibz-kpoint
    /// @param ikibz  [in] ibz-kpoint index
    /// @param sloc_ikibz [in] local overlap matrices of current ibz-kpoint: S(gk)
    /// @param sloc_ik [in] $\nu$-rearranged local overlap matrices of each k in current kstar : S(k)
    /// @param nbasis [in] global number of basis
    /// @param pv [in] 2d-block-cyclic info (only for matrix)
    std::vector<std::vector<std::complex<double>>> cal_invSkrot_Sgk_scalapack(
        const int& ikibz,
        const std::vector<std::complex<double>>& sloc_ikibz,
        const std::vector<std::vector<std::complex<double>>>& sloc_ik,
        const int& nbasis,
        const Parallel_2D& pv);

    /// @brief restore c_k from c_gk: $c_k=\tilde{S}^{-1}(k)S(gk)c_{gk}$ for one ibz-kpoint
    /// @param ikibz  [in] ibz-kpoint index
    /// @param ikfull_start [in] start index of k in all the kstars
    /// @param psi_ikibz [in] c_gk: wavefunction of ibz-kpoint
    /// @param sloc_ik $\nu$-rearranged local overlap matrices of each k in current kstar : S(k)
    /// @param nbasis [in] global number of basis
    /// @param nbands [in] global number of bands
    /// @param pv [in] parallel orbitals (for both matrix and wavefunction)
    /// @param col_inside [in] whether the matrix is column-major (major means memory continuity)
    /// @param psi_full: [out]wavefunction of each k in kstars[ikibz]
    void restore_psik_scalapack(
        const int& ikibz,
        const int& ikfull_start,
        const psi::Psi<std::complex<double>, psi::DEVICE_CPU>& psi_ikibz,
        const std::vector<std::complex<double>>& sloc_ikibz,
        const std::vector<std::vector<std::complex<double>>>& sloc_ik,
        const int& nbasis,
        const int& nbands,
        const Parallel_Orbitals& pv,
        psi::Psi<std::complex<double>, psi::DEVICE_CPU>* psi_full);

    /// @brief restore c_k from c_gk: $c_k=\tilde{S}^{-1}(k)S(gk)c_{gk}$ for one ibz-kpoint
    /// @param ikibz  [in] ibz-kpoint index
    /// @param ikfull_start [in] start index of k in all the kstars
    /// @param psi_ikibz [in] c_gk: wavefunction of ibz-kpoint
    /// @param invSkrot_Sgk [in] $S^{-1}(k)S(gk)$ of each k in current kstar
    /// @param nbasis [in] global number of basis
    /// @param nbands [in] global number of bands
    /// @param pv [in] parallel orbitals (for both matrix and wavefunction)
    /// @param col_inside [in] whether the matrix is column-major (major means memory continuity)
    /// @param psi_full: [out]wavefunction of each k in kstars[ikibz]
    void restore_psik_scalapack(
        const int& ikibz,
        const int& ikfull_start,
        const psi::Psi<std::complex<double>, psi::DEVICE_CPU>& psi_ikibz,
        const std::vector<std::vector<std::complex<double>>>& invSkrot_Sgk,
        const int& nbasis,
        const int& nbands,
        const Parallel_Orbitals& pv,
        psi::Psi<std::complex<double>, psi::DEVICE_CPU>* psi_full);
#endif

    std::vector<std::vector<std::complex<double>>> cal_invSkrot_Sgk_lapack(
        const int& ikibz,
        const std::vector<std::complex<double>>& sloc_ikibz,
        const std::vector<std::vector<std::complex<double>>>& sloc_ik,
        const int& nbasis);
    void restore_psik_lapack(
        const int& ikibz,
        const int& ikfull_start,
        const psi::Psi<std::complex<double>, psi::DEVICE_CPU>& psi_ikibz,
        const std::vector<std::complex<double>>& sloc_ikibz,
        const std::vector<std::vector<std::complex<double>>>& sloc_ik,
        const int& nbasis,
        const int& nbands,
        psi::Psi<std::complex<double>, psi::DEVICE_CPU>* psi_full);
    void restore_psik_lapack(
        const int& ikibz,
        const int& ikfull_start,
        const psi::Psi<std::complex<double>, psi::DEVICE_CPU>& psi_ikibz,
        const std::vector<std::vector<std::complex<double>>>& invSkrot_Sgk,
        const int& nbasis,
        const int& nbands,
        psi::Psi<std::complex<double>, psi::DEVICE_CPU>* psi_full);
    std::vector<std::complex<double>> get_full_smat(
        const std::vector<std::complex<double>>& locmat,
        const int& nbasis,
        const Parallel_2D& p2d,
        const bool col_inside);


    std::vector<ModuleBase::Vector3<double>> matfunc(
        const std::vector<ModuleBase::Vector3<double>>& vec,
        const ModuleBase::Matrix3& gmat,
        const std::function<ModuleBase::Vector3<double>(const ModuleBase::Vector3<double>&, const ModuleBase::Matrix3&)>& f);
}