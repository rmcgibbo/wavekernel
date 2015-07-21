#ifndef _psi_wavekernel_mosig_h
#define _psi_wavekernel_mosig_h
#include <libmints/mints.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <liboptions/liboptions.h>
#include <libscf_solver/sad.h>
#include "matrixutils.hpp"

using namespace boost;
namespace psi {
namespace wavekernel {

class MOSignature {
private:
    // Wavefunction
    shared_ptr<Wavefunction> wfn_;
    // R^3 point grid for spatial integrals
    shared_ptr<DFTGrid> grid_;
    shared_ptr<PointFunctions> wfn_properties_;
    shared_ptr<UKSFunctions> sad_properties_;

    // Global options
    Options& options_;

    // Superposition of atomic densities guess
    shared_ptr<scf::SADGuess> sad_guess_;

    // Gaussian blur of the orbitals by energy
    SharedMatrix orbital_mixing_a_;
    SharedMatrix orbital_mixing_b_;
    // Computed point signature vectors (only for one block of points)
    SharedMatrix v_;
    // Computed point classification (integers, only for one block of points)
    SharedVector s_;
    // MO energies, both alpha and beta together.
    SharedVector epsilon_;
    SharedVector epsilon_a_;
    SharedVector epsilon_b_;

    const int num_curve_;
    const int num_electrons_;

    void initialize_orbital_mixing_by_temperature();
    void initialize_orbital_mixing_by_chemical_potential();
    std::vector<std::vector<size_t> > sample_block_subset_indices(size_t n_samples);
    void check_basis(const SharedMatrix& basis);


public:
    void dump_npz(const std::string& filename) {
        save_npz(filename, orbital_mixing_a_, "orbtail_mixing_a", "w"); // "w" overwrites any existing file
        save_npz(filename, orbital_mixing_b_, "orbital_mixing_b", "a");
        save_npz(filename, epsilon_a_, "epsilon_a", "a");
        save_npz(filename, epsilon_b_, "epsilon_b", "a");
    }

    MOSignature(Options& options);

    void compute_v(int block);
    void compute_s(const SharedMatrix& basis, int block);
    SharedVector get_x(const SharedMatrix& basis);
    SharedMatrix sample_v(size_t n_samples, SharedMatrix coords);

    SharedMatrix v() const {
        return v_;
    }
    SharedVector s() const {
        return s_;
    }

    int nblocks() const {
        return grid_->blocks().size();
    }
    const std::vector<shared_ptr<BlockOPoints> >& blocks() const {
        return grid_->blocks();
    }
    shared_ptr<DFTGrid> grid() const {
        return grid_;
    }

};

}
}
#endif
