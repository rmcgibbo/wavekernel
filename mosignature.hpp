#ifndef _psi_wavekernel_mosig_h
#define _psi_wavekernel_mosig_h
#include <libmints/mints.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <liboptions/liboptions.h>

using namespace boost;
namespace psi{ namespace wavekernel {

class MOSignature {
private:
    // Wavefunction
    shared_ptr<Wavefunction> wfn_;
    // R^3 point grid for spatial integrals
    shared_ptr<DFTGrid> grid_;
    shared_ptr<PointFunctions> properties_;
    // Global options
    Options& options_;

    // Gaussian blur of the orbitals by energy
    SharedMatrix orbital_blur_;
    // Computed point signature vectors (only for one block of points)
    SharedMatrix v_;

    const int num_energies_;

    void initialize_orbital_blur(SharedVector epsilon);
    std::vector<std::vector<size_t> > sample_block_subset_indices(size_t n_samples);

public:
    MOSignature(Options& options);

    void compute_v(int block);
    SharedVector get_x(SharedMatrix basis);
    SharedMatrix sample_v(size_t n_samples);

    int nblocks() const { return grid_->blocks().size(); }
    SharedMatrix v() const { return v_; }

    const std::vector<shared_ptr<BlockOPoints> >& blocks() const { return grid_->blocks(); }
};

}}
#endif
