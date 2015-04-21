#include <boost/random.hpp>

#include "mosignature.hpp"
#include "matrixutils.hpp"


using namespace psi;
using namespace boost;
namespace psi{ namespace wavekernel{


MOSignature::MOSignature(Options& options) :
    options_(options),
    num_energies_(options_.get_int("E_DESCRIPTOR_NUM")),
    wfn_(Process::environment.wavefunction()),
    molecule_(Process::environment.molecule())
{
    if (!wfn_) {
        outfile->Printf("SCF has not been run yet!\n");
        throw PSIEXCEPTION("SCF has not been run yet!\n");
    }

    grid_ = shared_ptr<DFTGrid>(new DFTGrid(molecule_, wfn_->basisset(), options_));
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    // [TODO: deal with unrestricted wavefunctions]
    properties_ = shared_ptr<PointFunctions>(
        new RKSFunctions(wfn_->basisset(), max_points, max_functions));
    properties_->set_ansatz(0);
    properties_->set_Cs(wfn_->Ca());

    orbital_blur_ = SharedMatrix(new Matrix("Orbital Blur", num_energies_, wfn_->nmo()));
    v_ = SharedMatrix(new Matrix("V_BLOCK", num_energies_, max_points));

    orbital_blur_->zero();
    // [TODO]: Deal with unrestricted wavefunctions.
    initialize_orbital_blur(wfn_->epsilon_a());

}


void MOSignature::initialize_orbital_blur(SharedVector epsilon) {
    const double e_min = options_.get_double("E_DESCRIPTOR_MIN");
    const double e_max = options_.get_double("E_DESCRIPTOR_MAX");
    const double sigma = options_.get_double("E_DESCRIPTOR_SIGMA");
    const double prefactor = 1 / (sigma * std::sqrt(2*M_PI));

    // [TODO] assert that len(epsilon_a) == wfn->nmo();
    // Do we want to always include all of the orbitals?
    // Maybe we should include the virtual orbitals?


    for (int i = 0; i < num_energies_; i++) {
        double e_i = e_min + ((e_max - e_min) / (num_energies_-1)) * i;

        for (int j = 0; j < wfn_->nmo(); j++) {
            double exponent = -std::pow(epsilon->get(j) - e_i, 2) / (2*sigma*sigma);
            double overlap = prefactor * std::exp(exponent);
            orbital_blur_->set(i, j, orbital_blur_->get(i,j)+overlap);
        }
    }
}


void MOSignature::compute_v(int Q) {
    shared_ptr<BlockOPoints> block = blocks()[Q];
    properties_->compute_orbitals(block);
    SharedMatrix psi_a = properties_->orbital_value("PSI_A");
    inplace_element_square(psi_a);
    v_->zero();
    v_->accumulate_product(orbital_blur_, psi_a);
}


SharedMatrix MOSignature::sample_v(size_t n_samples) {
    SharedMatrix v_samples = SharedMatrix(
        new Matrix("V_BLOCK", num_energies_, n_samples));
    std::vector<std::vector<size_t> > block_indices = \
        sample_block_subset_indices(n_samples);

    size_t j = 0;
    for (int Q = 0; Q < nblocks(); Q++) {
        compute_v(Q);
        for (size_t i = 0; i < block_indices[Q].size(); i++, j++) {
            SharedVector col = v_->get_column(0, block_indices[Q][i]);
            v_samples->set_column(0, j, col);
        }
    }

    return v_samples;
}


/*
    Generate n_samples random indices for points on the molecular grid,
    broken into blocks.

    The return value is a list of vectors of length equal to nblocks
    where each element is a vector of the selected indices in that block, (in
    the internal indexing scheme of the block, not the global indexing)
*/
std::vector<std::vector<size_t> >
MOSignature::sample_block_subset_indices(size_t n_samples)
{
    size_t npoints = 0;  // total number of points in the grid
    std::vector<size_t> block_start;  // starting index of each block
    for (size_t Q = 0; Q < blocks().size(); Q++) {
        block_start.push_back(npoints);
        npoints += blocks()[Q]->npoints();
    }
    block_start.push_back(npoints);

    boost::random::mt19937 gen(time(0));
    boost::random::uniform_int_distribution<size_t> dist(0, npoints-1);
    std::vector<size_t> samples;
    for (size_t i = 0; i < n_samples; i++) {
        samples.push_back(dist(gen));
    }
    std::sort(samples.begin(), samples.end());
    std::vector<std::vector<size_t> > block_indices(blocks().size());

    size_t j = 0;
    for (size_t i = 0; i < n_samples; i++) {
        if (samples[i] >= block_start[j+1]) {
            j += 1;
        }
        block_indices[j].push_back(samples[i] - block_start[j]);
    }

    return block_indices;
}



}} /* namespace */
