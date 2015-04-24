#include <boost/random.hpp>

#include "mosignature.hpp"
#include "matrixutils.hpp"


using namespace psi;
using namespace boost;
namespace psi{ namespace wavekernel{


MOSignature::MOSignature(Options& options) :
    options_(options),
    num_energies_(options_.get_int("E_DESCRIPTOR_NUM")),
    wfn_(Process::environment.wavefunction())
{
    if (!wfn_) {
        outfile->Printf("SCF has not been run yet!\n");
        throw PSIEXCEPTION("SCF has not been run yet!\n");
    }

    shared_ptr<Molecule> mol = Process::environment.molecule();
    grid_ = shared_ptr<DFTGrid>(new DFTGrid(mol, wfn_->basisset(), options_));
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    // [TODO: deal with unrestricted wavefunctions]
    properties_ = shared_ptr<PointFunctions>(
        new RKSFunctions(wfn_->basisset(), max_points, max_functions));
    properties_->set_ansatz(0);
    properties_->set_Cs(wfn_->Ca());
    properties_->set_pointers(wfn_->Da());

    orbital_blur_ = SharedMatrix(new Matrix("Orbital Blur", num_energies_, wfn_->nmo()));
    v_ = SharedMatrix(new Matrix("V_BLOCK", num_energies_, max_points));
    s_ = SharedVector(new Vector("S_BLOCK", max_points));


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
        new Matrix("V_BLOCK", n_samples, num_energies_));
    std::vector<std::vector<size_t> > block_indices = \
        sample_block_subset_indices(n_samples);

    size_t j = 0;
    for (int Q = 0; Q < nblocks(); Q++) {
        compute_v(Q);
        for (size_t i = 0; i < block_indices[Q].size(); i++, j++) {
            SharedVector vv = v_->get_column(0, block_indices[Q][i]);
            v_samples->set_row(0, j, vv);
        }
    }

    return v_samples;
}

void MOSignature::check_basis(const SharedMatrix& basis) {
    if (basis->nirrep() != 1)
        throw PSIEXCEPTION("Basis must have 1 irrep.\n");
    if (basis->colspi(0) != num_energies_)
        throw PSIEXCEPTION("Basis must have n_cols match num_energies.\n");
}

void MOSignature::compute_s(const SharedMatrix& basis, int Q) {
    check_basis(basis);
    compute_v(Q);
    shared_ptr<BlockOPoints> block = blocks()[Q];
    properties_->compute_points(block);

    for (int i = 0; i < block->npoints(); i++) {
        int k = assign(v_->get_column(0, i), basis);
        s_->set(0, i, static_cast<double>(k));
    }
}

SharedVector MOSignature::get_x(const SharedMatrix& basis) {
    check_basis(basis);
    int n_basis = basis->rowspi(0);
    SharedVector x = SharedVector(new Vector(n_basis));

    for (int Q = 0; Q < nblocks(); Q++) {
        compute_s(basis, Q);
        shared_ptr<BlockOPoints> block = blocks()[Q];
        for (int i = 0; i < block->npoints(); i++) {
            double xx = properties_->point_value("RHO_A")->get(i) * block->w()[i];

            int s_i = static_cast<int>(s_->get(i));
            x->set(s_i, xx + x->get(s_i));
        }
    }

    return x;
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
