#include <boost/random.hpp>

#include "mosignature.hpp"
#include "matrixutils.hpp"
#include "fermilevel.hpp"


using namespace psi;
using namespace boost;
namespace psi{ namespace wavekernel{

/* Boltzmann constant in [Hartree / K] */
static const double BOLTZMANN = 3.166811429e-6;

MOSignature::MOSignature(Options& options) :
    options_(options),
    num_curve_(options_.get_int("NUM_CURVE")),
    wfn_(Process::environment.wavefunction()),
    num_electrons_(wfn_->nalpha() + wfn_->nbeta())
{
    if (!wfn_) {
        throw PSIEXCEPTION("SCF has not been run yet!\n");
    }

    if (wfn_->nirrep() != 1) {
        throw PSIEXCEPTION("This plugin is only implemented in the C1 symmetry group.");
    }

    shared_ptr<Molecule> mol = Process::environment.molecule();
    grid_ = shared_ptr<DFTGrid>(new DFTGrid(mol, wfn_->basisset(), options_));
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();

    std::string ref = options_.get_str("REFERENCE");
    if ((ref == "RKS") || (ref == "RHF")) {
        wfn_properties_ = shared_ptr<PointFunctions>(
            new RKSFunctions(wfn_->basisset(), max_points, max_functions));
        wfn_properties_->set_ansatz(0);
        wfn_properties_->set_Cs(wfn_->Ca());
        wfn_properties_->set_pointers(wfn_->Da());
    } else if ((ref == "UKS") || (ref == "UHF")) {
        wfn_properties_ = shared_ptr<PointFunctions>(
            new UKSFunctions(wfn_->basisset(), max_points, max_functions));
        wfn_properties_->set_ansatz(0);
        wfn_properties_->set_Cs(wfn_->Ca(), wfn_->Cb());
        wfn_properties_->set_pointers(wfn_->Da(), wfn_->Db());
    } else {
        throw PSIEXCEPTION("Unknown reference: " + ref + "\n");
    }

    sad_guess_ = shared_ptr<scf::SADGuess>(new scf::SADGuess(wfn_->basisset(), wfn_->nalpha(), wfn_->nbeta(), options_));
    sad_guess_->compute_guess();
    sad_properties_ = shared_ptr<UKSFunctions>(new UKSFunctions(wfn_->basisset(), max_points, max_functions));
    sad_properties_->set_ansatz(0);
    sad_properties_->set_Cs(sad_guess_->Ca(), sad_guess_->Cb());
    sad_properties_->set_pointers(sad_guess_->Da(), sad_guess_->Db());

    orbital_mixing_a_ = SharedMatrix(new Matrix("Orbital mixing A", num_curve_, wfn_->nmo()));
    orbital_mixing_b_ = SharedMatrix(new Matrix("Orbital mixing B", num_curve_, wfn_->nmo()));

    v_ = SharedMatrix(new Matrix("V_BLOCK", num_curve_, max_points));
    s_ = SharedVector(new Vector("S_BLOCK", max_points));

    // copy wfn_->epsilon_a() and wfn_->epsilon_b() into a single vector
    epsilon_ = SharedVector(new Vector(wfn_->epsilon_a()->dim() + wfn_->epsilon_b()->dim()));
    for (int i = 0; i < wfn_->epsilon_a()->dim(); i++) {
        epsilon_->set(i, wfn_->epsilon_a()->get(i));
    }
    for (int i = 0; i < wfn_->epsilon_b()->dim(); i++) {
        epsilon_->set(i + wfn_->epsilon_a()->dim(), wfn_->epsilon_b()->get(i));
    }

    if (options_.get_str("CURVE").compare("TEMP") == 0) {
        initialize_orbital_mixing_by_temperature();
    } else if (options_.get_str("CURVE").compare("MU") == 0) {
        initialize_orbital_mixing_by_chemical_potential();
    } else {
        throw PSIEXCEPTION("Unknown curve: " + options_.get_str("CURVE") + "\n");
    }
}

void MOSignature::initialize_orbital_mixing_by_chemical_potential() {
    const double mu_min = options_.get_double("CURVE_MIN");
    const double mu_max = options_.get_double("CURVE_MAX");
    const double temp = options_.get_double("TEMP");
    const double beta = 1 / (BOLTZMANN * temp);
    outfile->Printf("\nInitializing Orbital mixing by chemical potential\n");

    for (int i = 0; i < num_curve_; i++) {
        double n_electrons = 0;
        double mu = mu_min + ((mu_max - mu_min) / (num_curve_-1)) * i;
        for (int j = 0; j < wfn_->nmo(); j++) {
            double na = n_occ(wfn_->epsilon_a()->get(j), mu, beta);
            double nb = n_occ(wfn_->epsilon_b()->get(j), mu, beta);
            n_electrons += (na+nb);
            orbital_mixing_a_->set(i, j, na);
            orbital_mixing_b_->set(i, j, nb);
        }
        outfile->Printf("  mu=%f n_electrons=%f (%d)\n", mu, n_electrons, (wfn_->nalpha() + wfn_->nbeta()));
    }
}

void MOSignature::initialize_orbital_mixing_by_temperature() {
    const double T_min = options_.get_double("CURVE_MIN");
    const double T_max = options_.get_double("CURVE_MAX");
    outfile->Printf("\nInitializing orbital mixing by temperature\n");

    // orbital_mixing_a_ is of dimension (num_curve_ x num_molecular_orbitals)

    for (int i = 0; i < num_curve_; i++) {
        double T = T_min + ((T_max - T_min) / (num_curve_-1)) * i;
        double beta = 1 / (BOLTZMANN * T);
        // compute the fermi level at this temperature
        double mu = calculate_mu(num_electrons_, 1 / (BOLTZMANN * T), epsilon_);

        outfile->Printf("Temp:  %.3f K\n", T);
        outfile->Printf("beta:  %.3f 1/h\n", beta);
        outfile->Printf("fermi: %.3f h\n\n", mu);

        for (int j = 0; j < wfn_->nmo(); j++) {
            double na = n_occ(wfn_->epsilon_a()->get(j), mu, beta);
            double nb = n_occ(wfn_->epsilon_b()->get(j), mu, beta);
            //printf("%f ", na);

            orbital_mixing_a_->set(i, j, na);
            orbital_mixing_b_->set(i, j, nb);
        }
    }
    // orbital_mixing_a_->print();
}


// void MOSignature::initialize_orbital_mixing(SharedVector epsilon) {
//     const double e_min = options_.get_double("E_DESCRIPTOR_MIN");
//     const double e_max = options_.get_double("E_DESCRIPTOR_MAX");
//     const double sigma = options_.get_double("E_DESCRIPTOR_SIGMA");
//     const double prefactor = 1 / (sigma * std::sqrt(2*M_PI));
//
//     // [TODO] assert that len(epsilon_a) == wfn->nmo();
//     // Do we want to always include all of the orbitals?
//     // Maybe we should include the virtual orbitals?
//
//
//     for (int i = 0; i < num_energies_; i++) {
//         double e_i = e_min + ((e_max - e_min) / (num_energies_-1)) * i;
//
//         for (int j = 0; j < wfn_->nmo(); j++) {
//             double exponent = -std::pow(epsilon->get(j) - e_i, 2) / (2*sigma*sigma);
//             double overlap = prefactor * std::exp(exponent);
//             orbital_blur_->set(i, j, orbital_blur_->get(i,j)+overlap);
//         }
//     }
// }


/* Compute the local descriptor vectors for each point in the molecular grid's
   Qth block.
*/
void MOSignature::compute_v(int Q) {
    shared_ptr<BlockOPoints> block = blocks()[Q];
    wfn_properties_->compute_orbitals(block);

    // [number of MOs x n_max_points]
    SharedMatrix psi_a = wfn_properties_->orbital_value("PSI_A");
    SharedMatrix psi_b = wfn_properties_->orbital_value("PSI_B");

    inplace_element_square(psi_a);
    if (psi_a != psi_b) {
        // unrestricted wavefunction
        inplace_element_square(psi_b);
    }

    v_->accumulate_product(orbital_mixing_a_, psi_a);
    v_->accumulate_product(orbital_mixing_b_, psi_b);

    SharedVector sad_rho_a = sad_properties_->point_value("RHO_A");
    SharedVector sad_rho_b = sad_properties_->point_value("RHO_B");

    // subtract out the SAD density baseline from every point
    double** p = v_->pointer(0);
    for (int j = 0; j < v_->colspi(0); j++) {
        double offset = sad_rho_a->get(j) + sad_rho_b->get(j);
        for (int i = 0; i < v_->rowspi(0); i++) {
            p[i][j] = p[i][j] - offset;
        }
    }
}


SharedMatrix MOSignature::sample_v(size_t n_samples) {
    SharedMatrix v_samples = SharedMatrix(
        new Matrix("V_BLOCK", n_samples, num_curve_));
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
    if (basis->colspi(0) != num_curve_)
        throw PSIEXCEPTION("Basis must have n_cols match num_curve_.\n");
}


void MOSignature::compute_s(const SharedMatrix& basis, int Q) {
    check_basis(basis);
    compute_v(Q);
    shared_ptr<BlockOPoints> block = blocks()[Q];
    wfn_properties_->compute_points(block);

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
        double* rho_a = wfn_properties_->point_value("RHO_A")->pointer();
        double* rho_b = wfn_properties_->point_value("RHO_A")->pointer();
        double* w = block->w();
        for (int i = 0; i < block->npoints(); i++) {
            double xx = (rho_a[i] + rho_b[i]) * w[i];
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

    static uint32_t seed = 0;
    if (seed == 0) {
      FILE *dev_urandom;
      if((dev_urandom = fopen("/dev/urandom", "r")) == NULL)
	throw PSIEXCEPTION("/dev/urandom read error");
      if(fread(&seed, sizeof(seed), 1, dev_urandom) == 0)
	throw PSIEXCEPTION("/dev/unrandom read error");
      if(fclose(dev_urandom) != 0)
	throw PSIEXCEPTION("/dev/unrandom read error");
    }
    static boost::random::mt19937 gen(seed);
    //boost::random::mt19937 gen(0.0);

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
