#define _USE_MATH_DEFINES
#include <cmath>
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include <libfock/cubature.h>
#include <libfock/points.h>

#include "matrixutils.hpp"
#include "utils.hpp"

INIT_PLUGIN

using std::string;
using boost::shared_ptr;
using namespace boost;

namespace psi{ namespace wavekernel {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "WAVEKERNEL"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        options.add_double("E_DESCRIPTOR_MIN", -20);
        options.add_double("E_DESCRIPTOR_MAX", 0);
        options.add_int("E_DESCRIPTOR_NUM", 20);
        options.add_double("E_DESCRIPTOR_SIGMA", 2);
        options.add_str("MODE", "", "Wave kernel plugin mode");
        options.add_int("NUM_SAMPLE_DESCRIPTORS", 100);
        options.add_str("DESCRIPTOR_FN", "descriptors.dat");
    }

    return true;
}


SharedMatrix build_orbital_blur(SharedVector epsilon, Options& options) {
    shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

    const int e_num = options.get_int("E_DESCRIPTOR_NUM");
    const double e_min = options.get_double("E_DESCRIPTOR_MIN");
    const double e_max = options.get_double("E_DESCRIPTOR_MAX");
    const double sigma = options.get_double("E_DESCRIPTOR_SIGMA");
    const double prefactor = 1 / (sigma * std::sqrt(2*M_PI));

    SharedMatrix mixing = SharedMatrix(new Matrix("Orbital Mixing", e_num, wfn->nmo()));

    // [TODO] assert that len(epsilon_a) == wfn->nmo();
    // Do we want to always include all of the orbitals?
    // Maybe we should include the virtual orbitals?
    // [TODO]: Deal with unrestricted wavefunctions.

    for (int i = 0; i < e_num; i++) {
        double e_i = e_min + ((e_max - e_min) / (e_num-1)) * i;
        // outfile->Printf("e_i[%d]=%f\n", i, e_i);

        for (int j = 0; j < wfn->nmo(); j++) {
            double exponent = -std::pow(epsilon->get(j) - e_i, 2) / (2*sigma*sigma);
            double overlap = prefactor * std::exp(exponent);
            mixing->set(i, j, overlap);
        }

    }
    return mixing;
}


void sample_descriptors(Options& options) {
    outfile->Printf("WAVE KERNEL PLUGIN\n");
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    shared_ptr<BasisSet> primary = wfn->basisset();

    if (!wfn) {
        outfile->Printf("SCF has not been run yet!\n");
        throw PSIEXCEPTION("SCF has not been run yet!\n");
    }

    shared_ptr<DFTGrid> grid = shared_ptr<DFTGrid>(new DFTGrid(molecule, primary, options));
    grid->print("outfile", 1);
    int max_points = grid->max_points();
    int max_functions = grid->max_functions();

    // [TODO: deal with unrestricted wavefunctions]
    shared_ptr<PointFunctions> properties = \
         shared_ptr<PointFunctions>(new RKSFunctions(primary, max_points, max_functions));
    properties->set_ansatz(0);
    SharedMatrix Ca = wfn->Ca();
    properties->set_Cs(Ca);

    SharedMatrix blur = build_orbital_blur(wfn->epsilon_a(), options);
    const std::vector<shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    std::vector<std::vector<size_t> > block_indices = \
        blocks_rand_subset(blocks, options.get_int("NUM_SAMPLE_DESCRIPTORS"));

    SharedMatrix vm = shared_ptr<Matrix>(new Matrix("V_BLOCK", blur->rowspi(0), max_points));

    FILE* fh = fopen("vectors.dat", "a");

    for (size_t Q = 0; Q < blocks.size(); Q++) {
        shared_ptr<BlockOPoints> block = blocks[Q];

        properties->compute_orbitals(block);
        SharedMatrix PsiA = properties->orbital_value("PSI_A");
        inplace_element_square(PsiA);
        vm->zero();
        vm->accumulate_product(blur, PsiA);

        for (size_t i = 0; i < block_indices[Q].size(); i++) {
            SharedVector v = vm->get_column(0, block_indices[Q][i]);

            for (size_t j = 0; j < v->dim(0); j++) {
                fprintf(fh, "%12.8f ", v->get(j));
            }
            fprintf(fh, "\n");
        }
    }
    fclose(fh);
}


extern "C"
PsiReturnType wavekernel(Options& options)
{
    PsiReturnType status;
    try {
        int print = options.get_int("PRINT");
        string mode = options.get_str("MODE");
        if (mode == "SAMPLE_DESCRIPTORS") {
            sample_descriptors(options);
            status = Success;
        } else {
            status = Failure;
        }
    }
    catch (psi::PsiException e) {
        outfile->Printf( "====== ERROR ====\n%s %d %s", e.what(), e.line(), e.location());
        outfile->Printf("%s", e.what());
        throw e;
    }

    outfile->Printf("\nFINISHED WAVEKERNEL WITHOUT DEATH!\n");
    return status;
}

}} // End namespaces

