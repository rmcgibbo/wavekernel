#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include <libfock/cubature.h>
#include <libfock/points.h>

INIT_PLUGIN

using std::string;
using namespace boost;

namespace psi{ namespace wavekernel {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "WAVEKERNEL"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        options.add_str("MODE", "", "Wave kernel plugin mode");
    }

    return true;
}

void sample_descriptors(Options& options) {
    outfile->Printf("WAVE KERNEL PLUGIN\n");
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    boost::shared_ptr<BasisSet> basisset = wfn->basisset();

    if(!wfn) {
        outfile->Printf("SCF has not been run yet!\n");
        throw PSIEXCEPTION("SCF has not been run yet!\n");
    }

    outfile->Printf("Building grid\n");
    boost::shared_ptr<DFTGrid> grid = boost::shared_ptr<DFTGrid>(new DFTGrid(molecule, basisset, options));
    // fprintf(outfile, "Done with grid!\n");
    grid->print("outfile", 1);

    int max_points = grid->max_points();
    int max_functions = grid->max_functions();
    boost::shared_ptr<RKSFunctions> properties = \
         boost::shared_ptr<RKSFunctions>(new RKSFunctions(basisset, max_points, max_functions));
    properties->set_ansatz(0);
    SharedMatrix Ca = wfn->Ca();
    properties->set_Cs(Ca);

    // outfile->Printf("nirrep = %d", Ca->nirrep());
    // for (int h = 0; h < Ca->nirrep(); h++) {
    //     int nso = Ca->rowspi()[h];
    //     int nmo = Ca->colspi()[h];
    //     outfile->Printf("n_so=%d  n_mo=%d\n", nso, nmo);
    // }
    //
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    // outfile->Printf("Number of blocks: %lu\n", blocks.size());
    for (size_t Q = 0; Q < blocks.size(); Q++) {
        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        // int npoints = block->npoints();
        // double *restrict x = block->x();
        // double *restrict y = block->y();
        // double *restrict z = block->z();
        // double *restrict w = block->w();
        // outfile->Printf("Number of points in block: %d\n", npoints);
        // outfile->Printf("[%f %f %f %f]\n", x[0], y[0], z[0], w[0]);

        // const std::vector<int>& function_map = block->functions_local_to_global();
        // int nglobal = max_functions;
        // int nlocal  = function_map.size();
        // outfile->Printf("n_local=%d n_global=%d\n", nlocal, nglobal);

        properties->compute_orbitals(block);

        // this is the object we want!
        properties->orbital_value("PSI_A")->print();


        break;
    }


    // properties->point_value

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

