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
// #include "utils.hpp"
#include "mosignature.hpp"


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
        options.add_str("MODE", "SAMPLE_V", "SAMPLE_V CALCULATE_X");
        options.add_int("NUM_SAMPLE_DESCRIPTORS", 100);
        options.add_str("FILENAME", "descriptors.npy");
    }

    return true;
}


extern "C"
PsiReturnType wavekernel(Options& options)
{
    PsiReturnType status;
    shared_ptr<MOSignature> mosig = shared_ptr<MOSignature>(new MOSignature(options));

    std::string fn = options.get_str("FILENAME");
    std::transform(fn.begin(), fn.end(), fn.begin(), ::tolower);

    string mode = options.get_str("MODE");
    outfile->Printf("\n ======= Molecular Orbital Signature Plugin =======\n");
    outfile->Printf("    mode     = %s\n", mode.c_str());
    outfile->Printf("    filename = %s", fn.c_str());
    outfile->Printf("\n");

    if (mode == "SAMPLE_V") {
        int num_samples = options.get_int("NUM_SAMPLE_DESCRIPTORS");
        SharedMatrix v_samples = mosig->sample_v(num_samples);
        outfile->Printf("Writing %d wavekernel descriptors, v, to %s\n", num_samples, fn.c_str());
        save_npy(fn, v_samples);
    } else if (mode == "CALCULATE_X") {
        SharedMatrix basis = load_npy(fn);
        SharedVector x = mosig->get_x(basis);
        x->print();
    } else {
        return Failure;
    }




    outfile->Printf("\nFINISHED WAVEKERNEL WITHOUT DEATH!\n");
    return Success;
}


}} // End namespaces

