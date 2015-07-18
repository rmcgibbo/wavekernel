#include "Python.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <boost/filesystem.hpp>
#include <psi4-dec.h>
#include <libplugin/plugin.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include <libfock/cubature.h>
#include <libfock/points.h>

#include "matrixutils.hpp"
#include "mosignature.hpp"
#include "fermilevel.hpp"
#include "gitversion.hpp"


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
/
        options.add_double("TEMP", 1000);
        options.add_double("CURVE_MIN", 0);
        options.add_double("CURVE_MAX", 5000);
        options.add_int("NUM_CURVE", 10);

        options.add_str("MODE", "SAMPLE_V", "SAMPLE_V CALCULATE_X");
        options.add_str("CURVE", "TEMP", "TEMP MU");
        options.add_int("NUM_SAMPLE_DESCRIPTORS", 100);
        options.add_str("FILENAME_IN", "descriptors.npy");
        options.add_str("FILENAME_OUT", "NULL");
    }

    return true;
}


extern "C"
PsiReturnType wavekernel(Options& options)
{
    PsiReturnType status;
    shared_ptr<MOSignature> mosig = shared_ptr<MOSignature>(new MOSignature(options));

    std::string fn_in = options.get_str("FILENAME_IN");
    std::transform(fn_in.begin(), fn_in.end(), fn_in.begin(), ::tolower);
    std::string fn_out = options.get_str("FILENAME_OUT");
    std::transform(fn_out.begin(), fn_out.end(), fn_out.begin(), ::tolower);

    string mode = options.get_str("MODE");
    outfile->Printf("\n ======= Molecular Orbital Signature Plugin =======\n");
    outfile->Printf("    gitversion   = %s\n", GIT_VERSION);
    outfile->Printf("    mode         = %s\n", mode.c_str());
    outfile->Printf("    filename_in  = %s\n", fn_in.c_str());
    outfile->Printf("    filename_out = %s\n\n", fn_out.c_str());

    if (mode == "SAMPLE_V") {
        if (fn_out == "null") {
            outfile->Printf("Error: must set filename_out option");
            return Failure;
        }

        int num_samples = options.get_int("NUM_SAMPLE_DESCRIPTORS");
        SharedMatrix v_samples = mosig->sample_v(num_samples);
        outfile->Printf("Writing %d wavekernel descriptors, v, to %s\n", num_samples, fn_out.c_str());
        save_npy(fn_out, v_samples);
    } else if (mode == "DUMP_S") {
        SharedMatrix output = SharedMatrix(new Matrix("output", mosig->grid()->npoints(), 4));
        if (!boost::filesystem::exists(fn_in)) {
            outfile->Printf("Error: file does not exist: %s", fn_in.c_str());
            return Failure;
        }
        SharedMatrix basis = load_npy(fn_in);

        int kk = 0;
        for (int Q = 0; Q < mosig->nblocks(); Q++) {
            shared_ptr<BlockOPoints> block = mosig->blocks()[Q];
            mosig->compute_s(basis, Q);
            for (int k = 0; k < block->npoints(); k++, kk++) {
                output->set(kk, 0, block->x()[k]);
                output->set(kk, 1, block->y()[k]);
                output->set(kk, 2, block->z()[k]);
                output->set(kk, 3, mosig->s()->pointer()[k]);
            }
        }
        save_npy(fn_out, output);

    } else if (mode == "CALCULATE_X") {
        if (fn_out == "null") {
            outfile->Printf("Error. Must set filename_out option");
            return Failure;
        }
        if (!boost::filesystem::exists(fn_in)) {
            outfile->Printf("Error: file does not exist: %s", fn_in.c_str());
            return Failure;
        }

        SharedMatrix basis = load_npy(fn_in);
        SharedVector x = mosig->get_x(basis);
        outfile->Printf("Writing molecular vector, x, to %s\n", fn_out.c_str());
        save_npy(fn_out, x);
    } else {
        return Failure;
    }

    outfile->Printf("\nFINISHED WAVEKERNEL WITHOUT DEATH!\n");
    return Success;
}


}} // End namespaces
