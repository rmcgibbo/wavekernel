#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

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

void sample_descriptors() {
    fprintf(outfile, "WAVE KERNEL PLUGIN\n");
    shared_ptr<Molecule> molecule = Process::environment.molecule();
}


extern "C"
PsiReturnType wavekernel(Options& options)
{
    int print = options.get_int("PRINT");
    string mode = options.get_str("MODE");
    if (mode == "SAMPLE_DESCRIPTORS") {
        sample_descriptors();
        return Success;
    } else {
        return Failure;
    }




}

}} // End namespaces

