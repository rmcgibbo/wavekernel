#ifndef _psi_wavekernel_matrixutils_h
#define _psi_wavekernel_matrixutils_h

#include <libmints/mints.h>

using namespace boost;
using namespace psi;
namespace psi{ namespace wavekernel {

void inplace_element_square(SharedMatrix& x) {
    for (int h=0; h < x->nirrep(); ++h) {
        for (int i=0; i < x->rowspi()[h]; ++i) {
            for (int j=0; j < x->colspi()[h]; ++j) {
                double e = x->get(h, i, j);
                x->set(h, i, j, e*e);
            }
        }
    }
}

}} /* namespace */
#endif
