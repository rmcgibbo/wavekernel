#include "matrixutils.hpp"
#include "cnpy/cnpy.h"

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


void save_npy(const std::string& file, SharedMatrix& arr) {
    if (arr->nirrep() != 1) {
        throw PSIEXCEPTION("Only matrices with 1 irrep are supported!\n");
    }

    const unsigned int shape[] = {
        static_cast<unsigned int>(arr->rowspi(0)),
        static_cast<unsigned int>(arr->colspi(0))};
    cnpy::npy_save(file, arr->get_pointer(), shape, 2, "w");
}

SharedMatrix load_npy(const std::string& file) {
    cnpy::NpyArray arr = cnpy::npy_load(file);
    if (arr.word_size != sizeof(double)) {
        throw PSIEXCEPTION("Only matrices of type double are supported.\n");
    }
    if (arr.shape.size() != 2) {
        throw PSIEXCEPTION("Only 2D matrices are supported.");
    }
    double* data = reinterpret_cast<double*>(arr.data);


    SharedMatrix out = SharedMatrix(
        new Matrix(file, arr.shape[0], arr.shape[1]));
    memcpy(out->get_pointer(), data, sizeof(double)*arr.shape[0]*arr.shape[1]);
    return out;
}


}}
