#include "matrixutils.hpp"
#include "cnpy/cnpy.h"

using namespace boost;
using namespace psi;
namespace psi{ namespace wavekernel {


void inplace_element_square(const SharedMatrix& x) {
    for (int h=0; h < x->nirrep(); ++h) {
        for (int i=0; i < x->rowspi()[h]; ++i) {
            for (int j=0; j < x->colspi()[h]; ++j) {
                double e = x->get(h, i, j);
                x->set(h, i, j, e*e);
            }
        }
    }
}

int assign(const SharedVector& v, const SharedMatrix& X) {
    if (v->dim(0) != X->colspi(0)) {
         throw PSIEXCEPTION("len(v) must be equal to the 2nd dimension of X");
    }

    int closest_i = 0;
    double closest_d2 = std::numeric_limits<double>::max();

    for (int i = 0; i < X->rowspi(0); i++) {
        double d2 = 0;
        for (int j = 0; j < v->dim(0); j++) {
            d2 += std::pow(X->get(0, i, j) - v->get(0, j), 2);
        }
        if (d2 < closest_d2) {
            closest_d2 = d2;
            closest_i = i;
        }
    }

    // if v is so large that it has distance inf to everything,
    // it'll just get assigned to 0.
    return closest_i;
}


void save_npy(const std::string& file, const SharedMatrix& arr) {
    if (arr->nirrep() != 1) {
        throw PSIEXCEPTION("Only matrices with 1 irrep are supported!\n");
    }

    const unsigned int shape[] = {
        static_cast<unsigned int>(arr->rowspi(0)),
        static_cast<unsigned int>(arr->colspi(0))};
    cnpy::npy_save(file, arr->get_pointer(), shape, 2, "w");
}

void save_npy(const std::string& file, const SharedVector& arr) {
    if (arr->nirrep() != 1) {
        throw PSIEXCEPTION("Only vectors with 1 irrep are supported!\n");
    }

    const unsigned int shape[] = {
        static_cast<unsigned int>(arr->dimpi()[0]),};

    cnpy::npy_save(file, arr->pointer(0), shape, 1, "w");
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
