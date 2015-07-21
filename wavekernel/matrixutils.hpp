#ifndef _psi_wavekernel_matrixutils_h
#define _psi_wavekernel_matrixutils_h

#include <libmints/mints.h>

using namespace boost;
using namespace psi;
namespace psi {
namespace wavekernel {

/*
    Square every entry of a matrix inplace
*/
void inplace_element_square(const SharedMatrix& x);
int assign(const SharedVector& v, const SharedMatrix& X);
void save_npy(const std::string& file, const SharedMatrix& arr);
void save_npy(const std::string& file, const SharedVector& arr);
void save_npz(const std::string& file, const SharedMatrix& arr, const std::string& name, std::string mode="w");
void save_npz(const std::string& file, const SharedVector& arr, const std::string& name, std::string mode="w");
SharedMatrix load_npy(const std::string& file);

}
} /* namespace */
#endif
