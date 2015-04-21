#ifndef _psi_wavekernel_matrixutils_h
#define _psi_wavekernel_matrixutils_h

#include <libmints/mints.h>

using namespace boost;
using namespace psi;
namespace psi{ namespace wavekernel {

/*
    Square every entry of a matrix inplace
*/
void inplace_element_square(SharedMatrix& x);
void save_npy(const std::string& file, SharedMatrix& arr);
SharedMatrix load_npy(const std::string& file);

}} /* namespace */
#endif
