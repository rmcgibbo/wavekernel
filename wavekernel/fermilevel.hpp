#ifndef _psi_wavekernel_fermilevel_h
#define _psi_wavekernel_fermilevel_h

#include <vector>
#include <libmints/mints.h>
#include "brent/brent.hpp"

namespace psi {
namespace wavekernel {

/*
Fermi-Dirac occupation number of a single energy level at energy `eps_j`
at fermi level `mu` and inverse temperature `beta`.

Returns n, in (0, 1).
 */
double n_occ(double eps_j, double mu, double beta);


/*
Derivative with respect to `mu` of the the Fermi-Dirac occupation of
a single energy level at `eps_j` at fermi level `mu` and inverse
temperature `beta`.
*/
double n_occ_prime(double eps_j, double mu, double beta);

/*
Calculate the Fermi level of an ensemble containing `N` electrons (on average)
at an inverse temperature `beta` with orbitals of energies `epsilon`.

Returns: mu, the fermi level.
*/
double calculate_mu(double N, double beta, const SharedVector& epsilon);

}
} // namespace
#endif
