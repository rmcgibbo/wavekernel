#include "fermilevel.hpp"
using namespace psi;
namespace psi {
namespace wavekernel {


double n_occ(double eps_j, double mu, double beta) {
    if (std::isinf(beta)) {
        if (eps_j < mu) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    return 1 / (1 + std::exp(beta * (eps_j - mu)));
}

double n_occ_prime(double eps_j, double mu, double beta) {
    if (std::isinf(beta)) {
        // not differentiable at eps_j == mu if beta is inf...
        return 0;
    }

    double temp = (std::exp(beta*eps_j) + std::exp(beta*mu));
    return beta * std::exp(beta*(eps_j + mu)) / (temp * temp);
}



double calculate_mu(double N, double beta, const SharedVector& epsilon) {
    class occ_functor : public brent::func_base {
    public:
        occ_functor(double N, double beta, const SharedVector& epsilon) :
            N_(N), beta_(beta), epsilon_(epsilon) {}
        double operator()(double mu) {
            double total_occupation = 0;
            for (size_t i = 0; i < epsilon_->dim(); i++) {
                total_occupation += n_occ(epsilon_->get(i), mu, beta_);
            }
            return total_occupation - N_;
        }
    private:
        const SharedVector& epsilon_;
        const double beta_;
        double N_;
    };

    /* find the largest and smallest energies to use for bracketing the search.
    */
    double min = std::numeric_limits<double>::min();
    double max = std::numeric_limits<double>::max();
    for (size_t i = 0; i < epsilon->dim(); i++) {
        double e_i = epsilon->get(i);
        if (e_i < min) {
            min = e_i;
        }
        if (e_i > max) {
            max = e_i;
        }
    }

    occ_functor occ(N, beta, epsilon);
    return brent::zero(min, max, brent::r8_epsilon(), occ);
}

}
} // namespace
