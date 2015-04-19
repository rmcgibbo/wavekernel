#ifndef _psi_wavekernel_utils_h
#define _psi_wavekernel_utils_h

#include "boost/random.hpp"

using namespace boost;
using namespace psi;
namespace psi{ namespace wavekernel {


/*
    Generate n_samples random indices for points on the molecular grid,
    broken into blocks.

    The return value is a list of vectors of length equal to the number of blocks,
    where each element is a vector of the selected indices in that block, (in
    the internal indexing scheme of the block, not the global indexing)
*/
std::vector<std::vector<size_t> >
blocks_rand_subset(const std::vector<shared_ptr<BlockOPoints> >& blocks, size_t n_samples)
{
    size_t npoints = 0;  // total number of points in the grid
    std::vector<size_t> block_start;  // starting index of each block
    for (size_t Q = 0; Q < blocks.size(); Q++) {
        block_start.push_back(npoints);
        npoints += blocks[Q]->npoints();
    }
    block_start.push_back(npoints);

    boost::random::mt19937 gen(time(0));
    boost::random::uniform_int_distribution<size_t> dist(0, npoints-1);
    std::vector<size_t> samples;
    for (size_t i = 0; i < n_samples; i++) {
        samples.push_back(dist(gen));
    }
    std::sort(samples.begin(), samples.end());
    std::vector<std::vector<size_t> > block_indices(blocks.size());

    size_t j = 0;
    for (size_t i = 0; i < n_samples; i++) {
        if (samples[i] >= block_start[j+1]) {
            j += 1;
        }
        block_indices[j].push_back(samples[i] - block_start[j]);
    }

    return block_indices;
}

}}
#endif
