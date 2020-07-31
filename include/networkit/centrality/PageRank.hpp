/*
 * PageRank.h
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

// networkit-format

#ifndef NETWORKIT_CENTRALITY_PAGE_RANK_HPP_
#define NETWORKIT_CENTRALITY_PAGE_RANK_HPP_

#include <limits>

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Compute PageRank as node centrality measure.
 * NOTE: There is an inconsistency in the definition in Newman's book (Ch. 7) regarding
 * directed graphs; we follow the verbal description, which requires to sum over the incoming
 * edges (as opposed to outgoing ones).
 */
class PageRank final : public Centrality {

public:
    enum Norm { L1Norm, L2Norm };

    static constexpr double DEFAULT_DAMP = 0.85;
    static constexpr double DEFAULT_TOL = 1e-8;

    /**
     * Constructs the PageRank class for the Graph @a G
     *
     * @param[in] G Graph to be processed.
     * @param[in] damp Damping factor of the PageRank algorithm.
     * @param[in] tol Error tolerance for PageRank iteration.
     * @param[in] personalization Personalization vector for PageRank iteration.
     * @param[in] mask Vector of Nodes to ignore for PageRank Calculation
     */
    PageRank(const Graph &G, double damp = DEFAULT_DAMP, double tol = DEFAULT_TOL,
             const std::vector<node> &personalization = std::vector<node>(),
             const std::vector<node> &mask = std::vector<node>());

    /**
     * Computes page rank on the graph passed in constructor.
     */
    void run() override;

    /**
     * Returns upper bound on the page rank: 1.0. This could be tighter by assuming e.g. a star
     * graph with n nodes.
     */
    double maximum() override;

    /**
     * Return the number of iterations performed by the algorithm.
     *
     * @return Number of iterations performed by the algorithm.
     */
    count numberOfIterations() const {
        assureFinished();
        return iterations;
    }

    // Maximum number of iterations allowed
    count maxIterations = std::numeric_limits<count>::max();

    // Norm used as stopping criterion
    Norm norm = Norm::L2Norm;

private:
    double damp;
    double tol;
    const std::vector<node> personalization;
    const std::vector<node> mask;
    count iterations;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_PAGE_RANK_HPP_
