/*
 * PageRank.cpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

// networkit-format

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/PageRank.hpp>

namespace NetworKit {

PageRank::PageRank(const Graph &G, double damp, double tol,
                   const std::vector<node> &personalization)
    : Centrality(G, true), damp(damp), tol(tol), personalization(personalization) {}

void PageRank::run() {
    Aux::SignalHandler handler;
    const auto n = G.numberOfNodes();
    const auto z = G.upperNodeIdBound();

    auto teleportProb = (1.0 - damp) / static_cast<double>(n);
    scoreData.resize(z, 1.0 / static_cast<double>(n));
    std::vector<double> pr = scoreData;

    // Verify the personalization vector
    std::vector<node>::const_iterator it;
    auto p = personalization.size();
    for (it = personalization.begin(); it != personalization.end(); ++it) {
        if (!G.hasNode(*it)) {
            --p;
            WARN("Node ", *it, " is in the personalization vector but doesn't exist in the Graph");
        }
    }

    std::vector<double> deg(z, 0.0);
    G.parallelForNodes([&](const node u) { deg[u] = static_cast<double>(G.weightedDegree(u)); });

    iterations = 0;

    bool isPersonalized = (p > 0);
    std::vector<bool> inPersonalization;
    std::function<void(node)> teleport;
    if (isPersonalized) {
        inPersonalization.resize(z, false);
        std::vector<node>::const_iterator it;
        for (it = personalization.begin(); it != personalization.end(); ++it) {
            if (G.hasNode(*it)) {
                inPersonalization[*it] = true;
            }
        }

        teleportProb = (1.0 - damp) / static_cast<double>(p);
        teleport = [&](const node u) {
            if (inPersonalization[u]) {
                pr[u] += teleportProb;
            }
        };
    } else {
        teleport = [&](const node u) { pr[u] += teleportProb; };
    }

    auto sumL1Norm = [&](const node u) { return std::abs(scoreData[u] - pr[u]); };

    auto sumL2Norm = [&](const node u) {
        const auto d = scoreData[u] - pr[u];
        return d * d;
    };

    auto converged([&]() {
        if (iterations >= maxIterations) {
            return true;
        }

        if (norm == Norm::L2Norm) {
            return std::sqrt(G.parallelSumForNodes(sumL2Norm)) <= tol;
        }

        return G.parallelSumForNodes(sumL1Norm) <= tol;
    });

    bool isConverged = false;
    do {
        handler.assureRunning();
        G.balancedParallelForNodes([&](const node u) {
            pr[u] = 0.0;
            G.forInEdgesOf(u, [&](const node u, const node v, const edgeweight w) {
                // note: inconsistency in definition in Newman's book (Ch. 7) regarding directed
                // graphs we follow the verbal description, which requires to sum over the incoming
                // edges
                pr[u] += scoreData[v] * w / deg[v];
            });
            pr[u] *= damp;
            teleport(u);
        });

        ++iterations;
        isConverged = converged();
        std::swap(pr, scoreData);
    } while (!isConverged);

    handler.assureRunning();

    const auto sum = G.parallelSumForNodes([&](const node u) { return scoreData[u]; });

    // make sure scoreData sums up to 1
    assert(!Aux::NumericTools::equal(sum, 0.0, 1e-15));
    G.parallelForNodes([&](const node u) { scoreData[u] /= sum; });

    hasRun = true;
}

double PageRank::maximum() {
    return 1.0; // upper bound, could be tighter by assuming e.g. a star graph with n nodes
}

} /* namespace NetworKit */
