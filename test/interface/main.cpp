#include "leidenInterface.h"
#include "ModularityVertexPartition.h"
#include "Optimiser.h"
#include <iostream>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>

int main()
{
    // Initialize random seed
    std::srand(std::time(nullptr));
    size_t numberVertices = 500;
    size_t numberEdges = numberVertices/5;
    // Create a graph with 30 vertices
    mpl2::GraphForLeidenAlgorithm graph_(numberVertices);

    // Add edges with random weights
    for (int i = 0; i < numberVertices; ++i)
    {
        for (int j = i + 1; j < numberVertices; ++j)
        {
            if (std::rand() % numberEdges == 0)
            {                                                              // Randomly decide whether to add an edge (10% chance)
                float weight = static_cast<float>(std::rand()) / RAND_MAX; // Random weight between 0 and 1
                graph_.addEdge(i, j, weight);
            }
        }
    }

    // Print the adjacency list
    for (int i = 0; i < numberVertices; ++i)
    {
        std::cout << "Vertex " << i << " is connected to: ";
        for (const auto &neighbor : graph_.get_neighbours(i))
        {
            std::cout << neighbor << " (weight: " << graph_.getEdgeWeight(i, neighbor) << "), ";
        }
        std::cout << std::endl;
    }
    std::cout << "Finish Creating Graph" << std::endl;

    mpl2::ModularityVertexPartition part(&graph_);
    mpl2::Optimiser optimiser;
    optimiser.optimise_partition(&part);
    std::cout << "Node\tCommunity" << std::endl;
    for (int i = 0; i < graph_.numVertices(); i++)
        std::cout << i << "\t" << part.membership(i) << std::endl;
    return 0;
}