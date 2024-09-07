///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (c) 2021, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#include "leidenClusterEngine.h"

#include <queue>

#include "clusterEngine.h"
#include "db_sta/dbNetwork.hh"
#include "object.h"
#include "sta/Liberty.hh"

namespace mpl2 {

// Since init function has convert all the dbInst* to connect to the root
// cluster, leidenAlgorithm can get all the std cell, and create a hypergraph,
// there we get igraph from block and network
LeidenClusterEngine::LeidenClusterEngine(Cluster* root,
                                         Metrics* metrics,
                                         odb::dbBlock* block,
                                         sta::dbNetwork* network)
{
  // set up leidenClusterEngine for hypergraph
  root_ = root;
  design_metrics_ = metrics;
  block_ = block;
  network_ = network;
}

bool LeidenClusterEngine::isIgnoredMaster(odb::dbMaster* master)
{
  // IO corners are sometimes marked as end caps
  return master->isPad() || master->isCover() || master->isEndCap();
}

void LeidenClusterEngine::run()
{
  // setup igraph graph
  igraph_t graph;

  createGraph(graph);
  
  GraphHelper g(&graph);

  // get partition
  ModularityVertexPartition part(&g);

  // run leiden algorithm
  Optimiser o;
  o.optimise_partition(&part);

  // assign the partition to cluster
  resetCluster(g, part)

  // free graph
  igraph_destroy(&graph);
}

void LeidenClusterEngine::createGraph(igraph_t &graph)
{
  // Initialize an empty undirected graph with the number of standard cells
  igraph_empty(graph, block_->getTopModule()->getModInstCount(), IGRAPH_UNDIRECTED);

  for (odb::dbNet* net : block_->getNets()) {
    if (net->getSigType().isSupply()) {
      continue;
    }

    int driver_id = -1;
    std::set<int> loads_id;
    bool ignore = false;

    // Process all ITerms (pins) of the net
    for (odb::dbITerm* iterm : net->getITerms()) {
      odb::dbInst* inst = iterm->getInst();
      odb::dbMaster* master = inst->getMaster();

      if (isIgnoredMaster(master)) {
        ignore = true;
        break;
      }

      int vertex_id = -1;
      if (master->isBlock()) {
        vertex_id = odb::dbIntProperty::find(iterm, "vertex_id")->getValue();
      } else {
        odb::dbIntProperty* int_prop = odb::dbIntProperty::find(inst, "vertex_id");

        if (int_prop) {
          vertex_id = int_prop->getValue();
        } else {
          ignore = true;
          break;
        }
      }

      if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
        driver_id = vertex_id;
      } else {
        loads_id.insert(vertex_id);
      }
    }

    if (ignore) {
      continue;
    }

    // Process all BTerms (pins connected to block-level IOs) of the net
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      const int vertex_id = odb::dbIntProperty::find(bterm, "vertex_id")->getValue();
      if (bterm->getIoType() == odb::dbIoType::INPUT) {
        driver_id = vertex_id;
      } else {
        loads_id.insert(vertex_id);
      }
    }

    // Skip nets without a valid driver or without any loads
    if (driver_id < 0 || loads_id.empty()) {
      continue;
    }

    // Create edges between driver and loads
    for (int load_id : loads_id) {
      if (load_id != driver_id) {
        igraph_add_edge(graph, driver_id, load_id);
      }
    }
  }
}

void LeidenClusterEngine::resetCluster(GraphHelper &graph, ModularityVertexPartition &part)
{
  // transform partition to cluster
  std::vector<Cluster*> clusters;
  for (int i = 0; i < part.size(); i++) {
    Cluster* cluster = new Cluster();
    clusters.push_back(cluster);
  }

  for (int i = 0; i < part.size(); i++) {
    int clusterIndex = part[i];
    Cluster* cluster = clusters[clusterIndex];
    cluster->addVertex(i);
  }

  root_->clearClusters();
  for (Cluster* cluster : clusters) {
    root_->addCluster(cluster);
  }
}

}  // namespace mpl2
