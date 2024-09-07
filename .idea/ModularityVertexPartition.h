#ifndef MODULARITYVERTEXPARTITION_H
#define MODULARITYVERTEXPARTITION_H

#include "MutableVertexPartition.h"

namespace mpl2 {

class LIBLEIDENALG_EXPORT ModularityVertexPartition : public MutableVertexPartition
{
  public:
    ModularityVertexPartition(GraphHelper* graph,
        vector<size_t> const& membership);
    ModularityVertexPartition(GraphHelper* graph);
    virtual ~ModularityVertexPartition();
    virtual ModularityVertexPartition* create(GraphHelper* graph);
    virtual ModularityVertexPartition* create(GraphHelper* graph, vector<size_t> const& membership);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();

  protected:
  private:
};

} // namespace mpl2
#endif // MODULARITYVERTEXPARTITION_H
