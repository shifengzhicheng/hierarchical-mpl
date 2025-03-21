#include "GraphHelperForLeidenHelper.h"

#ifdef DEBUG
  using std::cerr;
  using std::endl;
#endif

vector<size_t> range(size_t n)
{
  vector<size_t> range_vec(n);
  for(size_t i = 0; i < n; i++)
    range_vec[i] = i;
  return range_vec;
}

bool orderCSize(const size_t* A, const size_t* B)
{

  if (A[1] == B[1])
  {
    if (A[2] == B[2])
      return A[0] < B[0];
    else
      return A[2] > B[2];
  }
  else
    return A[1] > B[1];
}

void shuffle(vector<size_t>& v, igraph_rng_t* rng)
{
  size_t n = v.size();
  if (n > 0)
  {
    for (size_t idx = n - 1; idx > 0; idx--)
    {
      size_t rand_idx = get_random_int(0, idx, rng);
      size_t tmp = v[idx];
      v[idx] = v[rand_idx];
      v[rand_idx] = tmp;
    }
  }
}

/****************************************************************************
  The binary Kullback-Leibler divergence.
****************************************************************************/
double KL(double q, double p)
{
  double KL = 0.0;
  if (q > 0.0 && p > 0.0)
    KL += q*log(q/p);
  if (q < 1.0 && p < 1.0)
    KL += (1.0-q)*log((1.0-q)/(1.0-p));
  return KL;
}

double KLL(double q, double p)
{
  double KL = 0.0;
  if (q > 0.0 && p > 0.0)
    KL += q*log(q/p);
  if (q < 1.0 && p < 1.0)
    KL += (1.0-q)*log((1.0-q)/(1.0-p));
  if (q < p)
    KL *= -1;
  return KL;
}

GraphHelperForLeiden::GraphHelperForLeiden(igraph_t* graph,
  vector<double> const& edge_weights,
  vector<double> const& node_sizes,
  vector<double> const& node_self_weights, int correct_self_loops)
{
  this->_graph = graph;
  this->_remove_graph = false;

  if (edge_weights.size() != this->ecount())
    throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
  this->_edge_weights = edge_weights;
  this->_is_weighted = true;

  if (node_sizes.size() != this->vcount())
    throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
  this->_node_sizes = node_sizes;

  if (node_self_weights.size() != this->vcount())
    throw Exception("Node self weights vector inconsistent length with the vertex count of the graph.");
  this->_node_self_weights = node_self_weights;

  this->_correct_self_loops = correct_self_loops;
  igraph_vector_int_init(&this->_temp_igraph_vector, this->vcount());
  this->init_admin();
}

GraphHelperForLeiden::GraphHelperForLeiden(igraph_t* graph,
  vector<double> const& edge_weights,
  vector<double> const& node_sizes,
  vector<double> const& node_self_weights)
{
  this->_graph = graph;
  this->_remove_graph = false;

  if (edge_weights.size() != this->ecount())
    throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
  this->_edge_weights = edge_weights;
  this->_is_weighted = true;

  if (node_sizes.size() != this->vcount())
    throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
  this->_node_sizes = node_sizes;

  this->_correct_self_loops = this->has_self_loops();

  this->_node_self_weights = node_self_weights;
  igraph_vector_int_init(&this->_temp_igraph_vector, this->vcount());
  this->init_admin();
}

GraphHelperForLeiden::GraphHelperForLeiden(igraph_t* graph,
  vector<double> const& edge_weights,
  vector<double> const& node_sizes, int correct_self_loops)
{
  this->_graph = graph;
  this->_remove_graph = false;

  if (edge_weights.size() != this->ecount())
    throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
  this->_edge_weights = edge_weights;
  this->_is_weighted = true;

  if (node_sizes.size() != this->vcount())
    throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
  this->_node_sizes = node_sizes;

  this->_correct_self_loops = correct_self_loops;
  igraph_vector_int_init(&this->_temp_igraph_vector, this->vcount());
  this->init_admin();
  this->set_self_weights();
}

GraphHelperForLeiden::GraphHelperForLeiden(igraph_t* graph,
  vector<double> const& edge_weights,
  vector<double> const& node_sizes)
{
  this->_graph = graph;
  this->_remove_graph = false;
  if (edge_weights.size() != this->ecount())
    throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
  this->_edge_weights = edge_weights;
  this->_is_weighted = true;

  if (node_sizes.size() != this->vcount())
    throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
  this->_node_sizes = node_sizes;

  this->_correct_self_loops = this->has_self_loops();

  igraph_vector_int_init(&this->_temp_igraph_vector, this->vcount());
  this->init_admin();
  this->set_self_weights();
}

GraphHelperForLeiden* GraphHelperForLeiden::GraphHelperForLeidenFromEdgeWeights(igraph_t* graph, vector<double> const& edge_weights, int correct_self_loops)
{
  GraphHelperForLeiden* g = new GraphHelperForLeiden(graph, correct_self_loops);

  if (edge_weights.size() != g->ecount())
    throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
  g->_edge_weights = edge_weights;
  g->_is_weighted = true;
  g->set_default_node_size();
  igraph_vector_int_init(&g->_temp_igraph_vector, g->vcount());
  g->init_admin();
  g->set_self_weights();

  return g;
}

GraphHelperForLeiden* GraphHelperForLeiden::GraphHelperForLeidenFromEdgeWeights(igraph_t* graph, vector<double> const& edge_weights)
{
  GraphHelperForLeiden* g = new GraphHelperForLeiden(graph);

  if (edge_weights.size() != g->ecount())
    throw Exception("Edge weights vector inconsistent length with the edge count of the graph.");
  g->_edge_weights = edge_weights;
  g->_is_weighted = true;
  g->set_default_node_size();
  igraph_vector_int_init(&g->_temp_igraph_vector, g->vcount());
  g->init_admin();
  g->set_self_weights();

  return g;
}

GraphHelperForLeiden* GraphHelperForLeiden::GraphHelperForLeidenFromNodeSizes(igraph_t* graph, vector<double> const& node_sizes, int correct_self_loops)
{
  GraphHelperForLeiden* g = new GraphHelperForLeiden(graph, correct_self_loops);

  if (node_sizes.size() != g->vcount())
    throw Exception("Node size vector inconsistent length with the vertex count of the graph.");
  g->_node_sizes = node_sizes;

  g->set_default_edge_weight();
  g->_is_weighted = false;
  igraph_vector_int_init(&g->_temp_igraph_vector, g->vcount());
  g->init_admin();
  g->set_self_weights();
  
  return g;
}

GraphHelperForLeiden* GraphHelperForLeiden::GraphHelperForLeidenFromNodeSizes(igraph_t* graph, vector<double> const& node_sizes)
{
  GraphHelperForLeiden* g = new GraphHelperForLeiden(graph);

  g->_graph = graph;
  g->_remove_graph = false;
  g->set_defaults();
  g->_is_weighted = false;

  if (node_sizes.size() != g->vcount())
    throw Exception("Node size vector inconsistent length with the vertex count of the graph.");

  g->_node_sizes = node_sizes;

  g->_correct_self_loops = g->has_self_loops();

  igraph_vector_int_init(&g->_temp_igraph_vector, g->vcount());
  g->init_admin();
  g->set_self_weights();

  return g;
}

GraphHelperForLeiden::GraphHelperForLeiden(igraph_t* graph, int correct_self_loops)
{
  this->_graph = graph;
  this->_remove_graph = false;
  this->_correct_self_loops = correct_self_loops;
  this->set_defaults();
  this->_is_weighted = false;
  igraph_vector_int_init(&this->_temp_igraph_vector, this->vcount());
  this->init_admin();
  this->set_self_weights();
}

GraphHelperForLeiden::GraphHelperForLeiden(igraph_t* graph)
{
  this->_graph = graph;
  this->_remove_graph = false;
  this->set_defaults();
  this->_is_weighted = false;

  this->_correct_self_loops = this->has_self_loops();

  igraph_vector_int_init(&this->_temp_igraph_vector, this->vcount());
  this->init_admin();
  this->set_self_weights();
}

GraphHelperForLeiden::~GraphHelperForLeiden()
{
  if (this->_remove_graph)
  {
    igraph_destroy((igraph_t*)this->_graph);
    delete this->_graph;
  }
  igraph_vector_int_destroy(&this->_temp_igraph_vector);
}

int GraphHelperForLeiden::has_self_loops()
{
  igraph_bool_t has_self_loops;  
  igraph_has_loop(this->_graph, &has_self_loops);
  return has_self_loops;
}

double GraphHelperForLeiden::possible_edges()
{
  return this->possible_edges(this->vcount());
}

double GraphHelperForLeiden::possible_edges(double n)
{
  double possible_edges = n*(n-1);
  if (!this->is_directed())
    possible_edges /= 2;
  if (this->correct_self_loops())
    possible_edges += n;

  return possible_edges;
}

void GraphHelperForLeiden::set_defaults()
{
  this->set_default_edge_weight();
  this->set_default_node_size();
}

void GraphHelperForLeiden::set_default_edge_weight()
{
  size_t m = this->ecount();

  // Set default edge weight of 1.0
  this->_edge_weights.clear(); this->_edge_weights.resize(m);
  fill(this->_edge_weights.begin(), this->_edge_weights.end(), 1.0);
  this->_is_weighted = false;
}

void GraphHelperForLeiden::set_default_node_size()
{
  size_t n = this->vcount();

  // Set default node size of 1
  this->_node_sizes.clear(); this->_node_sizes.resize(n);
  fill(this->_node_sizes.begin(), this->_node_sizes.end(), 1);
}

void GraphHelperForLeiden::set_self_weights()
{
  size_t n = this->vcount();

  // Set default self_weights of the total weight of any possible self-loops
  this->_node_self_weights.clear(); this->_node_self_weights.resize(n);
  for (size_t v = 0; v < n; v++)
  {
    #ifdef DEBUG
      cerr << "\t" << "Size node " << v << ": " << this->node_size(v) << endl;
    #endif
    double self_weight = 0.0;
    // There should be only one self loop
    igraph_integer_t eid;
    // Get edge id for self loop
    igraph_get_eid(this->_graph, &eid, v, v, this->is_directed(), false);
    if (eid >= 0)
      self_weight = this->edge_weight(eid);

    this->_node_self_weights[v] = self_weight;
    #ifdef DEBUG
      cerr << "\t" << "Self weight node " << v << ": " << self_weight << endl;
    #endif
  }
}

void GraphHelperForLeiden::init_admin()
{

  size_t m = this->ecount();
  size_t n = this->vcount();
  this->_is_directed = igraph_is_directed(this->_graph);

  this->_strength_in.clear();
  this->_strength_in.resize(n, 0.0);

  this->_degree_in.clear();
  this->_degree_in.resize(n, 0.0);

  if (this->_is_directed) {
    this->_strength_out.clear();
    this->_strength_out.resize(n, 0.0);

    this->_degree_out.clear();
    this->_degree_out.resize(n, 0);

    this->_degree_all.clear();
    this->_degree_all.resize(n, 0);
  }

  // Determine total weight in the graph.
  this->_total_weight = 0.0;
  for (size_t e = 0; e < m; e++) {
    double w = this->edge_weight(e);
    this->_total_weight += w;

    size_t from, to;
    this->edge(e, from, to);

    if (this->is_directed()) {
      this->_strength_in[to] += w;
      this->_strength_out[from] += w;

      this->_degree_in[to]++;
      this->_degree_out[from]++;
      this->_degree_all[to]++;
      this->_degree_all[from]++;
    } else {
      // we only compute strength_in and degree_in for undirected graphs
      this->_strength_in[to] += w;
      this->_strength_in[from] += w;

      // recall that igraph ignores the mode for undirected graphs
      this->_degree_in[to]++;
      this->_degree_in[from]++;
    }
  }

  // Make sure to multiply by 2 for undirected graphs
  //if (!this->is_directed())
  //  this->_total_weight *= 2.0;

  this->_total_size = 0;
  for (size_t v = 0; v < n; v++)
    this->_total_size += this->node_size(v);

  // Calculate density;
  double w = this->total_weight();
  double n_size = this->total_size();

  // For now we default to not correcting self loops.
  // this->_correct_self_loops = false; (remove this as this is set in the constructor)

  double normalise = 0.0;
  if (this->_correct_self_loops)
    normalise = n_size*n_size;
  else
    normalise = n_size*(n_size - 1);

  if (this->is_directed())
    this->_density = w/normalise;
  else
    this->_density = 2*w/normalise;

  this->_current_node_cache_neigh_edges_from = n + 1;
  this->_current_node_cache_neigh_edges_to = n + 1;
  this->_current_node_cache_neigh_edges_all = n + 1;

  this->_current_node_cache_neigh_from = n + 1;
  this->_current_node_cache_neigh_to = n + 1;
  this->_current_node_cache_neigh_all = n + 1;
}

void GraphHelperForLeiden::cache_neighbour_edges(size_t v, igraph_neimode_t mode)
{
  #ifdef DEBUG
    cerr << "void GraphHelperForLeiden::cache_neighbour_edges(" << v << ", " << mode << ");" << endl;
  #endif
  size_t degree = this->degree(v, mode);
  #ifdef DEBUG
    cerr << "Degree: " << degree << endl;
  #endif

  igraph_vector_int_t *incident_edges = &this->_temp_igraph_vector;
  igraph_incident(this->_graph, incident_edges, v, mode);

  vector<size_t>* _cached_neigh_edges = NULL;
  switch (mode)
  {
    case IGRAPH_IN:
      this->_current_node_cache_neigh_edges_from = v;
      _cached_neigh_edges = &(this->_cached_neigh_edges_from);
      break;
    case IGRAPH_OUT:
      this->_current_node_cache_neigh_edges_to = v;
      _cached_neigh_edges = &(this->_cached_neigh_edges_to);
      break;
    case IGRAPH_ALL:
      this->_current_node_cache_neigh_edges_all = v;
      _cached_neigh_edges = &(this->_cached_neigh_edges_all);
      break;
  }
  _cached_neigh_edges->assign(igraph_vector_int_get_ptr(incident_edges, 0),
                              igraph_vector_int_get_ptr(incident_edges, degree));
  #ifdef DEBUG
    cerr << "Number of edges: " << _cached_neigh_edges->size() << endl;
  #endif

  #ifdef DEBUG
    cerr << "exit void GraphHelperForLeiden::cache_neighbour_edges(" << v << ", " << mode << ");" << endl;
  #endif
}

vector<size_t> const& GraphHelperForLeiden::get_neighbour_edges(size_t v, igraph_neimode_t mode)
{
  if (!this->is_directed())
    mode = IGRAPH_ALL; // igraph ignores mode for undirected graphs

  switch (mode)
  {
    case IGRAPH_IN:
      if (this->_current_node_cache_neigh_edges_from != v)
      {
        cache_neighbour_edges(v, mode);
        this->_current_node_cache_neigh_edges_from = v;
      }
      return this->_cached_neigh_edges_from;
    case IGRAPH_OUT:
      if (this->_current_node_cache_neigh_edges_to != v)
      {
        cache_neighbour_edges(v, mode);
        this->_current_node_cache_neigh_edges_to = v;
      }
      return this->_cached_neigh_edges_to;
    case IGRAPH_ALL:
      if (this->_current_node_cache_neigh_edges_all != v)
      {
        cache_neighbour_edges(v, mode);
        this->_current_node_cache_neigh_edges_all = v;
      }
      return this->_cached_neigh_edges_all;
  }
  throw Exception("Incorrect model for getting neighbour edges.");
}

void GraphHelperForLeiden::cache_neighbours(size_t v, igraph_neimode_t mode)
{
  #ifdef DEBUG
    cerr << "void GraphHelperForLeiden::cache_neighbours(" << v << ", " << mode << ");" << endl;
  #endif
  size_t degree = this->degree(v, mode);
  #ifdef DEBUG
    cerr << "Degree: " << degree << endl;
  #endif

  igraph_vector_int_t *neighbours = &this->_temp_igraph_vector;
  igraph_neighbors(this->_graph, neighbours, v, mode);

  vector<size_t>* _cached_neighs = NULL;
  switch (mode)
  {
    case IGRAPH_IN:
      this->_current_node_cache_neigh_from = v;
      _cached_neighs = &(this->_cached_neighs_from);
      break;
    case IGRAPH_OUT:
      this->_current_node_cache_neigh_to = v;
      _cached_neighs = &(this->_cached_neighs_to);
      break;
    case IGRAPH_ALL:
      this->_current_node_cache_neigh_all = v;
      _cached_neighs = &(this->_cached_neighs_all);
      break;
  }
  _cached_neighs->assign(igraph_vector_int_get_ptr(neighbours, 0),
                         igraph_vector_int_get_ptr(neighbours, degree));

  #ifdef DEBUG
    cerr << "Number of edges: " << _cached_neighs->size() << endl;
  #endif

  #ifdef DEBUG
    cerr << "exit void GraphHelperForLeiden::cache_neighbours(" << v << ", " << mode << ");" << endl;
  #endif
}

vector< size_t > const& GraphHelperForLeiden::get_neighbours(size_t v, igraph_neimode_t mode)
{
  if (!this->is_directed())
    mode = IGRAPH_ALL; // igraph ignores mode for undirected graphs

  switch (mode)
  {
    case IGRAPH_IN:
      if (this->_current_node_cache_neigh_from != v)
      {
        cache_neighbours(v, mode);
        this -> _current_node_cache_neigh_from = v;
      }
      #ifdef DEBUG
        cerr << "Returning " << this->_cached_neighs_from.size() << " incoming neighbours" << endl;
      #endif
      return this->_cached_neighs_from;
    case IGRAPH_OUT:
      if (this->_current_node_cache_neigh_to != v)
      {
        cache_neighbours(v, mode);
        this -> _current_node_cache_neigh_to = v;
      }
      #ifdef DEBUG
        cerr << "Returning " << this->_cached_neighs_to.size() << " incoming neighbours" << endl;
      #endif
      return this->_cached_neighs_to;
    case IGRAPH_ALL:
      if (this->_current_node_cache_neigh_all != v)
      {
        cache_neighbours(v, mode);
        this->_current_node_cache_neigh_all = v;
      }
      #ifdef DEBUG
        cerr << "Returning " << this->_cached_neighs_all.size() << " incoming neighbours" << endl;
      #endif
      return this->_cached_neighs_all;
  }
  throw Exception("Invalid mode for getting neighbours.");
}

/********************************************************************************
 * This should return a random neighbour in O(1)
 ********************************************************************************/
size_t GraphHelperForLeiden::get_random_neighbour(size_t v, igraph_neimode_t mode, igraph_rng_t* rng)
{
  size_t node=v;
  size_t rand_neigh = -1;

  if (this->degree(v, mode) <= 0)
    throw Exception("Cannot select a random neighbour for an isolated node.");

  if (this->is_directed() && mode != IGRAPH_ALL)
  {
    if (mode == IGRAPH_OUT)
    {
      // Get indices of where neighbours are
      size_t cum_degree_this_node = (size_t) VECTOR(this->_graph->os)[node];
      size_t cum_degree_next_node = (size_t) VECTOR(this->_graph->os)[node+1];
      // Get a random index from them
      size_t rand_neigh_idx = get_random_int(cum_degree_this_node, cum_degree_next_node - 1, rng);
      // Return the neighbour at that index
      #ifdef DEBUG
        cerr << "Degree: " << this->degree(node, mode) << " diff in cumulative: " << cum_degree_next_node - cum_degree_this_node << endl;
      #endif
      rand_neigh = VECTOR(this->_graph->to)[ (size_t)VECTOR(this->_graph->oi)[rand_neigh_idx] ];
    }
    else if (mode == IGRAPH_IN)
    {
      // Get indices of where neighbours are
      size_t cum_degree_this_node = (size_t) VECTOR(this->_graph->is)[node];
      size_t cum_degree_next_node = (size_t) VECTOR(this->_graph->is)[node+1];
      // Get a random index from them
      size_t rand_neigh_idx = get_random_int(cum_degree_this_node, cum_degree_next_node - 1, rng);
      #ifdef DEBUG
        cerr << "Degree: " << this->degree(node, mode) << " diff in cumulative: " << cum_degree_next_node - cum_degree_this_node << endl;
      #endif
      // Return the neighbour at that index
      rand_neigh = VECTOR(this->_graph->from)[ (size_t)VECTOR(this->_graph->ii)[rand_neigh_idx] ];
    }
  }
  else
  {
    // both in- and out- neighbors in a directed graph.
    size_t cum_outdegree_this_node = (size_t)VECTOR(this->_graph->os)[node];
    size_t cum_indegree_this_node  = (size_t)VECTOR(this->_graph->is)[node];

    size_t cum_outdegree_next_node = (size_t)VECTOR(this->_graph->os)[node+1];
    size_t cum_indegree_next_node  = (size_t)VECTOR(this->_graph->is)[node+1];

    size_t total_outdegree = cum_outdegree_next_node - cum_outdegree_this_node;
    size_t total_indegree = cum_indegree_next_node - cum_indegree_this_node;

    size_t rand_idx = get_random_int(0, total_outdegree + total_indegree - 1, rng);

    #ifdef DEBUG
      cerr << "Degree: " << this->degree(node, mode) << " diff in cumulative: " << total_outdegree + total_indegree << endl;
    #endif
    // From among in or out neighbours?
    if (rand_idx < total_outdegree)
    { // From among outgoing neighbours
      size_t rand_neigh_idx = cum_outdegree_this_node + rand_idx;
      rand_neigh = VECTOR(this->_graph->to)[ (size_t)VECTOR(this->_graph->oi)[rand_neigh_idx] ];
    }
    else
    { // From among incoming neighbours
      size_t rand_neigh_idx = cum_indegree_this_node + rand_idx - total_outdegree;
      rand_neigh = VECTOR(this->_graph->from)[ (size_t)VECTOR(this->_graph->ii)[rand_neigh_idx] ];
    }
  }

  return rand_neigh;
}

/****************************************************************************
  Creates a graph with communities as node and links as weights between
  communities.

  The weight of the edges in the new graph is simply the sum of the weight
  of the edges between the communities. The self weight of a node (i.e. the
  weight of its self loop) is the internal weight of a community. The size
  of a node in the new graph is simply the size of the community in the old
  graph.
*****************************************************************************/
GraphHelperForLeiden* GraphHelperForLeiden::collapse_graph(MutableVertexPartition* partition)
{
  size_t n_collapsed = partition->n_communities();
  vector<vector<size_t> > community_memberships = partition->get_communities();

  vector<double> collapsed_weights;
  double total_collapsed_weight = 0.0;

  vector<double> edge_weight_to_community(n_collapsed, 0.0);
  vector<bool> neighbour_comm_added(n_collapsed, false);

  // collapsed edges for new graph
  igraph_vector_int_t edges;
  igraph_vector_int_init(&edges, 0);

  for (size_t v_comm = 0; v_comm < n_collapsed; v_comm++) {
    vector<size_t> neighbour_communities;
    for (size_t v : community_memberships[v_comm]) {
        for (size_t e : this->get_neighbour_edges(v, IGRAPH_OUT)) {
            size_t from, to;
            this->edge(e, from, to);

            if ((size_t) from != v) {
                // need to skip because IGRAPH_OUT is ignored for undirected graphs
                continue;
            }

            size_t u_comm = partition->membership(to);

            double w = this->edge_weight(e);
            // Self loops appear twice here if the graph is undirected, so divide by 2.0 in that case.
            if (from == to && !this->is_directed())
                w /= 2.0;

            if (!neighbour_comm_added[u_comm]) {
                neighbour_comm_added[u_comm] = true;
                neighbour_communities.push_back(u_comm);
            }
            edge_weight_to_community[u_comm] += w;
        }
    }

    for (size_t u_comm : neighbour_communities) {
        igraph_vector_int_push_back(&edges, v_comm);
        igraph_vector_int_push_back(&edges, u_comm);
        collapsed_weights.push_back(edge_weight_to_community[u_comm]);
        total_collapsed_weight += edge_weight_to_community[u_comm];

        // reset edge_weight_to_community to all 0.0 and neighbour_comm_added to all false
        edge_weight_to_community[u_comm] = 0.0;
        neighbour_comm_added[u_comm] = false;
    }
  }

  // Create graph based on edges
  igraph_t* graph = new igraph_t();
  igraph_create(graph, &edges, n_collapsed, this->is_directed());
  igraph_vector_int_destroy(&edges);

  if ((size_t) igraph_vcount(graph) != partition->n_communities())
    throw Exception("Something went wrong with collapsing the graph.");

  // Calculate new node sizes
  vector<double> csizes(n_collapsed, 0);
  for (size_t c = 0; c < partition->n_communities(); c++)
    csizes[c] = partition->csize(c);

  GraphHelperForLeiden* G = new GraphHelperForLeiden(graph, collapsed_weights, csizes, this->_correct_self_loops);
  G->_remove_graph = true;
  return G;
}
