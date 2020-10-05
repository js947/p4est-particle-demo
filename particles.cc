#include <array>
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <p4est_vtk.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_iterate.h>

using point_t = std::array<double,P4EST_DIM>;
using data_t = std::vector<point_t>;

data_t* particles(p4est_quadrant_t *q) {
  return (data_t*) q->p.user_data;
}

bool point_in(point_t x, point_t xmin, point_t xmax) {
  static_assert(P4EST_DIM == 2);
  return (x[0] >= xmin[0]) &&
         (x[0] <  xmax[0]) &&
         (x[1] >= xmin[1]) &&
         (x[1] <  xmax[1]);
}

struct read_initial_ctx_t {
  std::string filename;
};

std::pair<point_t, point_t> quad_boundaries(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *q) {
  point_t xmin, xmax;
  p4est_qcoord_t sz = P4EST_QUADRANT_LEN (q->level);
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x, q->y, xmin.data());
  //printf("quad boundaries xmin[%zu] %f %f\n", xmin.size(), xmin[0], xmin[1]);
  p4est_qcoord_to_vertex (p4est->connectivity, which_tree, q->x+sz, q->y+sz, xmax.data());
  //printf("quad boundaries xmax[%zu] %f %f\n", xmax.size(), xmax[0], xmax[1]);
  return {xmin, xmax};
}

void read_initial(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *q) {
  read_initial_ctx_t *ctx = (read_initial_ctx_t*) p4est->user_pointer;
  if (!q->p.user_data)
    q->p.user_data = (void*) new data_t();
  data_t *particles = (data_t*) q->p.user_data;

  std::pair<point_t, point_t> xs;
  //auto [xmin, xmax] = quad_boundaries(p4est, which_tree, q);
  xs = quad_boundaries(p4est, which_tree, q);
  point_t xmin = xs.first;
  point_t xmax = xs.second;

  std::ifstream input(ctx->filename);
  std::string str;
  while (std::getline(input, str)) {
    point_t x;
    size_t ir = 0;
    x[0] = std::stod(str, &ir);
    x[1] = std::stod(str.substr(ir));

    if (point_in(x, xmin, xmax))
      particles->push_back(x);
  }
  printf("read x [%f %f] y [%f %f] => %zu\n", xmin[0], xmax[0], xmin[1], xmax[1], particles->size());
}

void replace_quads(p4est_t * p4est, p4est_topidx_t which_tree,
                     int num_outgoing, p4est_quadrant_t * outgoing[],
                     int num_incoming, p4est_quadrant_t * incoming[]) {
  for (int j = 0; j < num_incoming; ++j) {
    auto [xmin, xmax] = quad_boundaries(p4est, which_tree, incoming[j]);
    if (!particles(incoming[j]))
      incoming[j]->p.user_data = (void*) new data_t();
    assert(particles(incoming[j]));
    for (int i = 0; i < num_outgoing; ++i) {
      assert(particles(outgoing[i]));
      for (auto x : *particles(outgoing[i]))
        if (point_in(x, xmin, xmax))
          particles(incoming[j])->push_back(x);
    }
  }

  size_t outgoing_total = 0;
  for (int i = 0; i < num_outgoing; ++i) {
    outgoing_total += particles(outgoing[i])->size();
    auto [xmin, xmax] = quad_boundaries(p4est, which_tree, outgoing[i]);
    printf("replace_quads old[%d] x [%f %f] y [%f %f] => %zu\n", i,
        xmin[0], xmax[0],
        xmin[1], xmax[1],
        particles(outgoing[i])->size());
  }
  size_t incoming_total = 0;
  for (int i = 0; i < num_incoming; ++i) {
    incoming_total += particles(incoming[i])->size();
    auto [xmin, xmax] = quad_boundaries(p4est, which_tree, incoming[i]);
    printf("replace_quads new[%d] x [%f %f] y [%f %f] => %zu\n", i,
        xmin[0], xmax[0],
        xmin[1], xmax[1],
        particles(incoming[i])->size());
  }
  printf("replace_quads total %d %zu -> %d %zu\n",
      num_outgoing, outgoing_total,
      num_incoming, incoming_total);

  assert(outgoing_total == incoming_total);

  for (int i = 0; i < num_outgoing; ++i) {
    particles(outgoing[i])->clear();
  }
}

int refine_flag(p4est_t *, p4est_topidx_t, p4est_quadrant_t * q) {
  return particles(q)->size() > 8;
}

void count_particles(p4est_iter_volume_info_t * info, void *user_data) {
  sc_array_t         *n_particles = (sc_array_t *) user_data;
  p4est_t            *p4est = info->p4est;
  p4est_quadrant_t   *q = info->quad;
  p4est_topidx_t      which_tree = info->treeid;
  p4est_locidx_t      local_id = info->quadid; 
  p4est_tree_t       *tree;
  p4est_locidx_t      arrayoffset;

  tree = p4est_tree_array_index (p4est->trees, which_tree);
  local_id += tree->quadrants_offset;   /* now the id is relative to the MPI process */
  arrayoffset = P4EST_CHILDREN * local_id;      /* each local quadrant has 2^d (P4EST_CHILDREN) values in u_interp */

  for (size_t i = 0; i < P4EST_CHILDREN; i++)
    *((double *) sc_array_index (n_particles, arrayoffset + i)) = particles(q)->size();
}

void write_output(p4est_t *p4est, char const *filename) {
  p4est_vtk_context_t *context = p4est_vtk_context_new(p4est, filename);
  p4est_vtk_context_set_scale (context, 0.99);
  context = p4est_vtk_write_header (context);
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing vtk header");

  sc_array_t *n_particles = sc_array_new_size (sizeof (double), p4est->local_num_quadrants * P4EST_CHILDREN);

  p4est_iterate(p4est, NULL,   /* we don't need any ghost quadrants for this loop */
               (void *) n_particles,     /* pass in n_particles so that we can fill it */
               count_particles,    /* callback function that counts particles in a cell */
               NULL,          /* there is no callback for the faces between quadrants */
#ifdef P4_TO_P8
               NULL,          /* there is no callback for the edges between quadrants */
#endif
               NULL);         /* there is no callback for the corners between quadrants */

  context = p4est_vtk_write_cell_dataf (context, 0, 1,  /* do write the refinement level of each quadrant */
                                        1,      /* do write the mpi process id of each quadrant */
                                        0,      /* do not wrap the mpi rank (if this were > 0, the modulus of the rank relative to this number would be written instead of the rank) */
                                        0,      /* there is no custom cell scalar data. */
                                        0,      /* there is no custom cell vector data. */
                                        context);       /* mark the end of the variable cell data. */

  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing cell info");

  /* write one scalar field: the solution value */
  context = p4est_vtk_write_point_dataf (context, 1, 0, /* write no vector fields */
                                         "particles", n_particles, context);        /* mark the end of the variable cell data. */
  SC_CHECK_ABORT (context != NULL, P4EST_STRING "_vtk: Error writing cell data");

  int retval = p4est_vtk_write_footer (context);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");

  sc_array_destroy(n_particles);
}

int main(int argc, char **argv) {
  SC_CHECK_MPI(sc_MPI_Init(&argc, &argv));
  sc_MPI_Comm mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 1, 1, NULL, SC_LP_ESSENTIAL);
  p4est_init (NULL, SC_LP_PRODUCTION);

  p4est_connectivity_t *conn = p4est_connectivity_new_periodic();

  read_initial_ctx_t read_initial_ctx{(argc > 1) ? argv[1] : "input_data_1k"};

  p4est_t *p4est = p4est_new_ext (mpicomm, conn, 0, 3, 1, sizeof(data_t), read_initial, (void *)&read_initial_ctx);

  p4est_refine_ext (p4est, 1, 5, refine_flag, NULL, replace_quads);
  p4est_balance_ext (p4est, P4EST_CONNECT_CORNER, NULL, replace_quads);
  write_output(p4est, (argc > 2) ? argv[2] : "particles.vtk");

  p4est_connectivity_destroy (conn);

  SC_CHECK_MPI(sc_MPI_Finalize());
}
