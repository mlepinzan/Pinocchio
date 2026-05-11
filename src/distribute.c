/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"

#define BUFLEN 50000

// #define DEBUG

static product_data *comm_buffer;
#ifndef CLASSIC_FRAGMENTATION
static unsigned long long frag_offset;
#endif

int intersection(int *, int *, int *);
int send_data(int *, int);
int recv_data(int *, int);
int keep_data(int *, int *);
unsigned int fft_space_index(unsigned int, int *);
unsigned int subbox_space_index(unsigned int, int *);
int distribute_alltoall(void);
#ifndef CLASSIC_FRAGMENTATION
int get_distmap_bit(unsigned int *, unsigned int);
void set_distmap_bit(unsigned int *, unsigned int, int);
void build_distmap(unsigned int *, int *);
void update_distmap(unsigned int *, int *);
unsigned int distmap_length(unsigned int);
#endif

typedef struct
{
  unsigned long long pos;
  product_data prod;
} dist_data;

static int get_task_count_size(size_t *ntasks)
{
  if (NTasks <= 0)
  {
    printf("ERROR on task %d: invalid MPI task count %d in distribute\n", ThisTask, NTasks);
    fflush(stdout);
    return 1;
  }

  *ntasks = (size_t)NTasks;
  return 0;
}

#ifdef DEBUG
FILE *DBGFD;
#endif

int distribute(void)
{
  /* Distributes products from fft-space to sub-volumes */

  int my_fft_box[6], my_subbox[6];
  int log_ntask, bit, receiver, sender;

  /* this defines the box that the task possesses in the FFT space */
  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  /* this defines the box that the task possesses in the subbox space
     (the starting coordinate may be negative) */
  my_subbox[0] = subbox.stabl[_x_];
  my_subbox[1] = subbox.stabl[_y_];
  my_subbox[2] = subbox.stabl[_z_];
  my_subbox[3] = subbox.Lgwbl[_x_];
  my_subbox[4] = subbox.Lgwbl[_y_];
  my_subbox[5] = subbox.Lgwbl[_z_];

#ifdef DEBUG
  char fname[SBLENGTH];
  sprintf(fname, "Task%d.dbg", ThisTask);
  DBGFD = fopen(fname, "a");
#endif

  /* let's synchronize the tasks here */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* these are the already stored particles */
#ifndef CLASSIC_FRAGMENTATION
  frag_offset = subbox.Nneeded;
#endif

  /* stores the data relative to the intersection of the fft box and subbox */
  if (keep_data(my_fft_box, my_subbox))
    return 1;

  comm_buffer = (product_data *)calloc(BUFLEN, sizeof(product_data));
  if (comm_buffer == 0x0)
  {
    printf("ERROR on task %d: could not allocate comm_buffer in distribute\n", ThisTask);
    fflush(stdout);
    return 1;
  }

  /* hypercubic communication scheme */
  for (log_ntask = 0; log_ntask < 1000; log_ntask++)
    if (1 << log_ntask >= NTasks)
      break;

  /* loop on hypercube dimension */
  for (bit = 1; bit < 1 << log_ntask; bit++)
  {
    /* loop on tasks */
    for (sender = 0; sender < NTasks; sender++)
    {
      /* receiver task is computed with a bitwise xor */
      receiver = sender ^ bit;

      /* condition on sender and receiver */
      if (receiver < NTasks && sender < receiver)
      {

        /* the communication will be done
     first sender -> receiver and then receiver -> sender */

        if (ThisTask == sender)
          send_data(my_fft_box, receiver);
        else if (ThisTask == receiver)
        {
          if (recv_data(my_subbox, sender))
            return 1;
        }

        if (ThisTask == receiver)
          send_data(my_fft_box, sender);
        else if (ThisTask == sender)
        {
          if (recv_data(my_subbox, receiver))
            return 1;
        }
      }
    }
  }

  /* updates the number of stored particles */
#ifdef CLASSIC_FRAGMENTATION
  subbox.Nstored = subbox.Npart;
#else

  if (frag_offset > subbox.Nalloc)
    subbox.Nstored = subbox.Nalloc;
  else
    subbox.Nstored = frag_offset;

  subbox.Nneeded = frag_offset;

#endif

  free(comm_buffer);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
  fclose(DBGFD);
#endif

  return 0;
}

/* -------------------------------------------------------------------
 * Helper types and functions for distribute_alltoall
 * ------------------------------------------------------------------- */

#ifndef CLASSIC_FRAGMENTATION

typedef struct
{
  int nsend_int;    /* intersections: my_fft_box ∩ partner_subbox */
  int send_int[48]; /* up to 8 intersection boxes × 6 ints each */
  int nrecv_int;    /* intersections: partner_fft_box ∩ my_subbox */
  int recv_int[48];
  int bm_sendcount; /* bitmap words I send (as receiver) to this peer */
  int bm_recvcount; /* bitmap words I receive (as sender) from this peer */
  int sendcount;    /* product particles I send */
  int recvcount;    /* product particles I receive */
} dist_peer_plan;

static inline void dist_append_frag(unsigned int spos, const product_data *src)
{
  if (frag_offset < (unsigned long long)subbox.Nalloc)
  {
    frag_pos[frag_offset] = spos;
    memcpy(&frag[frag_offset], src, sizeof(product_data));
  }
  ++frag_offset;
}

#endif /* !CLASSIC_FRAGMENTATION */

static int dist_intersection_copy(const int lhs[6], const int rhs[6], int *ibox)
{
  int a[6], b[6];
  memcpy(a, lhs, 6 * sizeof(int));
  memcpy(b, rhs, 6 * sizeof(int));
  return intersection(a, b, ibox);
}

static inline unsigned int dist_box_size(const int *box)
{
  return (unsigned int)box[3] * (unsigned int)box[4] * (unsigned int)box[5];
}

/* ---- Helper: exchange payload with one peer ----
   Performs a two-phase blocking send/recv with the given partner using
   the cached plan, bitmaps, and chunk buffers.  Called from both the
   default and staggered hypercube paths. */
#ifndef CLASSIC_FRAGMENTATION
static void peer_exchange(int partner, int send_first,
                          dist_peer_plan *plan,
                          unsigned int *bm_recvbuf, int *bm_rdispls,
                          dist_data *send_chunk, dist_data *recv_chunk,
                          MPI_Datatype dist_data_type)
{
  for (int phase = 0; phase < 2; phase++)
  {
    int do_send = (phase == 0) ? send_first : !send_first;

    if (do_send)
    {
      /* ---- SEND to partner ---- */
      if (plan[partner].sendcount == 0)
        continue;

      unsigned int *bm_ptr = bm_recvbuf + bm_rdispls[partner];
      int bufcount = 0;

      for (int box = 0; box < plan[partner].nsend_int; box++)
      {
        int *b = plan[partner].send_int + 6 * box;
        unsigned int size = dist_box_size(b);
        unsigned int mapl = distmap_length(size);

        for (unsigned int i = 0; i < size; i++)
        {
          if (!get_distmap_bit(bm_ptr, i))
            continue;

          unsigned int pfft = fft_space_index(i, b);
          unsigned int ip, jp, kp;
          int g[3];

          INDEX_TO_COORD(i, ip, jp, kp, b + 3);
          g[_x_] = (ip + b[_x_]) % MyGrids[0].GSglobal[_x_];
          g[_y_] = (jp + b[_y_]) % MyGrids[0].GSglobal[_y_];
          g[_z_] = (kp + b[_z_]) % MyGrids[0].GSglobal[_z_];

          send_chunk[bufcount].pos =
              COORD_TO_INDEX(g[_x_], g[_y_], g[_z_], MyGrids[0].GSglobal);
          memcpy(&send_chunk[bufcount].prod, &products[pfft], sizeof(product_data));
          bufcount++;

          if (bufcount == BUFLEN)
          {
            MPI_Send(send_chunk, bufcount, dist_data_type, partner, 0, MPI_COMM_WORLD);
            bufcount = 0;
          }
        }
        bm_ptr += mapl;
      }

      if (bufcount > 0)
        MPI_Send(send_chunk, bufcount, dist_data_type, partner, 0, MPI_COMM_WORLD);
    }
    else
    {
      /* ---- RECEIVE from partner ---- */
      int received = 0;
      int nrecv = plan[partner].recvcount;

      while (received < nrecv)
      {
        int chunk = nrecv - received;
        if (chunk > BUFLEN)
          chunk = BUFLEN;

        MPI_Recv(recv_chunk, chunk, dist_data_type, partner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int ri = 0; ri < chunk; ri++)
        {
          int g[3], l[3];
          INDEX_TO_COORD(recv_chunk[ri].pos,
                         g[_x_], g[_y_], g[_z_], MyGrids[0].GSglobal);

          l[_x_] = (g[_x_] - subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_]) % MyGrids[0].GSglobal[_x_];
          l[_y_] = (g[_y_] - subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_]) % MyGrids[0].GSglobal[_y_];
          l[_z_] = (g[_z_] - subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_]) % MyGrids[0].GSglobal[_z_];

          unsigned int spos = COORD_TO_INDEX(l[_x_], l[_y_], l[_z_], subbox.Lgwbl);
          dist_append_frag(spos, &recv_chunk[ri].prod);
        }
        received += chunk;
      }
    }
  }
}
#endif /* !CLASSIC_FRAGMENTATION */

int distribute_alltoall(void)
{
#ifdef CLASSIC_FRAGMENTATION
  if (!ThisTask)
    printf("ERROR: distribute_alltoall no longer supports CLASSIC_FRAGMENTATION\n");
  fflush(stdout);
  return 1;
#else
  /* Distributes products from fft-space to sub-volumes.
     All intersections and bitmap sizes are cached in a per-peer plan.
     Bitmap exchange uses non-blocking point-to-point (MPI_Isend/Irecv)
     to avoid the MPI_Alltoallv collective that scaled poorly.
     Data exchange uses the proven hypercube schedule with BUFLEN chunks. */

  int my_fft_box[6], my_subbox[6];
  size_t ntasks = 0;
  unsigned long long before_recv = 0ULL;

  int *all_fftboxes = NULL;
  int *all_subboxes = NULL;
  dist_peer_plan *plan = NULL;

  int *bm_sdispls = NULL;
  int *bm_rdispls = NULL;
  unsigned int *bm_sendbuf = NULL;
  unsigned int *bm_recvbuf = NULL;
  MPI_Request *reqs = NULL;

  int *sendcounts = NULL;
  int *recvcounts = NULL;

  dist_data *send_chunk = NULL;
  dist_data *recv_chunk = NULL;
  MPI_Datatype dist_data_type = MPI_DATATYPE_NULL;
  int dtype_committed = 0;

#ifdef NODE_AWARE_DISTRIBUTE
  int *node_id = NULL; /* node identifier per rank (identical on all tasks) */
#endif

  double t_phase, t_keep, t_plan, t_bm_build, t_bm_xchg, t_filter, t_alltoall, t_hypercube;
  double t_min, t_max;

  if (get_task_count_size(&ntasks))
    return 1;

  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  my_subbox[0] = subbox.stabl[_x_];
  my_subbox[1] = subbox.stabl[_y_];
  my_subbox[2] = subbox.stabl[_z_];
  my_subbox[3] = subbox.Lgwbl[_x_];
  my_subbox[4] = subbox.Lgwbl[_y_];
  my_subbox[5] = subbox.Lgwbl[_z_];

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  frag_offset = subbox.Nneeded;

  t_phase = MPI_Wtime();
  if (keep_data(my_fft_box, my_subbox))
    goto fail;

  before_recv = frag_offset;
  t_keep = MPI_Wtime() - t_phase;

  /* ---- Gather all boxes ---- */
  t_phase = MPI_Wtime();
  all_fftboxes = (int *)malloc(6 * ntasks * sizeof(int));
  all_subboxes = (int *)malloc(6 * ntasks * sizeof(int));
  plan = (dist_peer_plan *)calloc(ntasks, sizeof(dist_peer_plan));

  if (!all_fftboxes || !all_subboxes || !plan)
  {
    printf("ERROR on task %d: alloc failed in distribute_alltoall (boxes/plan)\n", ThisTask);
    fflush(stdout);
    goto fail;
  }

  MPI_Allgather(my_fft_box, 6, MPI_INT, all_fftboxes, 6, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(my_subbox, 6, MPI_INT, all_subboxes, 6, MPI_INT, MPI_COMM_WORLD);

  /* ---- Build per-peer intersection plan (computed once, reused) ---- */
  for (int peer = 0; peer < NTasks; peer++)
  {
    if (peer == ThisTask)
      continue;

    plan[peer].nsend_int =
        dist_intersection_copy(my_fft_box, all_subboxes + 6 * peer, plan[peer].send_int);
    plan[peer].nrecv_int =
        dist_intersection_copy(all_fftboxes + 6 * peer, my_subbox, plan[peer].recv_int);

    /* Bitmap sizes computed locally — no communication needed.
       bm_sendcount: I send need-maps for recv_int (peer_fft ∩ my_subbox).
       bm_recvcount: I receive need-maps for send_int (my_fft ∩ peer_subbox).
       Both sides agree on intersection layout, so sizes match. */
    plan[peer].bm_sendcount = 0;
    for (int box = 0; box < plan[peer].nrecv_int; box++)
      plan[peer].bm_sendcount += (int)distmap_length(dist_box_size(plan[peer].recv_int + 6 * box));

    plan[peer].bm_recvcount = 0;
    for (int box = 0; box < plan[peer].nsend_int; box++)
      plan[peer].bm_recvcount += (int)distmap_length(dist_box_size(plan[peer].send_int + 6 * box));
  }

  free(all_fftboxes);
  all_fftboxes = NULL;
  free(all_subboxes);
  all_subboxes = NULL;
  t_plan = MPI_Wtime() - t_phase;

  /* ---- Allocate bitmap flat buffers ---- */
  bm_sdispls = (int *)calloc(ntasks, sizeof(int));
  bm_rdispls = (int *)calloc(ntasks, sizeof(int));

  if (!bm_sdispls || !bm_rdispls)
  {
    printf("ERROR on task %d: alloc failed in distribute_alltoall (bm_displs)\n", ThisTask);
    fflush(stdout);
    goto fail;
  }

  {
    int bm_total_send = 0, bm_total_recv = 0;
    for (int peer = 0; peer < NTasks; peer++)
    {
      bm_sdispls[peer] = bm_total_send;
      bm_total_send += plan[peer].bm_sendcount;
      bm_rdispls[peer] = bm_total_recv;
      bm_total_recv += plan[peer].bm_recvcount;
    }

    if (bm_total_send > 0)
    {
      bm_sendbuf = (unsigned int *)calloc((size_t)bm_total_send, sizeof(unsigned int));
      if (!bm_sendbuf)
      {
        printf("ERROR on task %d: alloc failed in distribute_alltoall (bm_sendbuf)\n", ThisTask);
        fflush(stdout);
        goto fail;
      }
    }

    if (bm_total_recv > 0)
    {
      bm_recvbuf = (unsigned int *)calloc((size_t)bm_total_recv, sizeof(unsigned int));
      if (!bm_recvbuf)
      {
        printf("ERROR on task %d: alloc failed in distribute_alltoall (bm_recvbuf)\n", ThisTask);
        fflush(stdout);
        goto fail;
      }
    }
  }

  /* ---- Build receiver-side need-bitmaps ---- */
  t_phase = MPI_Wtime();
  for (int peer = 0; peer < NTasks; peer++)
  {
    if (peer == ThisTask || plan[peer].bm_sendcount == 0)
      continue;

    int offset = bm_sdispls[peer];
    for (int box = 0; box < plan[peer].nrecv_int; box++)
    {
      int *b = plan[peer].recv_int + 6 * box;
      unsigned int mapl = distmap_length(dist_box_size(b));
      build_distmap(bm_sendbuf + offset, b);
      offset += (int)mapl;
    }
  }

  t_bm_build = MPI_Wtime() - t_phase;

  /* ---- Exchange bitmaps via non-blocking point-to-point ---- */
  t_phase = MPI_Wtime();
  reqs = (MPI_Request *)malloc(2 * ntasks * sizeof(MPI_Request));
  if (!reqs)
  {
    printf("ERROR on task %d: alloc failed in distribute_alltoall (reqs)\n", ThisTask);
    fflush(stdout);
    goto fail;
  }

  {
    int nreq = 0;

    /* Post all receives first */
    for (int peer = 0; peer < NTasks; peer++)
    {
      if (peer == ThisTask || plan[peer].bm_recvcount == 0)
        continue;
      MPI_Irecv(bm_recvbuf + bm_rdispls[peer], plan[peer].bm_recvcount,
                MPI_UNSIGNED, peer, 1, MPI_COMM_WORLD, &reqs[nreq++]);
    }

    /* Post all sends */
    for (int peer = 0; peer < NTasks; peer++)
    {
      if (peer == ThisTask || plan[peer].bm_sendcount == 0)
        continue;
      MPI_Isend(bm_sendbuf + bm_sdispls[peer], plan[peer].bm_sendcount,
                MPI_UNSIGNED, peer, 1, MPI_COMM_WORLD, &reqs[nreq++]);
    }

    MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
  }

  free(reqs);
  reqs = NULL;
  free(bm_sendbuf);
  bm_sendbuf = NULL;
  free(bm_sdispls);
  bm_sdispls = NULL;
  t_bm_xchg = MPI_Wtime() - t_phase;

  /* ---- Sender-side Fmax filtering + particle count ---- */
  t_phase = MPI_Wtime();
  sendcounts = (int *)calloc(ntasks, sizeof(int));
  recvcounts = (int *)calloc(ntasks, sizeof(int));

  if (!sendcounts || !recvcounts)
  {
    printf("ERROR on task %d: alloc failed in distribute_alltoall (counts)\n", ThisTask);
    fflush(stdout);
    goto fail;
  }

  for (int peer = 0; peer < NTasks; peer++)
  {
    if (peer == ThisTask || plan[peer].bm_recvcount == 0)
      continue;

    unsigned int *bm_ptr = bm_recvbuf + bm_rdispls[peer];
    int count = 0;

    for (int box = 0; box < plan[peer].nsend_int; box++)
    {
      int *b = plan[peer].send_int + 6 * box;
      unsigned int size = dist_box_size(b);
      unsigned int mapl = distmap_length(size);

      update_distmap(bm_ptr, b);
      for (unsigned int i = 0; i < size; i++)
        if (get_distmap_bit(bm_ptr, i))
          count++;
      bm_ptr += mapl;
    }

    plan[peer].sendcount = count;
    sendcounts[peer] = count;
  }

  t_filter = MPI_Wtime() - t_phase;

  t_phase = MPI_Wtime();
  MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);
  t_alltoall = MPI_Wtime() - t_phase;

  for (int peer = 0; peer < NTasks; peer++)
    plan[peer].recvcount = recvcounts[peer];

  /* ---- Create MPI datatype for dist_data ---- */
  MPI_Type_contiguous((int)sizeof(dist_data), MPI_BYTE, &dist_data_type);
  MPI_Type_commit(&dist_data_type);
  dtype_committed = 1;

#ifdef NODE_AWARE_DISTRIBUTE
  /* ---- Build a global node-id array so that every rank can test
         whether any two ranks (a, b) share a node, not just ThisTask
         vs. some peer.  node_id[r] is the smallest world-rank on the
         node that contains rank r. ---- */
  {
    MPI_Comm node_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, ThisTask,
                        MPI_INFO_NULL, &node_comm);

    /* Find the minimum world rank on this node = our node leader */
    int my_node_leader;
    MPI_Allreduce(&ThisTask, &my_node_leader, 1, MPI_INT, MPI_MIN, node_comm);

    node_id = (int *)malloc(ntasks * sizeof(int));
    if (!node_id)
    {
      printf("ERROR on task %d: alloc failed in distribute_alltoall (node_id)\n", ThisTask);
      fflush(stdout);
      MPI_Comm_free(&node_comm);
      goto fail;
    }

    /* Gather so every rank knows every other rank's node leader */
    MPI_Allgather(&my_node_leader, 1, MPI_INT, node_id, 1, MPI_INT, MPI_COMM_WORLD);

    int node_size;
    MPI_Comm_size(node_comm, &node_size);
    int num_local = node_size - 1; /* excluding self */
    int num_remote = NTasks - node_size;
    if (!ThisTask)
      printf("  distribute_alltoall node-aware: %d local peers, %d remote peers per task\n",
             num_local, num_remote);

    MPI_Comm_free(&node_comm);
  }
#endif /* NODE_AWARE_DISTRIBUTE */

#ifndef MULTILAYER_DISTRIBUTE
  /* ---- Allocate chunk buffers ---- */
  send_chunk = (dist_data *)malloc(BUFLEN * sizeof(dist_data));
  recv_chunk = (dist_data *)malloc(BUFLEN * sizeof(dist_data));

  if (!send_chunk || !recv_chunk)
  {
    printf("ERROR on task %d: alloc failed in distribute_alltoall (chunks)\n", ThisTask);
    fflush(stdout);
    goto fail;
  }

  /* ---- Hypercube data exchange using cached plan ---- */
  t_phase = MPI_Wtime();
  {
    int log_ntask;
    for (log_ntask = 0; log_ntask < 1000; log_ntask++)
      if ((1 << log_ntask) >= NTasks)
        break;

#ifdef NODE_AWARE_DISTRIBUTE
    /* Two passes: pass 0 = intra-node only, pass 1 = inter-node only */
    for (int node_pass = 0; node_pass < 2; node_pass++)
    {
      if (!ThisTask && node_pass == 0)
        printf("  distribute_alltoall node-aware: intra-node pass\n");
      else if (!ThisTask && node_pass == 1)
        printf("  distribute_alltoall node-aware: inter-node pass\n");
#endif

      for (int bit = 1; bit < (1 << log_ntask); bit++)
      {
        int partner = ThisTask ^ bit;
        if (partner >= NTasks)
          continue;

#ifdef NODE_AWARE_DISTRIBUTE
        /* Skip partner if it doesn't belong to this node-aware pass */
        if (node_pass == 0 && node_id[partner] != node_id[ThisTask])
          continue;
        if (node_pass == 1 && node_id[partner] == node_id[ThisTask])
          continue;
#endif

        int send_first = (ThisTask < partner);
        peer_exchange(partner, send_first, plan, bm_recvbuf, bm_rdispls,
                      send_chunk, recv_chunk, dist_data_type);
      }

#ifdef NODE_AWARE_DISTRIBUTE
      /* Strict phase boundary: no inter-node traffic begins until all
         intra-node exchanges are complete everywhere. */
      if (node_pass == 0)
        MPI_Barrier(MPI_COMM_WORLD);
    } /* end node_pass loop */
#endif
  }

  t_hypercube = MPI_Wtime() - t_phase;

#else /* MULTILAYER_DISTRIBUTE */

#ifndef MAX_CONCURRENT_PAIRS
#define MAX_CONCURRENT_PAIRS 16
#endif

  /* ---- Staggered hypercube data exchange ---- */
  /* Same blocking MPI_Send/MPI_Recv as the default path, but within each
     hypercube step pairs are split into sub-rounds so that at most
     MAX_CONCURRENT_PAIRS pairs communicate simultaneously.  A barrier
     separates sub-rounds to enforce the global concurrency limit.
     Without NODE_AWARE_DISTRIBUTE all pairs are staggered.
     With NODE_AWARE_DISTRIBUTE only the inter-node pass is staggered;
     intra-node pairs exchange at full concurrency.
     No extra memory beyond send_chunk/recv_chunk. */

  send_chunk = (dist_data *)malloc(BUFLEN * sizeof(dist_data));
  recv_chunk = (dist_data *)malloc(BUFLEN * sizeof(dist_data));

  if (!send_chunk || !recv_chunk)
  {
    printf("ERROR on task %d: alloc failed in distribute_alltoall (chunks/staggered)\n", ThisTask);
    fflush(stdout);
    goto fail;
  }

  t_phase = MPI_Wtime();
  {
    int log_ntask;
    for (log_ntask = 0; log_ntask < 1000; log_ntask++)
      if ((1 << log_ntask) >= NTasks)
        break;

#ifdef NODE_AWARE_DISTRIBUTE
    /* Two passes: pass 0 = intra-node (full concurrency),
                   pass 1 = inter-node (staggered with MAX_CONCURRENT_PAIRS) */
    for (int node_pass = 0; node_pass < 2; node_pass++)
    {
      if (node_pass == 0)
      {
        /* ---- Intra-node pass: no throttling ---- */
        if (!ThisTask)
          printf("  distribute_alltoall node-aware: intra-node pass (full concurrency)\n");

        for (int bit = 1; bit < (1 << log_ntask); bit++)
        {
          int partner = ThisTask ^ bit;
          if (partner >= NTasks)
            continue;
          if (node_id[partner] != node_id[ThisTask])
            continue;

          int send_first = (ThisTask < partner);
          peer_exchange(partner, send_first, plan, bm_recvbuf, bm_rdispls,
                        send_chunk, recv_chunk, dist_data_type);
        }

        /* Strict phase boundary: no inter-node traffic begins until all
           intra-node exchanges are complete everywhere. */
        MPI_Barrier(MPI_COMM_WORLD);
      }
      else
      {
        /* ---- Inter-node pass: staggered with MAX_CONCURRENT_PAIRS ---- */
        if (!ThisTask)
          printf("  distribute_alltoall node-aware: inter-node pass"
                 " (staggered, MAX_CONCURRENT_PAIRS=%d)\n",
                 MAX_CONCURRENT_PAIRS);
#endif /* NODE_AWARE_DISTRIBUTE */

        for (int bit = 1; bit < (1 << log_ntask); bit++)
        {
          int partner = ThisTask ^ bit;

          /* Count active pairs for this bit.  When NODE_AWARE_DISTRIBUTE
             is active we are in the inter-node pass, so count only
             inter-node pairs using the global node_id[]. */
          int num_active_pairs = 0;
          for (int t = 0; t < NTasks; t++)
          {
            int p = t ^ bit;
            if (p < NTasks && t < p)
            {
#ifdef NODE_AWARE_DISTRIBUTE
              if (node_id[t] == node_id[p])
                continue; /* skip intra-node pairs */
#endif
              num_active_pairs++;
            }
          }

          int num_sub_rounds = (num_active_pairs + MAX_CONCURRENT_PAIRS - 1) / MAX_CONCURRENT_PAIRS;

          /* Compute which sub-round this task's pair belongs to */
          int my_sub_round = -1;
          if (partner < NTasks)
          {
#ifdef NODE_AWARE_DISTRIBUTE
            /* Only proceed for remote partners (inter-node) */
            if (node_id[ThisTask] != node_id[partner])
            {
#endif
              int my_pair_id = 0;
              int low = (ThisTask < partner) ? ThisTask : partner;
              for (int t = 0; t < low; t++)
              {
                int p = t ^ bit;
                if (p < NTasks && t < p)
                {
#ifdef NODE_AWARE_DISTRIBUTE
                  if (node_id[t] == node_id[p])
                    continue;
#endif
                  my_pair_id++;
                }
              }
              my_sub_round = my_pair_id / MAX_CONCURRENT_PAIRS;
#ifdef NODE_AWARE_DISTRIBUTE
            }
#endif
          }

          for (int sr = 0; sr < num_sub_rounds; sr++)
          {
            if (sr == my_sub_round && partner < NTasks)
            {
              int send_first = (ThisTask < partner);
              peer_exchange(partner, send_first, plan, bm_recvbuf, bm_rdispls,
                            send_chunk, recv_chunk, dist_data_type);
            }

            /* Synchronize before next sub-round so that at most
               MAX_CONCURRENT_PAIRS pairs are ever active at once */
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }

#ifdef NODE_AWARE_DISTRIBUTE
      } /* end else (inter-node pass) */
    } /* end node_pass loop */
#endif
  }

  t_hypercube = MPI_Wtime() - t_phase;

#ifdef NODE_AWARE_DISTRIBUTE
  if (!ThisTask)
    printf("  distribute_alltoall: intra-node full concurrency,"
           " inter-node staggered (MAX_CONCURRENT_PAIRS=%d, NTasks=%d)\n",
           MAX_CONCURRENT_PAIRS, NTasks);
#else
  if (!ThisTask)
    printf("  distribute_alltoall: staggered exchange"
           " (MAX_CONCURRENT_PAIRS=%d, NTasks=%d)\n",
           MAX_CONCURRENT_PAIRS, NTasks);
#endif

#endif /* MULTILAYER_DISTRIBUTE */

  /* ---- Phase timing report ---- */
  if (!ThisTask)
    printf("  distribute_alltoall phase timings (task 0):\n"
           "    keep_data:     %10.3f s\n"
           "    plan+gather:   %10.3f s\n"
           "    bitmap build:  %10.3f s\n"
           "    bitmap xchg:   %10.3f s\n"
           "    Fmax filter:   %10.3f s\n"
           "    alltoall cnt:  %10.3f s\n"
           "    hypercube:     %10.3f s\n",
           t_keep, t_plan, t_bm_build, t_bm_xchg, t_filter, t_alltoall, t_hypercube);

  /* Report min/max across tasks for each phase to detect outliers */
  {
    double my_times[7] = {t_keep, t_plan, t_bm_build, t_bm_xchg, t_filter, t_alltoall, t_hypercube};
    double min_times[7], max_times[7];
    MPI_Reduce(my_times, min_times, 7, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(my_times, max_times, 7, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!ThisTask)
    {
      const char *labels[7] = {"keep_data", "plan+gather", "bitmap build",
                               "bitmap xchg", "Fmax filter", "alltoall cnt", "hypercube"};
      printf("  distribute_alltoall min/max across %d tasks:\n", NTasks);
      for (int ph = 0; ph < 7; ph++)
        printf("    %-14s  min=%10.3f  max=%10.3f  (spread=%.1fx)\n",
               labels[ph], min_times[ph], max_times[ph],
               min_times[ph] > 0.001 ? max_times[ph] / min_times[ph] : 0.0);
    }
  }

  /* ---- Overflow check + finalize ---- */
  {
    unsigned long long accepted = frag_offset - before_recv;
    unsigned long long cap_left =
        (before_recv >= (unsigned long long)subbox.Nalloc)
            ? 0ULL
            : (unsigned long long)subbox.Nalloc - before_recv;

    if (accepted > cap_left)
    {
      printf("WARNING on task %d: distribute_alltoall accepted %llu particles but only %llu slots available\n",
             ThisTask, accepted, cap_left);
      fflush(stdout);
    }
  }

  subbox.Nstored = (frag_offset > (unsigned long long)subbox.Nalloc)
                       ? subbox.Nalloc
                       : frag_offset;
  subbox.Nneeded = frag_offset;

  /* ---- Cleanup (success path) ---- */
  if (dtype_committed)
    MPI_Type_free(&dist_data_type);
  free(recv_chunk);
  free(send_chunk);
  free(recvcounts);
  free(sendcounts);
  free(bm_recvbuf);
  free(bm_rdispls);
  free(plan);
#ifdef NODE_AWARE_DISTRIBUTE
  free(node_id);
#endif

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;

fail:
  if (dtype_committed)
    MPI_Type_free(&dist_data_type);
  free(recv_chunk);
  free(send_chunk);
  free(recvcounts);
  free(sendcounts);
  free(bm_recvbuf);
  free(bm_sendbuf);
  free(bm_rdispls);
  free(bm_sdispls);
  free(reqs);
  free(plan);
  free(all_subboxes);
  free(all_fftboxes);
#ifdef NODE_AWARE_DISTRIBUTE
  free(node_id);
#endif
  return 1;
#endif /* CLASSIC_FRAGMENTATION */
}

int intersection(int *fbox, int *sbox, int *ibox)
{

  unsigned int istart[3], istop[3], istart2[3], istop2[3], dim, Nint,
      stop1, stop2, this, off;
  unsigned int ISTARTx, ISTOPx, ISTARTy, ISTOPy, ISTARTz, ISTOPz, ax, ay, az;

  /* intersection of the two boxes is considered dimension by dimension */
  Nint = 1;
  for (dim = 0; dim < 3; dim++)
  {
    /* in case, fix negative starting point */
    if (sbox[dim] < 0)
      sbox[dim] += MyGrids[0].GSglobal[dim];

    /* first stopping point for subbox is at most the global box edge */
    stop1 = fbox[dim] + fbox[dim + 3];
    if (sbox[dim] + sbox[dim + 3] > MyGrids[0].GSglobal[dim])
      stop2 = MyGrids[0].GSglobal[dim];
    else
      stop2 = sbox[dim] + sbox[dim + 3];

    /* intersection up to the global box edge */
    istart[dim] = (fbox[dim] > sbox[dim] ? fbox[dim] : sbox[dim]);
    istop[dim] = (stop1 < stop2 ? stop1 : stop2);

    /* if the subbox goes beyond the global box edge,
 apply PBCs to the other segment and check the intersection */
    if ((stop2 = sbox[dim] + sbox[dim + 3]) > MyGrids[0].GSglobal[dim])
    {
      stop2 = stop2 % MyGrids[0].GSglobal[dim];
      istart2[dim] = (fbox[dim] > 0 ? fbox[dim] : 0);
      istop2[dim] = (stop1 < stop2 ? stop1 : stop2);
    }
    else
    {
      istart2[dim] = 1;
      istop2[dim] = 0;
    }

    /* this dimension contributes 0, 1 or 2 */
    Nint *= (istart[dim] < istop[dim]) + (istart2[dim] < istop2[dim]);
  }

  /* store all intersections, looping on the two options for each dimension */
  if (Nint)
  {
    this = 0;
    for (ax = 0; ax < 2; ax++)
    {
      if (ax)
      {
        ISTARTx = istart[0];
        ISTOPx = istop[0];
      }
      else
      {
        ISTARTx = istart2[0];
        ISTOPx = istop2[0];
      }
      for (ay = 0; ay < 2; ay++)
      {
        if (ay)
        {
          ISTARTy = istart[1];
          ISTOPy = istop[1];
        }
        else
        {
          ISTARTy = istart2[1];
          ISTOPy = istop2[1];
        }
        for (az = 0; az < 2; az++)
        {
          if (az)
          {
            ISTARTz = istart[2];
            ISTOPz = istop[2];
          }
          else
          {
            ISTARTz = istart2[2];
            ISTOPz = istop2[2];
          }

          if ((ISTARTx < ISTOPx) & (ISTARTy < ISTOPy) & (ISTARTz < ISTOPz))
          {
            off = this * 6;
            ibox[off] = ISTARTx;
            ibox[1 + off] = ISTARTy;
            ibox[2 + off] = ISTARTz;
            ibox[3 + off] = ISTOPx - ISTARTx;
            ibox[4 + off] = ISTOPy - ISTARTy;
            ibox[5 + off] = ISTOPz - ISTARTz;
            this++;
          }
        }
      }
    }
  }

#ifdef DEBUG
  fprintf(DBGFD, "INTERSECTION\n");
  fprintf(DBGFD, "Task %d, fft box:        %d %d %d   %d %d %d\n",
          ThisTask, fbox[0], fbox[1], fbox[2], fbox[3], fbox[4], fbox[5]);
  fprintf(DBGFD, "         subvolume:      %d %d %d   %d %d %d\n",
          sbox[0], sbox[1], sbox[2], sbox[3], sbox[4], sbox[5]);
  for (this = 0; this < Nint; this++)
  {
    off = this * 6;
    fprintf(DBGFD, "         intersection %d: %d %d %d   %d %d %d\n",
            this, ibox[off], ibox[1 + off], ibox[2 + off], ibox[3 + off], ibox[4 + off], ibox[5 + off]);
  }
  if (!Nint)
    fprintf(DBGFD, "        no intersections\n");
#endif

  return Nint;
}

int send_data(int *mybox, int target)
{
  /* This routine sends the content of a box to a target task.
     Communication is divided in these stages:
     1) the sender communicates the start and length of its box,
     2) the sender receives the start and the length of the needed box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the sender receives a map of needed particles
     5) the sender loops on particles and sends them in chunks of size BUFLEN
  */

#ifndef CLASSIC_FRAGMENTATION
  unsigned int *map, mapl;
#endif
  unsigned int size, bufcount;
  int targetbox[6], interbox[48], i, Nint, box, off;
  MPI_Status status;

  MPI_Send(mybox, 6, MPI_INT, target, 0, MPI_COMM_WORLD);
  MPI_Recv(targetbox, 6, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

  /* up to 8 intersections of the two boxes */
  Nint = intersection(mybox, targetbox, interbox);
  for (box = 0; box < Nint; box++)
  {
    off = box * 6;
    size = interbox[3 + off] * interbox[4 + off] * interbox[5 + off];

#ifdef DEBUG
    fprintf(DBGFD, "Task %d will send this box: %d %d %d -- %d %d %d\n",
            ThisTask, interbox[0 + off], interbox[1 + off], interbox[2 + off],
            interbox[3 + off], interbox[4 + off], interbox[5 + off]);
#endif

#ifndef CLASSIC_FRAGMENTATION
    /* receive the map of needed particles */
    mapl = distmap_length(size);
    map = (unsigned int *)calloc(mapl, sizeof(unsigned int));

    /* the receiver has built its map and is sending it */
    MPI_Recv(map, mapl, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

    /* updates the map and sends it back to the target */
    update_distmap(map, interbox + off);
    MPI_Send(map, mapl, MPI_INT, target, 0, MPI_COMM_WORLD);
#endif

    /* load particles on the buffer and send them to the target */
    bufcount = 0;
    for (i = 0; i < size; i++)
    {
#ifndef CLASSIC_FRAGMENTATION
      if (get_distmap_bit(map, i))
#endif
      {
        memcpy(&comm_buffer[bufcount], &products[fft_space_index(i, interbox + off)], sizeof(product_data));

        bufcount++;
      }
      if (bufcount == BUFLEN)
      {
#ifdef DEBUG
        fprintf(DBGFD, "...sending %d products to Task %d... %d\n", bufcount, target, (int)sizeof(product_data));
        if (bufcount < 10)
        {
          for (int u = 0; u < bufcount; u++)
            fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
          fprintf(DBGFD, "\n");
        }
        else
        {
          for (int u = 0; u < 5; u++)
            fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
          fprintf(DBGFD, " ... ");
          for (int u = bufcount - 5; u < bufcount; u++)
            fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
          fprintf(DBGFD, "\n");
        }
#endif
        MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
        MPI_Send(comm_buffer, bufcount * sizeof(product_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
        bufcount = 0;
      }
    }

    if (bufcount)
    {
#ifdef DEBUG
      fprintf(DBGFD, "...sending %d products to Task %d...\n", bufcount, target);
      if (bufcount < 10)
      {
        for (int u = 0; u < bufcount; u++)
          fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
        fprintf(DBGFD, "\n");
      }
      else
      {
        for (int u = 0; u < 5; u++)
          fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
        fprintf(DBGFD, " ... ");
        for (int u = bufcount - 5; u < bufcount; u++)
          fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
        fprintf(DBGFD, "\n");
      }
#endif
      MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
      MPI_Send(comm_buffer, bufcount * sizeof(product_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
    }

#ifndef CLASSIC_FRAGMENTATION
    free(map);
#endif
  }

  /* done */
  return 0;
}

int recv_data(int *mybox, int sender)
{
  /* This routine receives the content of a box from a sender task.
     Communication is divided in these stages:
     1) the target receives the start and length of sender box,
     2) the target sends the start and the length of its box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the target sends a map of needed particles
     5) the target loops on particles and receives them in chunks of size BUFLEN
  */

#ifndef CLASSIC_FRAGMENTATION
  unsigned int *map, mapl;
  int nstore;
#endif
  unsigned int size, nsent, boff, expected, received;
  int senderbox[6], interbox[48], i, Nint, box;
  MPI_Status status;

  MPI_Recv(senderbox, 6, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
  MPI_Send(mybox, 6, MPI_INT, sender, 0, MPI_COMM_WORLD);

  /* up to 8 intersections of the two boxes */
  Nint = intersection(senderbox, mybox, interbox);
  for (box = 0; box < Nint; box++)
  {
    boff = box * 6;
    size = interbox[3 + boff] * interbox[4 + boff] * interbox[5 + boff];

#ifdef DEBUG
    fprintf(DBGFD, "Task %d will receive this box: %d %d %d -- %d %d %d\n",
            ThisTask, interbox[0 + boff], interbox[1 + boff], interbox[2 + boff],
            interbox[3 + boff], interbox[4 + boff], interbox[5 + boff]);
#endif

#ifndef CLASSIC_FRAGMENTATION
    /* send the map of needed particles */
    mapl = distmap_length(size);
    map = (unsigned int *)calloc(mapl, sizeof(unsigned int));

    /* constructs the map for the intersection and send it to the sender */
    build_distmap(map, interbox + boff);
    MPI_Send(map, mapl, MPI_INT, sender, 0, MPI_COMM_WORLD);

    /* map is nulled before getting it back */
    memset(map, 0, mapl * sizeof(unsigned int));

    /* the map is updated by the sender and sent back */
    MPI_Recv(map, mapl, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);

    /* record positions of particles that will be stored
 and count how many particles will be sent */
    int off = frag_offset;
    int count = frag_offset;

    for (i = 0; i < size; i++)
      if (get_distmap_bit(map, i))
      {
        ++count;
        if (off < subbox.Nalloc)
          frag_pos[off++] = subbox_space_index(i, interbox + boff);
      }

    expected = count - frag_offset;
#else
    expected = size;
#endif

    /* receive particles */
    received = 0;

    if (expected)
    {

      do
      {

        MPI_Recv(&nsent, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(comm_buffer, nsent * sizeof(product_data), MPI_BYTE, sender, 0, MPI_COMM_WORLD, &status);
#ifdef DEBUG
        fprintf(DBGFD, "...receiving %d products from Task %d... %d\n", nsent, sender, (int)sizeof(product_data));
        if (nsent < 10)
        {
          for (int u = 0; u < nsent; u++)
            fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
          fprintf(DBGFD, "\n");
        }
        else
        {
          for (int u = 0; u < 5; u++)
            fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
          fprintf(DBGFD, " ... ");
          for (int u = nsent - 5; u < nsent; u++)
            fprintf(DBGFD, "  %f  ", comm_buffer[u].Fmax);
          fprintf(DBGFD, "\n");
        }
#endif

#ifdef CLASSIC_FRAGMENTATION
        for (i = 0; i < nsent; i++)
          memcpy(&frag[subbox_space_index(i + received, interbox + boff)], &comm_buffer[i], sizeof(product_data));
#else
        if (frag_offset + nsent < subbox.Nalloc)
          nstore = nsent;
        else
          nstore = (long long)subbox.Nalloc - (long long)frag_offset;

        for (i = 0; i < nstore; i++)
          memcpy(&frag[frag_offset + i], &comm_buffer[i], sizeof(product_data));
        frag_offset += nsent;
#endif
        received += nsent;

      } while (received < expected);
    }

#ifndef CLASSIC_FRAGMENTATION
    free(map);
#endif
  }

  /* done */
  return 0;
}

int keep_data(int *fft_box, int *sub_box)
{
  /* this routine transfers products from the fft space to the subbox space */

#ifndef CLASSIC_FRAGMENTATION
  unsigned int *map, mapl;
#endif
  int interbox[48], i, Nint, box, off;
  unsigned int size;

  Nint = intersection(fft_box, sub_box, interbox);
  for (box = 0; box < Nint; box++)
  {
    off = box * 6;
    size = interbox[3 + off] * interbox[4 + off] * interbox[5 + off];

#ifdef DEBUG
    fprintf(DBGFD, "Task %d will keep this box: %d %d %d -- %d %d %d\n",
            ThisTask, interbox[0 + off], interbox[1 + off], interbox[2 + off], interbox[3 + off], interbox[4 + off], interbox[5 + off]);
#endif

#ifdef CLASSIC_FRAGMENTATION
    /* copy data */
    for (i = 0; i < size; i++)
      memcpy(&frag[subbox_space_index(i, interbox + off)],
             &products[fft_space_index(i, interbox + off)], sizeof(product_data));
#else

    /* map of needed particles */
    mapl = distmap_length(size);
    map = (unsigned int *)calloc(mapl, sizeof(unsigned int));
    build_distmap(map, interbox + off);
    update_distmap(map, interbox + off);

    /* copy data */
    for (i = 0; i < size; i++)
    {
      if (get_distmap_bit(map, i))
      {
        if (frag_offset < subbox.Nalloc)
        {
          frag_pos[frag_offset] = subbox_space_index(i, interbox + off);
          memcpy(&frag[frag_offset],
                 &products[fft_space_index(i, interbox + off)], sizeof(product_data));
        }
        ++frag_offset;
      }
    }
    free(map);
#endif
  }

  return 0;
}

unsigned int fft_space_index(unsigned int pos, int *box)
{
  /* here we move:
     (1) from index i to the relative position of the point in intersection,
     (2) from that to global position without PBCs,
     (3) then we impose PBCs
     (4) then we compute the position in the local FFT box
     (5) and finally we compute the local particle index
     (these are the logical steps, formulas are more compact)
  */

  unsigned int ip, jp, kp;

  INDEX_TO_COORD(pos, ip, jp, kp, (box + 3));

  return COORD_TO_INDEX((ip + box[0]) % MyGrids[0].GSglobal[_x_] - MyGrids[0].GSstart[_x_],
                        (jp + box[1]) % MyGrids[0].GSglobal[_y_] - MyGrids[0].GSstart[_y_],
                        (kp + box[2]) % MyGrids[0].GSglobal[_z_] - MyGrids[0].GSstart[_z_],
                        MyGrids[0].GSlocal);
}

unsigned int subbox_space_index(unsigned int pos, int *box)
{
  /* here we go
     (1) from index pos to relative position in intersection,
     (2) from that to global position, imposing PBCs,
     (3) then we compute the position in the local subbox, imposing PBCs again,
     (4) and finally the index */

  unsigned int ip, jp, kp;

  INDEX_TO_COORD(pos, ip, jp, kp, (box + 3)); /* coords within the intersection */

  return COORD_TO_INDEX(((ip + box[_x_]) % MyGrids[0].GSglobal[_x_] - subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_]) % MyGrids[0].GSglobal[_x_],
                        ((jp + box[_y_]) % MyGrids[0].GSglobal[_y_] - subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_]) % MyGrids[0].GSglobal[_y_],
                        ((kp + box[_z_]) % MyGrids[0].GSglobal[_z_] - subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_]) % MyGrids[0].GSglobal[_z_],
                        subbox.Lgwbl);
}

#ifndef CLASSIC_FRAGMENTATION
unsigned int distmap_length(unsigned int size)
{
  return size / UINTLEN + (size % UINTLEN != 0);
}

int get_distmap_bit(unsigned int *map, unsigned int pos)
{
  /* this operates on the map allocated for the distribution */
  /* gets the bit corresponding to position pos */
  unsigned int rem = pos % UINTLEN;
  return (map[pos / UINTLEN] & (1 << rem)) >> rem;
}

void set_distmap_bit(unsigned int *map, unsigned int pos, int value)
{
  /* this operates on the map allocated for the distribution */
  /* sets to 1 the bit corresponding to position pos */
  if (value)
    map[pos / UINTLEN] |= (1 << pos % UINTLEN);
  else
    map[pos / UINTLEN] &= ~(1 << pos % UINTLEN);
}

void build_distmap(unsigned int *map, int *box)
{
  /* the receiver builds the map used for distribution */
  unsigned int size = box[3] * box[4] * box[5];

  if (map_to_be_used)
    for (unsigned int i = 0; i < size; i++)
      set_distmap_bit(map, i, get_map_bit(subbox_space_index(i, box)));
  else
    for (unsigned int i = 0; i < size; i++)
      set_distmap_bit(map, i, get_mapup_bit(subbox_space_index(i, box)));
}

void update_distmap(unsigned int *map, int *box)
{

  unsigned int size = box[3] * box[4] * box[5];
  unsigned int bit, i;

  for (i = 0; i < size; i++)
  {

    bit = get_distmap_bit(map, i);
    set_distmap_bit(map, i, (bit & (products[fft_space_index(i, box)].Fmax >= outputs.Flast)));
  }
}

#endif

#ifdef SNAPSHOT
/* distribution of zacc back to the FFT space */

typedef struct
{
  unsigned int pos;
  PRODFLOAT zacc;
  int group_ID;
} back_data;
back_data *back_buffer;

int keep_data_back(int *);
int send_data_back(int);
int recv_data_back(int *, int);

int distribute_back_alltoall(void)
{
  /* Distributes accretion times from sub-volumes to FFT-space
     using a single MPI_Alltoallv instead of explicit pairwise sends.
   */

  int my_fft_box[6];
  int This = ThisTask;
  size_t ntasks;

  if (get_task_count_size(&ntasks))
    return 1;

  /* This defines the box that the task possesses in the FFT space */
  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* Local intersection work, as before */
  if (keep_data_back(my_fft_box))
    return 1;

  /* back_buffer is no longer used for communication, but may still be
     needed elsewhere, so keep allocation for compatibility. */
  back_buffer = (back_data *)calloc(BUFLEN, sizeof(back_data));
  if (back_buffer == NULL)
  {
    printf("ERROR on task %d: could not allocate back_buffer in distribute_back_alltoall\n", This);
    fflush(stdout);
    return 1;
  }

  /* ------------------------------------------------------------------
   * 0) Gather all FFT boxes from all tasks
   * ------------------------------------------------------------------ */

  int *all_boxes = (int *)malloc(6 * ntasks * sizeof(int));
  if (all_boxes == NULL)
  {
    printf("ERROR on task %d: could not allocate all_boxes\n", This);
    fflush(stdout);
    free(back_buffer);
    return 1;
  }

  MPI_Allgather(my_fft_box, 6, MPI_INT,
                all_boxes, 6, MPI_INT,
                MPI_COMM_WORLD);

  /* ------------------------------------------------------------------
   * 1) Count how many back_data entries go to each destination
   * ------------------------------------------------------------------ */

  int *sendcounts = (int *)calloc(ntasks, sizeof(int));
  int *recvcounts = (int *)calloc(ntasks, sizeof(int));
  int *sdispls = (int *)calloc(ntasks, sizeof(int));
  int *rdispls = (int *)calloc(ntasks, sizeof(int));

  if (!sendcounts || !recvcounts || !sdispls || !rdispls)
  {
    printf("ERROR on task %d: could not allocate count/displs arrays\n", This);
    fflush(stdout);
    free(back_buffer);
    free(all_boxes);
    free(sendcounts);
    free(recvcounts);
    free(sdispls);
    free(rdispls);
    return 1;
  }

  int iz, ibox, jbox, kbox;

  for (iz = 0; iz < subbox.Nstored; iz++)
  {
#ifdef CLASSIC_FRAGMENTATION
    INDEX_TO_COORD(iz, ibox, jbox, kbox, subbox.Lgwbl);
#else
    INDEX_TO_COORD(frag_pos[iz], ibox, jbox, kbox, subbox.Lgwbl);
#endif

    int core_particle =
        (ibox >= subbox.safe[_x_] && ibox < subbox.Lgwbl[_x_] - subbox.safe[_x_] &&
         jbox >= subbox.safe[_y_] && jbox < subbox.Lgwbl[_y_] - subbox.safe[_y_] &&
         kbox >= subbox.safe[_z_] && kbox < subbox.Lgwbl[_z_] - subbox.safe[_z_]);

#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
    int good_particle = (core_particle || frag[iz].zacc > (PRODFLOAT)-0.5);
#else
    int good_particle = core_particle;
#endif

    /* Global box frame, as in send_data_back */
    ibox = (ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_]) % MyGrids[0].GSglobal[_x_];
    jbox = (jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_]) % MyGrids[0].GSglobal[_y_];
    kbox = (kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_]) % MyGrids[0].GSglobal[_z_];

    if (!good_particle)
      continue;

    /* Find destination rank whose FFT box contains this cell */
    for (int dest = 0; dest < NTasks; dest++)
    {
      int *box = &all_boxes[6 * dest];

      if (ibox >= box[0] && ibox < box[0] + box[3] &&
          jbox >= box[1] && jbox < box[1] + box[4] &&
          kbox >= box[2] && kbox < box[2] + box[5])
      {
        sendcounts[dest] += 1;
        /* Assume non-overlapping domain, so break after first match */
        break;
      }
    }
  }

  /* ------------------------------------------------------------------
   * 2) Exchange counts so every rank knows recvcounts
   * ------------------------------------------------------------------ */

  MPI_Alltoall(sendcounts, 1, MPI_INT,
               recvcounts, 1, MPI_INT,
               MPI_COMM_WORLD);

  /* ------------------------------------------------------------------
   * 3) Build displacements and allocate send/recv buffers
   * ------------------------------------------------------------------ */

  int total_send = 0;
  int total_recv = 0;

  for (int r = 0; r < NTasks; r++)
  {
    sdispls[r] = total_send;
    total_send += sendcounts[r];

    rdispls[r] = total_recv;
    total_recv += recvcounts[r];
  }

  back_data *sendbuf = NULL;
  back_data *recvbuf = NULL;

  if (total_send > 0)
  {
    sendbuf = (back_data *)malloc((size_t)total_send * sizeof(back_data));
    if (sendbuf == NULL)
    {
      printf("ERROR on task %d: could not allocate sendbuf\n", This);
      fflush(stdout);
      free(back_buffer);
      free(all_boxes);
      free(sendcounts);
      free(recvcounts);
      free(sdispls);
      free(rdispls);
      return 1;
    }
  }

  if (total_recv > 0)
  {
    recvbuf = (back_data *)malloc((size_t)total_recv * sizeof(back_data));
    if (recvbuf == NULL)
    {
      printf("ERROR on task %d: could not allocate recvbuf\n", This);
      fflush(stdout);
      free(back_buffer);
      free(all_boxes);
      free(sendbuf);
      free(sendcounts);
      free(recvcounts);
      free(sdispls);
      free(rdispls);
      return 1;
    }
  }

  /* ------------------------------------------------------------------
   * 4) Pack data into sendbuf grouped by destination
   * ------------------------------------------------------------------ */

  int *offset = (int *)calloc(ntasks, sizeof(int));
  if (!offset)
  {
    printf("ERROR on task %d: could not allocate offset array\n", This);
    fflush(stdout);
    free(back_buffer);
    free(all_boxes);
    free(sendbuf);
    free(recvbuf);
    free(sendcounts);
    free(recvcounts);
    free(sdispls);
    free(rdispls);
    return 1;
  }

  for (int r = 0; r < NTasks; r++)
    offset[r] = sdispls[r];

  for (iz = 0; iz < subbox.Nstored; iz++)
  {
#ifdef CLASSIC_FRAGMENTATION
    INDEX_TO_COORD(iz, ibox, jbox, kbox, subbox.Lgwbl);
#else
    INDEX_TO_COORD(frag_pos[iz], ibox, jbox, kbox, subbox.Lgwbl);
#endif

    int core_particle =
        (ibox >= subbox.safe[_x_] && ibox < subbox.Lgwbl[_x_] - subbox.safe[_x_] &&
         jbox >= subbox.safe[_y_] && jbox < subbox.Lgwbl[_y_] - subbox.safe[_y_] &&
         kbox >= subbox.safe[_z_] && kbox < subbox.Lgwbl[_z_] - subbox.safe[_z_]);

#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
    int good_particle = (core_particle || frag[iz].zacc > (PRODFLOAT)-0.5);
#else
    int good_particle = core_particle;
#endif

    ibox = (ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_]) % MyGrids[0].GSglobal[_x_];
    jbox = (jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_]) % MyGrids[0].GSglobal[_y_];
    kbox = (kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_]) % MyGrids[0].GSglobal[_z_];

    if (!good_particle)
      continue;

    for (int dest = 0; dest < NTasks; dest++)
    {
      int *box = &all_boxes[6 * dest];

      if (ibox >= box[0] && ibox < box[0] + box[3] &&
          jbox >= box[1] && jbox < box[1] + box[4] &&
          kbox >= box[2] && kbox < box[2] + box[5])
      {
        int pos = COORD_TO_INDEX(ibox - box[0],
                                 jbox - box[1],
                                 kbox - box[2],
                                 (box + 3));

        int idx = offset[dest]++;

        sendbuf[idx].pos = pos;
        sendbuf[idx].zacc = frag[iz].zacc;
        sendbuf[idx].group_ID = frag[iz].group_ID;

        /* Assume non-overlapping domain */
        break;
      }
    }
  }

  free(offset);
  free(all_boxes);

  /* ------------------------------------------------------------------
   * 5) Perform the all-to-all variable-size exchange
   *    Use an MPI datatype matching back_data so counts/displs are in elements
   * ------------------------------------------------------------------ */

  MPI_Datatype back_data_type;
  MPI_Type_contiguous((int)sizeof(back_data), MPI_BYTE, &back_data_type);
  MPI_Type_commit(&back_data_type);

  MPI_Alltoallv(sendbuf, sendcounts, sdispls, back_data_type,
                recvbuf, recvcounts, rdispls, back_data_type,
                MPI_COMM_WORLD);

  MPI_Type_free(&back_data_type);

  /* ------------------------------------------------------------------
   * 6) Consume received data: same as recv_data_back inner loop
   * ------------------------------------------------------------------ */

  for (int i = 0; i < total_recv; i++)
  {
    unsigned int pos = recvbuf[i].pos;

    /* Existing local state (may have been filled by keep_data_back) */
    PRODFLOAT local_zacc = products[pos].zacc;
    int local_gid = products[pos].group_ID;

    /* Incoming state from another task */
    PRODFLOAT incoming_zacc = recvbuf[i].zacc;
    int incoming_gid = recvbuf[i].group_ID;

    /* ZACC:
       - Do not allow an incoming "unset" value (<= -0.5) to erase a local
         collapse tag (> -0.5).
       - Only fill a collapse time if the local value is still unset. */
    if (local_zacc <= (PRODFLOAT)-0.5 && incoming_zacc > (PRODFLOAT)-0.5)
      products[pos].zacc = incoming_zacc;

    /* group_ID:
       - Do not demote a local halo membership (>FILAMENT) to FILAMENT/none
         (<=FILAMENT) based on a ghost view.
       - Otherwise accept the incoming value. */
    if (!(local_gid > FILAMENT && incoming_gid <= FILAMENT))
      products[pos].group_ID = incoming_gid;
  }

  /* ------------------------------------------------------------------
   * Cleanup
   * ------------------------------------------------------------------ */

  free(recvbuf);
  free(sendbuf);
  free(sendcounts);
  free(recvcounts);
  free(sdispls);
  free(rdispls);

  free(back_buffer);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}

int my_distribute_back(void)
{
  /* Distributes accretion times from sub-volumes to fft-space */

  int my_fft_box[6];
  int This = ThisTask;

  /* this defines the box that the task possesses in the FFT space */
  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  /* synchronize the tasks here (kept for drop-in compatibility) */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* store the data relative to the intersection of the fft box and subbox */
  if (keep_data_back(my_fft_box))
    return 1;

  back_buffer = (back_data *)calloc(BUFLEN, sizeof(back_data));
  if (back_buffer == NULL)
  {
    printf("ERROR on task %d: could not allocate back_buffer in my_distribute_back\n", This);
    fflush(stdout);
    return 1;
  }

  /* ---- pairwise all-to-all schedule via round-robin algorithm ----
   *
   * We build a "tournament" schedule over total_teams slots.
   * If NTasks is odd, we add a dummy team that never corresponds to a real rank.
   * In each round, every real rank is paired with at most one other real rank.
   * Across all rounds, each unordered pair (i, j) appears exactly once.
   */

  int total_teams = NTasks;
  int have_dummy = 0;
  int dummy_rank = -1;

  if (total_teams % 2 == 1)
  {
    /* add a dummy slot */
    total_teams++;
    have_dummy = 1;
    dummy_rank = NTasks; /* no real task has this rank */
  }

  int *team = (int *)malloc(total_teams * sizeof(int));
  if (team == NULL)
  {
    printf("ERROR on task %d: could not allocate team array in my_distribute_back\n", This);
    fflush(stdout);
    free(back_buffer);
    return 1;
  }

  /* initial arrangement: 0..NTasks-1, plus dummy if needed */
  for (int i = 0; i < NTasks; i++)
    team[i] = i;

  if (have_dummy)
    team[total_teams - 1] = dummy_rank;
  else
    team[total_teams - 1] = total_teams - 1; /* last is also a real rank */

  int n_rounds = total_teams - 1; /* standard for round-robin */

  for (int round = 0; round < n_rounds; round++)
  {
    int partner = -1;

    /* determine my partner in this round, if any */
    for (int i = 0; i < total_teams / 2; i++)
    {
      int a = team[i];
      int b = team[total_teams - 1 - i];

      /* skip any pair involving the dummy slot */
      if (have_dummy && (a == dummy_rank || b == dummy_rank))
        continue;

      if (a == This)
      {
        partner = b;
        break;
      }
      else if (b == This)
      {
        partner = a;
        break;
      }
    }

    if (partner >= 0 && partner < NTasks)
    {
      /* symmetric two-way exchange between This and partner
         order chosen by rank to avoid deadlocks with blocking MPI calls */

      if (This < partner)
      {
        /* first send, then receive */
        if (send_data_back(partner))
        {
          free(team);
          free(back_buffer);
          return 1;
        }

        if (recv_data_back(my_fft_box, partner))
        {
          free(team);
          free(back_buffer);
          return 1;
        }
      }
      else
      {
        /* first receive, then send */
        if (recv_data_back(my_fft_box, partner))
        {
          free(team);
          free(back_buffer);
          return 1;
        }

        if (send_data_back(partner))
        {
          free(team);
          free(back_buffer);
          return 1;
        }
      }
    }

    /* rotate teams (round-robin "circle method"):
       keep team[0] fixed, rotate team[1..total_teams-1] to the right by one */
    int last = team[total_teams - 1];
    for (int i = total_teams - 1; i > 1; i--)
      team[i] = team[i - 1];
    team[1] = last;
  }

  free(team);
  free(back_buffer);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}

int distribute_back(void)
{
  /* Distributes accretion times from sub-volumes to fft-space */

  int my_fft_box[6];
  int log_ntask, bit, receiver, sender;

  /* this defines the box that the task possesses in the FFT space */
  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  /* let's synchronize the tasks here */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* stores the data relative to the intersection of the fft box and subbox */
  if (keep_data_back(my_fft_box))
    return 1;

  back_buffer = (back_data *)calloc(BUFLEN, sizeof(back_data));
  if (back_buffer == 0x0)
  {
    printf("ERROR on task %d: could not allocate back_buffer in distribute_back\n", ThisTask);
    fflush(stdout);
    return 1;
  }

  /* hypercubic communication scheme */
  for (log_ntask = 0; log_ntask < 1000; log_ntask++)
    if (1 << log_ntask >= NTasks)
      break;

  /* loop on hypercube dimension */
  for (bit = 1; bit < 1 << log_ntask; bit++)
  {
    /* loop on tasks */
    for (sender = 0; sender < NTasks; sender++)
    {
      /* receiver task is computed with a bitwise xor */
      receiver = sender ^ bit;

      /* condition on sender and receiver */
      if (receiver < NTasks && sender < receiver)
      {

        /* the communication will be done
     first sender -> receiver and then receiver -> sender */

        if (ThisTask == sender)
          send_data_back(receiver);
        else if (ThisTask == receiver)
        {
          if (recv_data_back(my_fft_box, sender))
            return 1;
        }

        if (ThisTask == receiver)
          send_data_back(sender);
        else if (ThisTask == sender)
        {
          if (recv_data_back(my_fft_box, receiver))
            return 1;
        }
      }
    }
  }

  free(back_buffer);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}

int keep_data_back(int *fft_box)
{
  /* this routine transfers products from the fft space to the subbox space */

  unsigned int ibox, jbox, kbox, good_particle, fftpos;

  /* find particles that are wanted by the receiver */
  for (int iz = 0; iz < subbox.Nstored; iz++)
  {
#ifdef CLASSIC_FRAGMENTATION
    INDEX_TO_COORD(iz, ibox, jbox, kbox, subbox.Lgwbl);
#else
    INDEX_TO_COORD(frag_pos[iz], ibox, jbox, kbox, subbox.Lgwbl);
#endif

    /* this is still in the subbox frame */
    unsigned int core_particle = (ibox >= subbox.safe[_x_] && ibox < subbox.Lgwbl[_x_] - subbox.safe[_x_] &&
                                  jbox >= subbox.safe[_y_] && jbox < subbox.Lgwbl[_y_] - subbox.safe[_y_] &&
                                  kbox >= subbox.safe[_z_] && kbox < subbox.Lgwbl[_z_] - subbox.safe[_z_]);

#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
    /* For MASS_MAPS_FILTER_UNCOLLAPSED we also need collapse times for halo
       particles that lie in the boundary layer; otherwise only the subset
       inside the core region would be filtered out of the mass maps while
       PLC catalogs use full halo masses, causing small mass mismatches when
       BoundaryLayerFactor > 0.  Treat any fragment with a valid zacc as
       eligible, even if outside the core. */
    good_particle = (core_particle || frag[iz].zacc > (PRODFLOAT)-0.5);
#else
    good_particle = core_particle;
#endif

    /* global box frame */
    ibox = (ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_]) % MyGrids[0].GSglobal[_x_];
    jbox = (jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_]) % MyGrids[0].GSglobal[_y_];
    kbox = (kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_]) % MyGrids[0].GSglobal[_z_];

    if (good_particle &&
        ibox >= fft_box[0] && ibox < fft_box[0] + fft_box[3] &&
        jbox >= fft_box[1] && jbox < fft_box[1] + fft_box[4] &&
        kbox >= fft_box[2] && kbox < fft_box[2] + fft_box[5])
    {
      /* position in the fft domain */
      fftpos = COORD_TO_INDEX(ibox - fft_box[0], jbox - fft_box[1], kbox - fft_box[2], (fft_box + 3));
      products[fftpos].zacc = frag[iz].zacc;
      products[fftpos].group_ID = frag[iz].group_ID;
    }
  }

  return 0;
}

int send_data_back(int target)
{
  /* This routine sends the content of a box to a target task.
     Communication is divided in these stages:
     1) the sender communicates the start and length of its box,
     2) the sender receives the start and the length of the needed box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the sender receives a map of needed particles
     5) the sender loops on particles and sends them in chunks of size BUFLEN
  */

  int ibox, jbox, kbox, good_particle;
  int recv_box[6];
  MPI_Status status;

  MPI_Recv(recv_box, 6, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

  int bufcount = 0;
  int sent = 0;

  /* find particles that are wanted by the receiver */
  for (int iz = 0; iz < subbox.Nstored; iz++)
  {
#ifdef CLASSIC_FRAGMENTATION
    INDEX_TO_COORD(iz, ibox, jbox, kbox, subbox.Lgwbl);
#else
    INDEX_TO_COORD(frag_pos[iz], ibox, jbox, kbox, subbox.Lgwbl);
#endif
    unsigned int core_particle = (ibox >= subbox.safe[_x_] && ibox < subbox.Lgwbl[_x_] - subbox.safe[_x_] &&
                                  jbox >= subbox.safe[_y_] && jbox < subbox.Lgwbl[_y_] - subbox.safe[_y_] &&
                                  kbox >= subbox.safe[_z_] && kbox < subbox.Lgwbl[_z_] - subbox.safe[_z_]);

#ifdef MASS_MAPS_FILTER_UNCOLLAPSED
    good_particle = (core_particle || frag[iz].zacc > (PRODFLOAT)-0.5);
#else
    good_particle = core_particle;
#endif

    /* global box frame */
    ibox = (ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_]) % MyGrids[0].GSglobal[_x_];
    jbox = (jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_]) % MyGrids[0].GSglobal[_y_];
    kbox = (kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_]) % MyGrids[0].GSglobal[_z_];

    if (good_particle &&
        ibox >= recv_box[0] && ibox < recv_box[0] + recv_box[3] &&
        jbox >= recv_box[1] && jbox < recv_box[1] + recv_box[4] &&
        kbox >= recv_box[2] && kbox < recv_box[2] + recv_box[5])
    {
      /* position in the fft domain */
      back_buffer[bufcount].pos = COORD_TO_INDEX(ibox - recv_box[0], jbox - recv_box[1], kbox - recv_box[2], (recv_box + 3));
      back_buffer[bufcount].zacc = frag[iz].zacc;
      back_buffer[bufcount].group_ID = frag[iz].group_ID;

      ++bufcount;
      ++sent;

      if (bufcount == BUFLEN)
      {
        MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
        MPI_Send(back_buffer, bufcount * sizeof(back_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
        bufcount = 0;
      }
    }
  }

  if (bufcount)
  {
    MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
    MPI_Send(back_buffer, bufcount * sizeof(back_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
  }

  bufcount = 0;
  MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);

  /* done */
  return 0;
}

int recv_data_back(int *mybox, int sender)
{
  /* This routine receives the content of a box from a sender task.
     Communication is divided in these stages:
     1) the target receives the start and length of sender box,
     2) the target sends the start and the length of its box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the target sends a map of needed particles
     5) the target loops on particles and receives them in chunks of size BUFLEN
  */

  unsigned int nsent, received;
  MPI_Status status;

  MPI_Send(mybox, 6, MPI_INT, sender, 0, MPI_COMM_WORLD);

  /* receive particles */
  received = 0;

  do
  {
    MPI_Recv(&nsent, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
    if (nsent)
      MPI_Recv(back_buffer, nsent * sizeof(back_data), MPI_BYTE, sender, 0, MPI_COMM_WORLD, &status);
    for (int iz = 0; iz < nsent; iz++)
    {
      unsigned int pos = back_buffer[iz].pos;

      PRODFLOAT local_zacc = products[pos].zacc;
      int local_gid = products[pos].group_ID;

      PRODFLOAT incoming_zacc = back_buffer[iz].zacc;
      int incoming_gid = back_buffer[iz].group_ID;

      /* ZACC: preserve any existing collapse tag; only fill if unset. */
      if (local_zacc <= (PRODFLOAT)-0.5 && incoming_zacc > (PRODFLOAT)-0.5)
        products[pos].zacc = incoming_zacc;

      /* group_ID: avoid demoting an existing halo membership to FILAMENT/none. */
      if (!(local_gid > FILAMENT && incoming_gid <= FILAMENT))
        products[pos].group_ID = incoming_gid;
    }
    received += nsent;
  } while (nsent);

  /* done */
  return 0;
}

#endif
