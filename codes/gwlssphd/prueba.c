#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>

#include "_gadgetconfig.h"
#include "allvars.h"
#include "fof.h"
#include "io.h"
#include "proto.h"
#include "version.h"
#include "subfind.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
#include "hdf5_util.h"
#endif

#define SELECT_MAIN_BRANCHES
#define PH_DECOMPOSITION 1
#define RANDOM_DECOMPOSITION 0
#define HASHTABLE_BITS_PER_DIMENSION 7

#if defined(FOF) && defined(SUBFIND)

static int GhostRegionThickness;
static int Ncells_per_iter;
static int Niter;
static int LastSnapShotNr;
static int LocNsubgroups;
int *TreeFirstHalo;
int *FileMap, *FilestoLoad;
int HashTabSize, NFiles, NumFilestoLoad;

double DomainCenterTree[3], DomainRadius;
extern struct struct_cell
{
  int Start;
  int Last;
} *Cell;

struct tree_data
{
  /* merger tree pointers */
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  /* properties of halo */
  int Len;
#ifndef LITE_MERGERTREE
  float M_Mean200, M_Crit200, M_TopHat200;
  float Pos[3];
  float Vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  // MODIFICACIÓN PARA L-GalaxiesBH: Variable CM eliminada 
  // float CM[3];
#endif

  MyIDType SubMostBoundID;
  /* original position in subfind output */
  int SnapNum, FileNr, SubIndex;
#ifndef LITE_MERGERTREE
  float RadMhalf;
#endif
};


struct halo_data
{
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  /* properties of halo */

  int Len;
  float Pos[3];
#ifndef LITE_MERGERTREE
  float M_Mean200, M_Crit200, M_TopHat200;
  float Vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  float CM[3];
#endif
  MyIDType SubMostBoundID;

  /* original position in subfind output */
  int SnapNum, FileNr, SubIndex;
  float RadMhalf;
  long long FoF_GrNr;

  /* descendant information */
  int Desc_FileNr, Desc_SubIndex, Desc_GrNr;
  float Desc_Weight;

  int TargetIndex;
  char UsedFlag, HaloFlag, GhostFlag;

}
 *Halo, *HaloList;



struct struct_list
{
  int FileNr;
  int next;
} *List;



struct haloaux_data
{
  long long MassInStem;
  int DoneFlag;
} *HaloAux;


#ifdef SELECT_MAIN_BRANCHES
static int count, countalt;
static int count_fof, countalt_fof;
#endif

static int NTrees, MaxNTrees;   /* number of trees, and allocated maximum */
static int NHalos, MaxNHalos;   /* number of halos in trees, and allocated maximum */
static int *TreeLen;            /* Table that gives the length of each tree */
static struct tree_data *TreeHalos, *ProcessHalo;       /* data for tree halos */
static int MaxTreeLen;          /* gives the maximum allocated tree length */
static int CountUsed;

static int TreeNimport, TreeNexport;

void set_desc_grnr(void);
int ShouldWeLoadCell(int iter, float x, float y, float z, double domainFac);
void tree_domain_decompose(int iter, int type);
int build_mergertree_single(int finalsnap, int iter);
int comp_halo_order(const void *a, const void *b);
int comp_halo_order_grnr(const void *a, const void *b);
void tree_load_subhalo_catalogue(int num, int lastsnap);
int get_target_task(int iter, float *pos, int type, char typePart);
void get_ghost_halos(void);
void walk_it(int i);
int check_tree_candidate(int iter);
void generate_trees(int iter);
void save_trees(int iter);
void set_descendant_pointers(void);
void set_progenitor_pointers(void);
void set_fof_pointers(void);

#ifdef SELECT_MAIN_BRANCHES
void process_tree(int tree);
void walk(int nr);
void modify_tree(int tree);
#endif


/* 
 * build_mergertree: Initializes spatial hash table grid parameters or 
 * single-iteration configurations, sets the ghost region thickness, 
 * and sequentially drives the execution of the merger tree building
 * routine across all defined iterations.
 */

void build_mergertree(int finalsnap)
{
  int total_number_of_cells = (((peanokey) 1) << HASHTABLE_BITS_PER_DIMENSION);
#ifdef HASHTABLE
  Ncells_per_iter = CUB(1);
  Niter = CUB(HASHTABLE_BITS_PER_DIMENSION) / Ncells_per_iter;
#else
  Niter = 1;
  Ncells_per_iter = CUB(total_number_of_cells);
#endif
  mpi_printf("number of cells per dimension: %d \n", total_number_of_cells);

  GhostRegionThickness = 10.;

  for(int iter = 0; iter < Niter; iter++)
    build_mergertree_single(finalsnap, iter);
}


/* 
 * build_mergertree_single: Executes a single iteration of the merger 
 * tree building pipeline, loading catalogues across snapshots,
 * performing domain decompositions and ghost halo exchanges, linking 
 * descendant, progenitor, and FOF structures, optional main branch 
 * corrections, and saving the final trees.
 */

int build_mergertree_single(int finalsnap, int iter)
{
  LastSnapShotNr = finalsnap;
  LocNsubgroups = 0;
  Halo = NULL;

  DomainFac = 1.0 / All.BoxSize * (((peanokey) 1) << (HASHTABLE_BITS_PER_DIMENSION));

  Send_count = (int *) mymalloc("Send_count", sizeof(int) * NTask);
  Send_offset = (int *) mymalloc("Send_offset", sizeof(int) * NTask);
  Recv_count = (int *) mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc("Recv_offset", sizeof(int) * NTask);

  mpi_printf("loading halo catalogues and descendant tree files...\n");

  for(int num = finalsnap; num >= 0; num--)
    {
      tree_load_subhalo_catalogue(num, finalsnap);
      tree_domain_decompose(iter, RANDOM_DECOMPOSITION);

      int maxload, minload, totload, totghost, LocNghost = 0;

      for(int i = 0; i < LocNsubgroups; i++)
        if(Halo[i].GhostFlag == 1)
          LocNghost++;

      MPI_Reduce(&LocNsubgroups, &maxload, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&LocNsubgroups, &minload, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&LocNsubgroups, &totload, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&LocNghost, &totghost, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

      mpi_printf(" ----------------------------------------------\n");
      mpi_printf(" Snap:=%d max:=%d  min:=%d total=%d\n", num, maxload, minload, totload);
      mpi_printf(" Ghost =%d \n", totghost);
      mpi_printf(" --------------------------------------------------\n");
    }

#ifdef HASHTABLE
  Halo = myrealloc_movable(Halo, LocNsubgroups * sizeof(struct halo_data));
#endif

  int minfile = 10000;
  int maxfile = -20;

  for(int i = 0; i < LocNsubgroups; i++)
    {
      if(minfile > Halo[i].Desc_FileNr)
        minfile = Halo[i].Desc_FileNr;

      if(maxfile < Halo[i].Desc_FileNr)
        maxfile = Halo[i].Desc_FileNr;
    }

  mpi_printf("loading done.\n");
  mpi_printf("ThisTask=%d minmax FileNr %d, %d \n", ThisTask, minfile, maxfile);

  tree_domain_decompose(iter, PH_DECOMPOSITION);

  printf("Task %d Number of haloes: %d \n", ThisTask, LocNsubgroups);

  get_ghost_halos();

  for(int i = 0; i < LocNsubgroups; i++)
    {
      Halo[i].GhostFlag = 1;
      if(ShouldWeLoadCell(iter, Halo[i].Pos[0], Halo[i].Pos[1], Halo[i].Pos[2], DomainFac) == 1
         && ThisTask == get_target_task(iter, Halo[i].Pos, 1, 0))
        Halo[i].GhostFlag = 0;
    }

  MPI_Barrier(MPI_COMM_WORLD);
  mpi_printf("subhalos distributed and ghost regions imported.\n");

  set_desc_grnr();

  MPI_Barrier(MPI_COMM_WORLD);
  mpi_printf("sorting of subhalos done. Start to find descendant pointers.\n");

  qsort(Halo, LocNsubgroups, sizeof(struct halo_data), comp_halo_order_grnr);

  set_descendant_pointers();

  int totndesc, ndesc = 0;

  for(int i = 0; i < LocNsubgroups; i++)
    if(Halo[i].Descendant != -1 && Halo[i].GhostFlag == 0)
      ndesc++;

  MPI_Reduce(&ndesc, &totndesc, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("decendant pointers found for %d halos\n", totndesc);
  set_progenitor_pointers();

  mpi_printf("progenitor pointers set\n");
  set_fof_pointers();

  mpi_printf("fof pointers done.\n");
  generate_trees(iter);

#ifdef SELECT_MAIN_BRANCHES
  mpi_printf("selecting main branches.\n");

  int glob_count, glob_countalt, glob_count_fof, glob_countalt_fof;

  CountUsed = 0;
  count = countalt = count_fof = countalt_fof = 0;

  TreeFirstHalo = mymalloc("TreeFirstHalo", sizeof(int) * NTrees);

  if(NTrees)
    TreeFirstHalo[0] = 0;

  for(int i = 1; i < NTrees; i++)
    {
      TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeLen[i - 1];
    }

  for(int tree = 0; tree < NTrees; tree++)
    {
      if((tree % 10000) == 0)
        printf("%d | %d\n", tree, NTrees);

      ProcessHalo = mymalloc("ProcessHalo", sizeof(struct tree_data) * TreeLen[tree]);
      memcpy(ProcessHalo, &TreeHalos[TreeFirstHalo[tree]], sizeof(struct tree_data) * TreeLen[tree]);

      int *done = mymalloc("done", sizeof(int) * TreeLen[tree]);
      HaloAux = mymalloc("HaloAux", sizeof(struct haloaux_data) * TreeLen[tree]);

      process_tree(tree);
      modify_tree(tree);

      memcpy(&TreeHalos[TreeFirstHalo[tree]], ProcessHalo, sizeof(struct tree_data) * TreeLen[tree]);

      myfree(HaloAux);
      myfree(done);
      myfree(ProcessHalo);
    }
  myfree(TreeFirstHalo);

  MPI_Reduce(&count, &glob_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&countalt, &glob_countalt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&count_fof, &glob_count_fof, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&countalt_fof, &glob_countalt_fof, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("count=%d countalt=%d  fraction=%g\n", glob_count, glob_countalt,
             ((double) glob_countalt) / glob_count);
  mpi_printf("FOF count=%d countalt=%d  fraction=%g\n", glob_count_fof, glob_countalt_fof,
             ((double) glob_countalt_fof) / glob_count_fof);
#endif

  save_trees(iter);

  return 0;
}

void set_desc_grnr()
{
  MPI_Barrier(MPI_COMM_WORLD);
  mpi_printf("Setting the Desc_GrNr key in subhalos catalogues\n");

  qsort(Halo, LocNsubgroups, sizeof(struct halo_data), comp_halo_order);
  set_descendant_pointers();

  for (int i =0; i < LocNsubgroups; i++)
      if (Halo[i].Descendant != -1)
        Halo[i].Desc_GrNr = Halo[Halo[i].Descendant].FoF_GrNr;
}


#ifdef SELECT_MAIN_BRANCHES
void process_tree(int tree)
{
  for(int halonr = 0; halonr < TreeLen[tree]; halonr++)
    {
      HaloAux[halonr].MassInStem = 0;
      HaloAux[halonr].DoneFlag = 0;
    }

  for(int halonr = 0; halonr < TreeLen[tree]; halonr++)
    {
      if(HaloAux[halonr].MassInStem == 0)
        walk(halonr);
    }
}

void walk(int nr)
{
  int first, next;
  long long maxstem;
  maxstem = 0;

  if(ProcessHalo[nr].FirstProgenitor >= 0)
    {
      first = ProcessHalo[nr].FirstProgenitor;

      if(HaloAux[first].MassInStem == 0)
        walk(first);

      if(maxstem < HaloAux[first].MassInStem)
        maxstem = HaloAux[first].MassInStem;

      next = ProcessHalo[first].NextProgenitor;

      while(next >= 0)
        {
          if(HaloAux[next].MassInStem == 0)
            walk(next);

          if(maxstem < HaloAux[next].MassInStem)
            maxstem = HaloAux[next].MassInStem;

          next = ProcessHalo[next].NextProgenitor;
        }
    }

  maxstem += ProcessHalo[nr].Len;
  HaloAux[nr].MassInStem = maxstem;
}

void modify_tree(int tree)
{
  int first, next, alternate_first;
  int previous, alternate_previous=-1;
  long long maxstem;

  for(int nr = 0; nr < TreeLen[tree]; nr++)
    {
      if(ProcessHalo[nr].FirstProgenitor >= 0)
        {
          count++;
          first = alternate_first = ProcessHalo[nr].FirstProgenitor;
          maxstem = HaloAux[first].MassInStem;
          previous = first;
          next = ProcessHalo[first].NextProgenitor;

          while(next >= 0)
            {
              if(maxstem < HaloAux[next].MassInStem)
                {
                  maxstem = HaloAux[next].MassInStem;
                  alternate_first = next;
                  alternate_previous = previous;
                }
              previous = next;
              next = ProcessHalo[next].NextProgenitor;
            }

          if(alternate_first != first)
            {
              countalt++;
              ProcessHalo[alternate_previous].NextProgenitor = ProcessHalo[alternate_first].NextProgenitor;
              ProcessHalo[alternate_first].NextProgenitor = ProcessHalo[nr].FirstProgenitor;
              ProcessHalo[nr].FirstProgenitor = alternate_first;
            }
        }
    }

  for(int nr = 0; nr < TreeLen[tree]; nr++)
    {
      if(ProcessHalo[nr].FirstHaloInFOFgroup == nr && HaloAux[nr].DoneFlag == 0)
        {
          count_fof++;
          first = alternate_first = ProcessHalo[nr].FirstHaloInFOFgroup;
          maxstem = HaloAux[first].MassInStem;
          previous = first;
          next = ProcessHalo[first].NextHaloInFOFgroup;

          while(next >= 0)
            {
              if(maxstem < HaloAux[next].MassInStem)
                {
                  maxstem = HaloAux[next].MassInStem;
                  alternate_first = next;
                  alternate_previous = previous;
                }
              previous = next;
              next = ProcessHalo[next].NextHaloInFOFgroup;
            }

          if(alternate_first != first)
            {
              countalt_fof++;

#ifndef LITE_MERGERTREE
              ProcessHalo[alternate_first].M_Mean200 = ProcessHalo[first].M_Mean200;
              ProcessHalo[alternate_first].M_Crit200 = ProcessHalo[first].M_Crit200;
              ProcessHalo[alternate_first].M_TopHat200 = ProcessHalo[first].M_TopHat200;
              ProcessHalo[first].M_Mean200 = 0;
              ProcessHalo[first].M_Crit200 = 0;
              ProcessHalo[first].M_TopHat200 = 0;
#endif

              next = first;
              while(next >= 0)
                {
                  ProcessHalo[next].FirstHaloInFOFgroup = alternate_first;
                  next = ProcessHalo[next].NextHaloInFOFgroup;
                }

              ProcessHalo[alternate_previous].NextHaloInFOFgroup =
                ProcessHalo[alternate_first].NextHaloInFOFgroup;
              ProcessHalo[alternate_first].NextHaloInFOFgroup = first;
              HaloAux[alternate_first].DoneFlag = 1;
            }
        }
    }
}
#endif

int comp_halo_order_grnr(const void *a, const void *b)
{
  if(((struct halo_data *) a)->SnapNum < ((struct halo_data *) b)->SnapNum)
    return -1;
  if(((struct halo_data *) a)->SnapNum > ((struct halo_data *) b)->SnapNum)
    return +1;
  if(((struct halo_data *) a)->FileNr < ((struct halo_data *) b)->FileNr)
    return -1;
  if(((struct halo_data *) a)->FileNr > ((struct halo_data *) b)->FileNr)
    return +1;
  if(((struct halo_data *) a)->FoF_GrNr < ((struct halo_data *) b)->FoF_GrNr)
    return -1;
  if(((struct halo_data *) a)->FoF_GrNr > ((struct halo_data *) b)->FoF_GrNr)
    return +1;
  if(((struct halo_data *) a)->SubIndex < ((struct halo_data *) b)->SubIndex)
    return -1;
  if(((struct halo_data *) a)->SubIndex > ((struct halo_data *) b)->SubIndex)
    return +1;
  return 0;
}

int comp_halo_order(const void *a, const void *b)
{
  if(((struct halo_data *) a)->SnapNum < ((struct halo_data *) b)->SnapNum)
    return -1;
  if(((struct halo_data *) a)->SnapNum > ((struct halo_data *) b)->SnapNum)
    return +1;
  if(((struct halo_data *) a)->FileNr < ((struct halo_data *) b)->FileNr)
    return -1;
  if(((struct halo_data *) a)->FileNr > ((struct halo_data *) b)->FileNr)
    return +1;
  if(((struct halo_data *) a)->SubIndex < ((struct halo_data *) b)->SubIndex)
    return -1;
  if(((struct halo_data *) a)->SubIndex > ((struct halo_data *) b)->SubIndex)
    return +1;
  return 0;
}

void tree_load_subhalo_catalogue(int num, int lastsnap)
{
#ifndef HASHTABLE
    int load_id_list = 0;
    fof_subfind_load_groups(num, load_id_list); 
    if (num != lastsnap)
        fof_subfind_load_descendants(num);

    if(Halo == NULL)
        Halo = mymalloc_movable(&Halo, "Halo", (LocNsubgroups + Nsubgroups) * sizeof(struct halo_data));
    else
        Halo = myrealloc_movable(Halo, (LocNsubgroups + Nsubgroups) * sizeof(struct halo_data));

    for(int i = 0; i < Nsubgroups; i++)
    {
      Halo[LocNsubgroups + i].SnapNum = num;
      Halo[LocNsubgroups + i].FileNr = ThisTask;
      Halo[LocNsubgroups + i].SubIndex = i;
      Halo[LocNsubgroups + i].FoF_GrNr = (int) SubGroup[i].GrNr;
      Halo[LocNsubgroups + i].Len = SubGroup[i].Len;

      for(int j = 0; j < 3; j++)
          Halo[LocNsubgroups + i].Pos[j] = SubGroup[i].Pos[j];
#ifndef LITE_MERGERTREE
      Halo[LocNsubgroups + i].GhostFlag = 0;
      Halo[LocNsubgroups + i].M_Crit200 = 0.0;
      Halo[LocNsubgroups + i].M_Mean200 = 0.0;
      Halo[LocNsubgroups + i].M_TopHat200 = 0.0;
      for(int j = 0; j < 3; j++)
        {
          Halo[LocNsubgroups + i].Vel[j] = SubGroup[i].Vel[j];
          Halo[LocNsubgroups + i].Spin[j] = SubGroup[i].Spin[j];
        }
      Halo[LocNsubgroups + i].VelDisp = SubGroup[i].SubVpeak;
      Halo[LocNsubgroups + i].Vmax = SubGroup[i].SubVmax;
#endif

      if (num == lastsnap)
      {
        Halo[LocNsubgroups + i].SubMostBoundID = SubGroup[i].SubMostBoundID;
        Halo[LocNsubgroups + i].Desc_FileNr = -1;
        Halo[LocNsubgroups + i].Desc_SubIndex = -1;
        Halo[LocNsubgroups + i].Desc_Weight = -1;
        Halo[LocNsubgroups + i].Desc_GrNr = -1;
      }
      else
      {
        Halo[LocNsubgroups + i].SubMostBoundID = CatA[i].SubMostBoundID;
        Halo[LocNsubgroups + i].Desc_FileNr = CatA[i].Desc_FileNr;
        Halo[LocNsubgroups + i].Desc_SubIndex = CatA[i].Desc_SubIndex;
        Halo[LocNsubgroups + i].Desc_Weight = CatA[i].Desc_Weight;
        Halo[LocNsubgroups + i].Desc_GrNr = -1;

        if(CatA[i].Desc_FileNr < -1 || CatA[i].Desc_FileNr >= NTask)
            Terminate("[task=%d] num=%d  desc_FileNr[%d]=%d, LastSnapShotNr=%d Nsubgroups=%d\n", ThisTask, num, i,
                  CatA[i].Desc_FileNr, LastSnapShotNr, Nsubgroups);
      }
    }

#ifndef LITE_MERGERTREE
    for(int i = 0; i < Ngroups; i++)
    {
      if(Group[i].Nsubs > 0)
        {
          Halo[LocNsubgroups + Group[i].FirstSub].M_Crit200 = Group[i].M_Crit200;
          Halo[LocNsubgroups + Group[i].FirstSub].M_Mean200 = Group[i].M_Mean200;
          Halo[LocNsubgroups + Group[i].FirstSub].M_TopHat200 = Group[i].M_TopHat200;
        }
    }
#endif
  LocNsubgroups = LocNsubgroups + Nsubgroups;
 
  if (num != lastsnap)
    myfree_movable(CatA);
  free_loaded_groups(load_id_list);
#else
  int i, j, k, nr, num_files, keep, hashkey, ip;
  double len, xmin[3], xmax[3], glob_xmin[3], glob_xmax[3];
  double dx, dy, dz, rr;

  if(ThisTask == 0)
    load_hash_table(num);

  if(num == LastSnapShotNr)
    {
      mark_hash_cells(iter, 1);
      for(int ii = ThisTask; ii < NumFilestoLoad; ii += NTask)
        {
          printf("Loading %d out of %d file %d\n", ii, NumFilestoLoad, FilestoLoad[ii]);
          load_subhalo_catalogue_single(num, FilestoLoad[ii]);
        }

      for(keep = 0, ip = 0; ip < LocNsubgroups; ip++)
        {
          if(ShouldWeLoadCell(iter, Halo[ip].Pos[0], Halo[ip].Pos[1], Halo[ip].Pos[2], DomainFac))
            {
              Halo[ip].GhostFlag = 0;
              memcpy(&Halo[keep], &Halo[ip], sizeof(struct halo_data));
              keep++;
            }
        }
      LocNsubgroups = keep;

      for(j = 0; j < 3; j++) { xmin[j] = 1.0e30; xmax[j] = -1.0e30; }
      for(i = 0; i < LocNsubgroups; i++)
        {
          for(j = 0; j < 3; j++)
            {
              if(xmin[j] > Halo[i].Pos[j]) xmin[j] = Halo[i].Pos[j];
              if(xmax[j] < Halo[i].Pos[j]) xmax[j] = Halo[i].Pos[j];
            }
        }

      MPI_Allreduce(xmin, glob_xmin, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(xmax, glob_xmax, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      len = 0;
      for(j = 0; j < 3; j++)
        if(glob_xmax[j] - glob_xmin[j] > len) len = glob_xmax[j] - glob_xmin[j];

      for(j = 0; j < 3; j++)
        DomainCenterTree[j] = 0.5 * (glob_xmin[j] + glob_xmax[j]);

      DomainRadius = 0.5 * len * sqrt(3);
      LocNsubgroups = 0;
    }

  mark_hash_cells(iter, 0);

  for(int ii = ThisTask; ii < NumFilestoLoad; ii += NTask)
    {
      printf("Loading %d out of %d file %d\n", ii, NumFilestoLoad, FilestoLoad[ii]);
      load_subhalo_catalogue_single(num, FilestoLoad[ii]);
    }

  for(keep = 0, ip = 0; ip < LocNsubgroups; ip++)
    {
      dx = periodic_wrap(Halo[ip].Pos[0] - DomainCenterTree[0]);
      dy = periodic_wrap(Halo[ip].Pos[1] - DomainCenterTree[1]);
      dz = periodic_wrap(Halo[ip].Pos[2] - DomainCenterTree[2]);

      rr = sqrt(dx * dx + dy * dy + dz * dz);
      hashkey = (i * HashTabSize + j) * HashTabSize + k;

      if(ShouldWeLoadCell(iter, Halo[ip].Pos[0], Halo[ip].Pos[1], Halo[ip].Pos[2], DomainFac))
        Halo[ip].GhostFlag = 0;

      if(rr < (DomainRadius + GhostRegionThickness + All.BoxSize / HashTabSize))
        {
          memcpy(&Halo[keep], &Halo[ip], sizeof(struct halo_data));
          keep++;
          if(num == LastSnapShotNr)
            Halo[ip].GhostFlag = 1;
        }
    }
  LocNsubgroups = keep;

  HighMarkBytes_load = 0;
  report_memory_usage(&HighMarkBytes_load, "LOAD");

  myfree_movable(FilestoLoad);

  if(ThisTask == 0)
    free_hash_table();
#endif
}

#ifdef HASHTABLE
void mark_hash_cells(int iter, int type)
{
  int i, j, k, next, hashkey, ncellsload, flag;
  double x, y, z, dx, dy, dz, rr;
  peanokey key;

  if(ThisTask == 0)
    {
      printf("Marking cells (HashTabSize = %d, BoxSize=%f, DomainFac=%f)\n", HashTabSize, All.BoxSize, DomainFac);
      ncellsload = 0;
      NumFilestoLoad = 0;
      for(i = 0; i < HashTabSize; i++)
        for(j = 0; j < HashTabSize; j++)
          for(k = 0; k < HashTabSize; k++)
            {
              x = ((i + 0.5) * All.BoxSize) / HashTabSize;
              y = ((j + 0.5) * All.BoxSize) / HashTabSize;
              z = ((k + 0.5) * All.BoxSize) / HashTabSize;

              dx = periodic_wrap(x - DomainCenterTree[0]);
              dy = periodic_wrap(y - DomainCenterTree[1]);
              dz = periodic_wrap(z - DomainCenterTree[2]);

              rr = sqrt(dx * dx + dy * dy + dz * dz);
              hashkey = (i * HashTabSize + j) * HashTabSize + k;

              flag = 0;
              if(type == 0)
                if(rr < (DomainRadius + GhostRegionThickness + All.BoxSize / HashTabSize)) flag = 1;

              if(type == 1)
                if(ShouldWeLoadCell(iter, x, y, z, DomainFac)) flag = 1;

              if(flag == 1)
                {
                  ncellsload++;
                  next = Cell[hashkey].Start;
                  while(next != -1)
                    {
                      FileMap[List[next].FileNr] = 1;
                      next = List[next].next;
                    }
                }
            }

      for(int i = 0; i < NFiles; i++)
        if(FileMap[i] == 1) NumFilestoLoad++;

      printf("Loading %d cells on %d files\n", ncellsload, NumFilestoLoad);
    }

  MPI_Bcast(&HashTabSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NumFilestoLoad, 1, MPI_INT, 0, MPI_COMM_WORLD);
  FilestoLoad = mymalloc_movable(&FilestoLoad, "FilestoLoad", sizeof(int) * NumFilestoLoad);

  int n = 0;
  if(ThisTask == 0)
    for(i = 0; i < NFiles; i++)
      if(FileMap[i] == 1) { FilestoLoad[n] = i; n++; }

  MPI_Bcast(FilestoLoad, NumFilestoLoad, MPI_INT, 0, MPI_COMM_WORLD);
}
#endif

int ShouldWeLoadCell(int iter, float x, float y, float z, double domainFac)
{
#ifndef HASHTABLE
  return 1;
#else
  peanokey key = peano_hilbert_key((int) (x * domainFac),
                                   (int) (y * domainFac),
                                   (int) (z * domainFac),
                                   HASHTABLE_BITS_PER_DIMENSION);

  if(key >= iter * Ncells_per_iter && key < (iter + 1) * Ncells_per_iter)
    return 1;

  return 0;
#endif
}

int get_target_task(int iter, float *pos, int type, char typePart)
{
  int task;
  if(type == PH_DECOMPOSITION)
    {
      if(typePart == 1)         
        {
          task = (int) floor((pos[0] - (int) pos[0]) * NTask);
        }
      else
        {
          peanokey key = peano_hilbert_key((int) (pos[0] * DomainFac),
                                           (int) (pos[1] * DomainFac),
                                           (int) (pos[2] * DomainFac),
                                           HASHTABLE_BITS_PER_DIMENSION);

          task = (key - Ncells_per_iter * iter) / ((float) (Ncells_per_iter) / NTask);
          if (task < 0) task = 0;
          if (task >= NTask) task = NTask-1;
        }
    }
  else if(type == RANDOM_DECOMPOSITION)
    {
      task = (int) floor((pos[0] - (int) pos[0]) * NTask);
    }
  else
  Terminate("get_target_task mode =%d not recognized\n", type) 
  if(task >= NTask)
    task = NTask - 1;
  assert(task >= 0 && task < NTask);
  return task;
}

void tree_domain_decompose(int iter, int type)
{
  int n, target, ngrp;
  long long i;

  for(n = 0; n < NTask; n++) { Send_count[n] = 0; Recv_count[n] = 0; }

  for(i = 0; i < LocNsubgroups; i++)
    {
      target = get_target_task(iter, Halo[i].Pos, type, Halo[i].GhostFlag);
      Send_count[target]++;
    }

  myMPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, TreeNexport = 0, TreeNimport = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      TreeNimport += Recv_count[i];
      TreeNexport += Send_count[i];
      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  for(n = 0; n < NTask; n++) Send_count[n] = 0;

  struct halo_data *Export_Halo_List =
    mymalloc_movable(&Export_Halo_List, "Export_Halo_List", TreeNexport * sizeof(struct halo_data));

  for(i = 0; i < LocNsubgroups; i++)
    {
      target = get_target_task(iter, Halo[i].Pos, type, Halo[i].GhostFlag);
      Export_Halo_List[Send_offset[target] + Send_count[target]] = Halo[i];
      Send_count[target]++;
    }

  LocNsubgroups = TreeNimport;

  if(Halo == NULL)
    Halo = mymalloc_movable(&Halo, "Halo", LocNsubgroups * sizeof(struct halo_data));
  else
    Halo = myrealloc_movable(Halo, LocNsubgroups * sizeof(struct halo_data));

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;
      if(target < NTask)
        {
          if(Send_count[target] > 0 || Recv_count[target] > 0)
            {
              MPI_Sendrecv(Export_Halo_List + Send_offset[target],
                           Send_count[target] * sizeof(struct halo_data), MPI_BYTE, target, TAG_PDATA,
                           Halo + Recv_offset[target],
                           Recv_count[target] * sizeof(struct halo_data), MPI_BYTE, target, TAG_PDATA,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
  myfree_movable(Export_Halo_List);
  MPI_Barrier(MPI_COMM_WORLD);
  mpi_printf("domain decomposition done.\n");
}

double periodic_wrap(double x)
{
  if(x > 0.5 * All.BoxSize) x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize) x += All.BoxSize;
  return x;
}

void get_ghost_halos(void)
{
  int i, j, ngrp;
  double len, xmin[3], xmax[3];

  for(j = 0; j < 3; j++) { xmin[j] = 1.0e30; xmax[j] = -1.0e30; }

  for(i = 0; i < LocNsubgroups; i++)
    {
      if(Halo[i].SnapNum == LastSnapShotNr)
        for(j = 0; j < 3; j++)
          {
            if(xmin[j] > Halo[i].Pos[j]) xmin[j] = Halo[i].Pos[j];
            if(xmax[j] < Halo[i].Pos[j]) xmax[j] = Halo[i].Pos[j];
          }
    }

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax[j] - xmin[j] > len) len = xmax[j] - xmin[j];

  double domainCenter[3], domainRadius;
  for(j = 0; j < 3; j++) domainCenter[j] = 0.5 * (xmin[j] + xmax[j]);
  domainRadius = 0.5 * len * sqrt(3);

  if(LocNsubgroups == 0) DomainRadius = -1;

  double *DomainCenterTable =
    mymalloc_movable(&DomainCenterTable, "DomainCenterTable", NTask * 3 * sizeof(double));
  double *DomainRadiusTable =
    mymalloc_movable(&DomainRadiusTable, "DomainRadiusTable", NTask * sizeof(double));

  MPI_Allgather(domainCenter, 3, MPI_DOUBLE, DomainCenterTable, 3, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(&domainRadius, 1, MPI_DOUBLE, DomainRadiusTable, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  struct halo_data *halobuf = mymalloc_movable(&halobuf, "halobuf", LocNsubgroups * sizeof(struct halo_data));
  int recv_off = LocNsubgroups;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int target = ThisTask ^ ngrp;
      if(target < NTask)
        {
          for(i = 0, TreeNexport = 0; i < LocNsubgroups; i++)
            {
              double dx = periodic_wrap(Halo[i].Pos[0] - DomainCenterTable[3 * target + 0]);
              double dy = periodic_wrap(Halo[i].Pos[1] - DomainCenterTable[3 * target + 1]);
              double dz = periodic_wrap(Halo[i].Pos[2] - DomainCenterTable[3 * target + 2]);
              double r = sqrt(dx * dx + dy * dy + dz * dz);

              if(DomainRadiusTable[target] >= 0)
                {
                  if(r < (DomainRadiusTable[target] + GhostRegionThickness))
                    {
                      halobuf[TreeNexport] = Halo[i];
                      TreeNexport++;
                    }
                }
            }

          MPI_Sendrecv(&TreeNexport, 1, MPI_INT, target, TAG_N,
                       &TreeNimport, 1, MPI_INT, target, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          if(Halo == NULL)
            Halo = mymalloc_movable(&Halo, "Halo", (recv_off + TreeNimport) * sizeof(struct halo_data));
          else
            Halo = myrealloc_movable(Halo, (recv_off + TreeNimport) * sizeof(struct halo_data));

          MPI_Sendrecv(halobuf, TreeNexport * sizeof(struct halo_data), MPI_BYTE, target, TAG_PDATA,
                       Halo + recv_off, TreeNimport * sizeof(struct halo_data), MPI_BYTE, target, TAG_PDATA,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          recv_off += TreeNimport;
        }
    }

  LocNsubgroups = recv_off;
  myfree(halobuf);
  myfree(DomainRadiusTable);
  myfree(DomainCenterTable);
}

void walk_it(int i)
{
  int p;
  Halo[i].UsedFlag = 1;
  Halo[i].TargetIndex = CountUsed;

  if(CountUsed >= MaxTreeLen)
    {
      MaxTreeLen *= 2;
      HaloList = myrealloc_movable(HaloList, MaxTreeLen * sizeof(struct halo_data));
    }

  HaloList[CountUsed] = Halo[i];
  CountUsed++;

  if(Halo[i].Descendant >= 0)
    if(Halo[Halo[i].Descendant].UsedFlag == 0)
      walk_it(Halo[i].Descendant);

  p = Halo[i].FirstProgenitor;
  while(p >= 0)
    {
      if(Halo[p].UsedFlag == 0) walk_it(p);
      p = Halo[p].NextProgenitor;
    }

  p = Halo[i].FirstHaloInFOFgroup;
  if(Halo[p].HaloFlag == 0)
    {
      Halo[p].HaloFlag = 1;
      while(p >= 0)
        {
          if(Halo[p].UsedFlag == 0) walk_it(p);
          p = Halo[p].NextHaloInFOFgroup;
        }
    }
}

int check_tree_candidate(int iter)
{
  int i, k, kmin, valid, target;
#ifndef LITE_MERGERTREE
  int j;
#endif

  for(k = 0, kmin = -1, valid = 0; k < CountUsed; k++)
    {
      if(HaloList[k].SnapNum == LastSnapShotNr)
        {
          target = get_target_task(iter, HaloList[k].Pos, 1, HaloList[k].GhostFlag);
          if(HaloList[k].GhostFlag == 0) valid = 1;          
          if(kmin < 0) kmin = k;
          else { if(HaloList[k].SubMostBoundID < HaloList[kmin].SubMostBoundID) kmin = k; }
        }
    }

  if(valid)
    {
      target = get_target_task(iter, HaloList[kmin].Pos, 1, HaloList[kmin].GhostFlag);
      if(HaloList[kmin].GhostFlag == 0)
        {
          if(NTrees >= MaxNTrees)
            {
              MaxNTrees *= 2;
              TreeLen = myrealloc_movable(TreeLen, MaxNTrees * sizeof(int));
            }

          if(NHalos + CountUsed > MaxNHalos)
            {
              MaxNHalos = NHalos + CountUsed + 10000;
              TreeHalos = myrealloc_movable(TreeHalos, MaxNHalos * sizeof(struct tree_data));
            }

          for(i = 0; i < CountUsed; i++)
            {
              TreeHalos[NHalos + i].Descendant = HaloList[i].Descendant;
              TreeHalos[NHalos + i].FirstProgenitor = HaloList[i].FirstProgenitor;
              TreeHalos[NHalos + i].NextProgenitor = HaloList[i].NextProgenitor;
              TreeHalos[NHalos + i].FirstHaloInFOFgroup = HaloList[i].FirstHaloInFOFgroup;
              TreeHalos[NHalos + i].NextHaloInFOFgroup = HaloList[i].NextHaloInFOFgroup;
              TreeHalos[NHalos + i].Len = HaloList[i].Len;

#ifndef LITE_MERGERTREE
              TreeHalos[NHalos + i].M_Mean200 = HaloList[i].M_Mean200;
              TreeHalos[NHalos + i].M_Crit200 = HaloList[i].M_Crit200;
              TreeHalos[NHalos + i].M_TopHat200 = HaloList[i].M_TopHat200;
              for(j = 0; j < 3; j++)
                {
                  TreeHalos[NHalos + i].Pos[j] = HaloList[i].Pos[j];
                  TreeHalos[NHalos + i].Vel[j] = HaloList[i].Vel[j];
                  TreeHalos[NHalos + i].Spin[j] = HaloList[i].Spin[j];
                }
              TreeHalos[NHalos + i].RadMhalf = HaloList[i].RadMhalf;
              TreeHalos[NHalos + i].Vmax = HaloList[i].Vmax;
              TreeHalos[NHalos + i].VelDisp = HaloList[i].VelDisp;
#endif
              TreeHalos[NHalos + i].SubMostBoundID = HaloList[i].SubMostBoundID;
              TreeHalos[NHalos + i].SnapNum = HaloList[i].SnapNum;
              TreeHalos[NHalos + i].FileNr = HaloList[i].FileNr;
              TreeHalos[NHalos + i].SubIndex = HaloList[i].SubIndex;
            }

          NHalos += CountUsed;
          TreeLen[NTrees] = CountUsed;
          NTrees++;
        }
    }
  return valid;
}

void generate_trees(int iter)
{
  NTrees = 0;
  NHalos = 0;
  MaxNTrees = 10;
  MaxNHalos = 10000;
  MaxTreeLen = 10;

  TreeLen = mymalloc_movable(&TreeLen, "TreeLen", MaxNTrees * sizeof(int));     
  TreeHalos = mymalloc_movable(&TreeHalos, "TreeHalos", MaxNHalos * sizeof(struct tree_data));
  HaloList = mymalloc_movable(&HaloList, "HaloList", MaxTreeLen * sizeof(struct halo_data));    

  for(int i = 0; i < LocNsubgroups; i++) Halo[i].UsedFlag = Halo[i].HaloFlag = 0;

  for(int i = 0; i < LocNsubgroups; i++)
    {
      if(Halo[i].SnapNum == LastSnapShotNr && Halo[i].UsedFlag == 0)
        {
          CountUsed = 0;
          walk_it(i);

          for(int k = 0; k < CountUsed; k++)
            {
              if(HaloList[k].Descendant >= 0)
                HaloList[k].Descendant = Halo[HaloList[k].Descendant].TargetIndex;
              if(HaloList[k].FirstProgenitor >= 0)
                HaloList[k].FirstProgenitor = Halo[HaloList[k].FirstProgenitor].TargetIndex;
              if(HaloList[k].NextProgenitor >= 0)
                HaloList[k].NextProgenitor = Halo[HaloList[k].NextProgenitor].TargetIndex;
              if(HaloList[k].FirstHaloInFOFgroup >= 0)
                HaloList[k].FirstHaloInFOFgroup = Halo[HaloList[k].FirstHaloInFOFgroup].TargetIndex;
              if(HaloList[k].NextHaloInFOFgroup >= 0)
                HaloList[k].NextHaloInFOFgroup = Halo[HaloList[k].NextHaloInFOFgroup].TargetIndex;
            }
          check_tree_candidate(iter);
        }
    }
}


/*=============================================================================
 * MODIFICACIONES PARA GADGET-4 HDF5 - RUTINAS DE GUARDADO
 *=============================================================================*/
#ifdef HAVE_HDF5
void hdf5_create_dataset_1d(hid_t group_id, const char* name, hid_t type, hsize_t size) {
    hid_t space = H5Screate_simple(1, &size, NULL);
    hid_t dset = H5Dcreate2(group_id, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dset);
    H5Sclose(space);
}

void hdf5_create_dataset_2d(hid_t group_id, const char* name, hid_t type, hsize_t size1, hsize_t size2) {
    hsize_t dims[2] = {size1, size2};
    hid_t space = H5Screate_simple(2, dims, NULL);
    hid_t dset = H5Dcreate2(group_id, name, type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose(dset);
    H5Sclose(space);
}

void hdf5_write_slab_1d(hid_t file_id, const char* name, hid_t type, hsize_t offset, hsize_t count, void* data) {
    hid_t dset = H5Dopen2(file_id, name, H5P_DEFAULT);
    hid_t memspace = H5Screate_simple(1, &count, NULL);
    hid_t filespace = H5Dget_space(dset);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, NULL);
    H5Dwrite(dset, type, memspace, filespace, H5P_DEFAULT, data);
    H5Dclose(dset);
    H5Sclose(memspace);
    H5Sclose(filespace);
}

void hdf5_write_slab_2d(hid_t file_id, const char* name, hid_t type, hsize_t offset, hsize_t count, hsize_t dim2, void* data) {
    hid_t dset = H5Dopen2(file_id, name, H5P_DEFAULT);
    hsize_t mem_dims[2] = {count, dim2};
    hid_t memspace = H5Screate_simple(2, mem_dims, NULL);
    hid_t filespace = H5Dget_space(dset);
    hsize_t start[2] = {offset, 0};
    hsize_t dcount[2] = {count, dim2};
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, dcount, NULL);
    H5Dwrite(dset, type, memspace, filespace, H5P_DEFAULT, data);
    H5Dclose(dset);
    H5Sclose(memspace);
    H5Sclose(filespace);
}
#endif

/* 
 * save_trees: Versión actualizada para escribir nativamente en HDF5. 
 * Organiza los datos en un esquema Structure of Arrays (SoA) dentro 
 * del grupo /TreeData para maximizar la compatibilidad con L-GalaxiesBH. 
 * Las escrituras se sincronizan y distribuyen en bloques lógicos (hyperslabs) 
 * iterando sobre las MPI tasks hacia un archivo HDF5 central.
 */
void save_trees(int iter)
{
#ifndef HAVE_HDF5
  if(ThisTask == 0) {
    printf("ERROR FATAL: El código requiere estar compilado con HDF5 (HAVE_HDF5) para esta rutina de guardado.\n");
    exit(1);
  }
#else
  char buf[2*MAXLEN_PATH];
  long long glob_NTrees = 0, glob_NHalos = 0;
  long long local_NTrees = NTrees, local_NHalos = NHalos;
  long long offset_NTrees = 0, offset_NHalos = 0;

  /* 1. Preparar directorio base (solo rank 0) */
  if(ThisTask == 0)
    {
      sprintf(buf, "%s/treedata/trees_%03d", All.OutputDir, LastSnapShotNr);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef LITE_MERGERTREE
  sprintf(buf, "%s/treedata/trees_%03d/%s_%03d.%d.hdf5", All.OutputDir, LastSnapShotNr, "tree", LastSnapShotNr, iter);
#else
  sprintf(buf, "%s/treedata/trees_%03d/%s_%03d.%d.hdf5", All.OutputDir, LastSnapShotNr, "lite_tree", LastSnapShotNr, iter);
#endif

  /* 2. Reducir estadísticas globales y calcular prefijos de offset (Exscan) */
  int maxlen = 0, len = 0;
  for(int i = 0; i < NTrees; i++)
    if(len < TreeLen[i])
      len = TreeLen[i];

  MPI_Reduce(&len, &maxlen, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_NTrees, &glob_NTrees, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_NHalos, &glob_NHalos, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Exscan(&local_NTrees, &offset_NTrees, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Exscan(&local_NHalos, &offset_NHalos, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  
  if(ThisTask == 0) {
    offset_NTrees = 0;
    offset_NHalos = 0;
  }

  /* 3. Nodo 0: Crea el esqueleto del HDF5 (archivo, atributos y datasets estáticos) */
  if(ThisTask == 0)
    {
      printf("Writing %lld trees containing %lld halos in GADGET-4 HDF5 format \n", glob_NTrees, glob_NHalos);
      printf("Longest tree has %d nodes\n", maxlen);

      hid_t file_id = H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
      /* Escribimos la cabecera (Group Header) */
      hid_t header_id = H5Gcreate(file_id, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      hid_t attr_space = H5Screate(H5S_SCALAR);
      
      hid_t attr1 = H5Acreate(header_id, "Ntrees", H5T_NATIVE_ULLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr1, H5T_NATIVE_ULLONG, &glob_NTrees);
      H5Aclose(attr1);
      
      hid_t attr2 = H5Acreate(header_id, "Nhalos", H5T_NATIVE_ULLONG, attr_space, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr2, H5T_NATIVE_ULLONG, &glob_NHalos);
      H5Aclose(attr2);
      
      H5Sclose(attr_space);
      H5Gclose(header_id);

      /* Inicialización de las bases de datos SoA vacías (Group TreeData) */
      hid_t tree_data_id = H5Gcreate(file_id, "/TreeData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      if(glob_NTrees > 0) hdf5_create_dataset_1d(tree_data_id, "TreeLen", H5T_NATIVE_INT, glob_NTrees);
      
      if(glob_NHalos > 0) {
          hdf5_create_dataset_1d(tree_data_id, "Descendant", H5T_NATIVE_INT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "FirstProgenitor", H5T_NATIVE_INT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "NextProgenitor", H5T_NATIVE_INT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "FirstHaloInFOFgroup", H5T_NATIVE_INT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "NextHaloInFOFgroup", H5T_NATIVE_INT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "Len", H5T_NATIVE_INT, glob_NHalos);
          
          hid_t id_type = (sizeof(MyIDType) == 8) ? H5T_NATIVE_ULLONG : H5T_NATIVE_UINT;
          hdf5_create_dataset_1d(tree_data_id, "SubMostBoundID", id_type, glob_NHalos);
          
          hdf5_create_dataset_1d(tree_data_id, "SnapNum", H5T_NATIVE_INT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "FileNr", H5T_NATIVE_INT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "SubIndex", H5T_NATIVE_INT, glob_NHalos);

#ifndef LITE_MERGERTREE
          hdf5_create_dataset_1d(tree_data_id, "M_Mean200", H5T_NATIVE_FLOAT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "M_Crit200", H5T_NATIVE_FLOAT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "M_TopHat200", H5T_NATIVE_FLOAT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "VelDisp", H5T_NATIVE_FLOAT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "Vmax", H5T_NATIVE_FLOAT, glob_NHalos);
          hdf5_create_dataset_1d(tree_data_id, "RadMhalf", H5T_NATIVE_FLOAT, glob_NHalos);

          hdf5_create_dataset_2d(tree_data_id, "Pos", H5T_NATIVE_FLOAT, glob_NHalos, 3);
          hdf5_create_dataset_2d(tree_data_id, "Vel", H5T_NATIVE_FLOAT, glob_NHalos, 3);
          hdf5_create_dataset_2d(tree_data_id, "Spin", H5T_NATIVE_FLOAT, glob_NHalos, 3);
#endif
      }
      H5Gclose(tree_data_id);
      H5Fclose(file_id);
    }
  MPI_Barrier(MPI_COMM_WORLD);

  /* 4. Escritura sincronizada (round-robin) para poblar los arrays */
  hsize_t off_T = offset_NTrees;
  hsize_t off_H = offset_NHalos;
  
  for(int it = 0; it < NTask; it++)
    {
      if(it == ThisTask && (NTrees > 0 || NHalos > 0))
        {
          hid_t file_id = H5Fopen(buf, H5F_ACC_RDWR, H5P_DEFAULT);

          if (NTrees > 0) {
              hdf5_write_slab_1d(file_id, "/TreeData/TreeLen", H5T_NATIVE_INT, off_T, NTrees, TreeLen);
          }

          if (NHalos > 0) {
              int *int_buf = mymalloc("int_buf", NHalos * sizeof(int));

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].Descendant;
              hdf5_write_slab_1d(file_id, "/TreeData/Descendant", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].FirstProgenitor;
              hdf5_write_slab_1d(file_id, "/TreeData/FirstProgenitor", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].NextProgenitor;
              hdf5_write_slab_1d(file_id, "/TreeData/NextProgenitor", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].FirstHaloInFOFgroup;
              hdf5_write_slab_1d(file_id, "/TreeData/FirstHaloInFOFgroup", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].NextHaloInFOFgroup;
              hdf5_write_slab_1d(file_id, "/TreeData/NextHaloInFOFgroup", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].Len;
              hdf5_write_slab_1d(file_id, "/TreeData/Len", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].SnapNum;
              hdf5_write_slab_1d(file_id, "/TreeData/SnapNum", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].FileNr;
              hdf5_write_slab_1d(file_id, "/TreeData/FileNr", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              for(int i=0; i<NHalos; i++) int_buf[i] = TreeHalos[i].SubIndex;
              hdf5_write_slab_1d(file_id, "/TreeData/SubIndex", H5T_NATIVE_INT, off_H, NHalos, int_buf);

              myfree(int_buf);

              /* Identificador de ligadura con el tipo adecuado para L-Galaxies */
              MyIDType *id_buf = mymalloc("id_buf", NHalos * sizeof(MyIDType));
              for(int i=0; i<NHalos; i++) id_buf[i] = TreeHalos[i].SubMostBoundID;
              hid_t id_type = (sizeof(MyIDType) == 8) ? H5T_NATIVE_ULLONG : H5T_NATIVE_UINT;
              hdf5_write_slab_1d(file_id, "/TreeData/SubMostBoundID", id_type, off_H, NHalos, id_buf);
              myfree(id_buf);

#ifndef LITE_MERGERTREE
              float *float_buf = mymalloc("float_buf", NHalos * sizeof(float));
              
              for(int i=0; i<NHalos; i++) float_buf[i] = TreeHalos[i].M_Mean200;
              hdf5_write_slab_1d(file_id, "/TreeData/M_Mean200", H5T_NATIVE_FLOAT, off_H, NHalos, float_buf);

              for(int i=0; i<NHalos; i++) float_buf[i] = TreeHalos[i].M_Crit200;
              hdf5_write_slab_1d(file_id, "/TreeData/M_Crit200", H5T_NATIVE_FLOAT, off_H, NHalos, float_buf);

              for(int i=0; i<NHalos; i++) float_buf[i] = TreeHalos[i].M_TopHat200;
              hdf5_write_slab_1d(file_id, "/TreeData/M_TopHat200", H5T_NATIVE_FLOAT, off_H, NHalos, float_buf);

              for(int i=0; i<NHalos; i++) float_buf[i] = TreeHalos[i].VelDisp;
              hdf5_write_slab_1d(file_id, "/TreeData/VelDisp", H5T_NATIVE_FLOAT, off_H, NHalos, float_buf);

              for(int i=0; i<NHalos; i++) float_buf[i] = TreeHalos[i].Vmax;
              hdf5_write_slab_1d(file_id, "/TreeData/Vmax", H5T_NATIVE_FLOAT, off_H, NHalos, float_buf);

              for(int i=0; i<NHalos; i++) float_buf[i] = TreeHalos[i].RadMhalf;
              hdf5_write_slab_1d(file_id, "/TreeData/RadMhalf", H5T_NATIVE_FLOAT, off_H, NHalos, float_buf);
              
              myfree(float_buf);

              /* Matrices 2D (Nx3) */
              float *vec_buf = mymalloc("vec_buf", NHalos * 3 * sizeof(float));
              
              for(int i=0; i<NHalos; i++) { vec_buf[i*3]=TreeHalos[i].Pos[0]; vec_buf[i*3+1]=TreeHalos[i].Pos[1]; vec_buf[i*3+2]=TreeHalos[i].Pos[2]; }
              hdf5_write_slab_2d(file_id, "/TreeData/Pos", H5T_NATIVE_FLOAT, off_H, NHalos, 3, vec_buf);

              for(int i=0; i<NHalos; i++) { vec_buf[i*3]=TreeHalos[i].Vel[0]; vec_buf[i*3+1]=TreeHalos[i].Vel[1]; vec_buf[i*3+2]=TreeHalos[i].Vel[2]; }
              hdf5_write_slab_2d(file_id, "/TreeData/Vel", H5T_NATIVE_FLOAT, off_H, NHalos, 3, vec_buf);

              for(int i=0; i<NHalos; i++) { vec_buf[i*3]=TreeHalos[i].Spin[0]; vec_buf[i*3+1]=TreeHalos[i].Spin[1]; vec_buf[i*3+2]=TreeHalos[i].Spin[2]; }
              hdf5_write_slab_2d(file_id, "/TreeData/Spin", H5T_NATIVE_FLOAT, off_H, NHalos, 3, vec_buf);
              
              myfree(vec_buf);
#endif
          }
          H5Fclose(file_id);
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
}

/* ========================================================================= */
/*  FIN DE LA RUTINA HDF5                                                    */
/* ========================================================================= */

void set_descendant_pointers(void)
{
  int i, j, dir;
  float minsep[3], maxsep[3], dx[3], glob_minsep[3], glob_maxsep[3];

  for(j = 0; j < 3; j++)
    {
      minsep[j] = All.BoxSize;
      maxsep[j] = -All.BoxSize;
    }

  for(i = 0; i < LocNsubgroups; i++) Halo[i].Descendant = -1;

  for(i = 0; i < LocNsubgroups; i++)
    {
      if((i % 100000) == 0 && ThisTask == 0) printf("%d|%d\n", i, LocNsubgroups);

      if(Halo[i].Desc_FileNr >= 0 && Halo[i].SnapNum < LastSnapShotNr)
        {
          int snapnum = Halo[i].SnapNum + 1;
          int filenr = Halo[i].Desc_FileNr;
          int subindex = Halo[i].Desc_SubIndex;
          int grnr = Halo[i].Desc_GrNr;

          int left = 0;
          int right = LocNsubgroups - 1;
          do
            {
              int mid = (left + right) / 2;

              if(snapnum < Halo[mid].SnapNum) dir = -1;
              else if(snapnum > Halo[mid].SnapNum) dir = +1;
              else if(filenr < Halo[mid].FileNr) dir = -1;
              else if(filenr > Halo[mid].FileNr) dir = +1;
              else if((grnr < Halo[mid].FoF_GrNr) && (grnr != -1)) dir = -1;
              else if((grnr > Halo[mid].FoF_GrNr) && (grnr != -1)) dir = +1;
              else if(subindex < Halo[mid].SubIndex) dir = -1;
              else if(subindex > Halo[mid].SubIndex) dir = +1;
              else dir = 0;

              if(dir == 0) { Halo[i].Descendant = mid; break; }
              if(left == right) break;

              if(dir < 0)
                {
                  if(right == mid) right--;
                  else right = mid;
                }
              if(dir > 0)
                {
                  if(left == mid) left++;
                  else left = mid;
                }
            }
          while(dir != 0);

          if(Halo[i].GhostFlag == 0)
            {
              if(Halo[i].Descendant == -1)      
                {
                  if(Halo[i].SnapNum != (LastSnapShotNr - 1))
                    {
                    for(int j = 0; j < LocNsubgroups; j++)
                        if ((Halo[j].SnapNum == snapnum) && (Halo[j].FileNr == filenr))
                        {
                            printf("i=%d snap=%d filenr=%d grnr=%lld subibdex=%d\n", j, Halo[j].SnapNum, Halo[j].FileNr, Halo[j].FoF_GrNr, Halo[j].SubIndex); 
                            if(Halo[j].SubIndex == subindex) printf("found it!!\n");
                        }
                     fflush(stdout);  
                      Terminate
                        ("[%d] strange: the descendant subhalo has not been found for a local subhalo (i=%d subindex=%d filenr=%d @ snap=%d lastsnap=%d, )!",
                         ThisTask, i, subindex, filenr, Halo[i].SnapNum, LastSnapShotNr);
                    }
                }
              else
                {
                  for(j = 0; j < 3; j++)
                    {
                      dx[j] = NEAREST(Halo[i].Pos[j] - Halo[Halo[i].Descendant].Pos[j]);
                      if(dx[j] > maxsep[j]) maxsep[j] = dx[j];
                      if(dx[j] < minsep[j]) minsep[j] = dx[j];
                    }
                }
            }
        }
    }

  MPI_Reduce(minsep, glob_minsep, 3, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(maxsep, glob_maxsep, 3, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    for(j = 0; j < 3; j++)
      printf("minsep[%d]=%f maxsep[%d]=%f \n", j, glob_minsep[j], j, glob_maxsep[j]);
}

void set_progenitor_pointers(void)
{
  int i, first, desc;

  for(i = 0; i < LocNsubgroups; i++) Halo[i].FirstProgenitor = Halo[i].NextProgenitor = -1;

  for(i = 0; i < LocNsubgroups; i++)
    {
      if((desc = Halo[i].Descendant) >= 0)
        {
          if((first = Halo[desc].FirstProgenitor) >= 0)
            {
              if(Halo[i].Desc_Weight >= Halo[first].Desc_Weight)
                {
                  Halo[i].NextProgenitor = first;
                  Halo[desc].FirstProgenitor = i;
                }
              else
                {
                  Halo[i].NextProgenitor = Halo[first].NextProgenitor;
                  Halo[first].NextProgenitor = i;
                }
            }
          else
            {
              Halo[desc].FirstProgenitor = i;
            }
        }
    }
}

void set_fof_pointers(void)
{
  int i, last;
  for(i = 0; i < LocNsubgroups; i++) Halo[i].FirstHaloInFOFgroup = Halo[i].NextHaloInFOFgroup = -1;

  if(LocNsubgroups)
    {
      Halo[0].FirstHaloInFOFgroup = 0;
      Halo[0].NextHaloInFOFgroup = -1;
    }

  for(i = 1, last = 0; i < LocNsubgroups; i++)
    {
      if(Halo[last].SnapNum == Halo[i].SnapNum && Halo[last].FoF_GrNr == Halo[i].FoF_GrNr)
        {
          Halo[last].NextHaloInFOFgroup = i;
          Halo[i].FirstHaloInFOFgroup = Halo[last].FirstHaloInFOFgroup;
        }
      else
        {
          Halo[last].NextHaloInFOFgroup = -1;
          Halo[i].FirstHaloInFOFgroup = i;
        }
      last = i;
    }
}
#endif