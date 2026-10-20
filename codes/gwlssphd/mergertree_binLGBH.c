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
  //if(total_number_of_cells % Ncells_per_iter != 0)
  //  Terminate("total_number_of_cells must be a multiple of the Number of cells per iteration");

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

      //if((num % 5) == 0)
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

  /* now redistribute the subhalos such that we have spatially disjoint domains */
  tree_domain_decompose(iter, PH_DECOMPOSITION);

  printf("Task %d Number of haloes: %d \n", ThisTask, LocNsubgroups);

  /* now extend the domains by ghost region that needs to be thick enough to
     enclose all local trees completely */
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


/*
 * set_desc_grnr: Sorts local subhalos, establishes temporary
 * descendant pointers via binary search, and assigns the
 * Friends-of-Friends group number of the descendant (Desc_GrNr) back
 * to each subhalo for subsequent tree operations.
 */

void set_desc_grnr()
{

  MPI_Barrier(MPI_COMM_WORLD);
  mpi_printf("Setting the Desc_GrNr key in subhalos catalogues\n");

  /* now sort all the local subhalos  */
  qsort(Halo, LocNsubgroups, sizeof(struct halo_data), comp_halo_order);

  set_descendant_pointers();

  for (int i =0; i < LocNsubgroups; i++)
      if (Halo[i].Descendant != -1)
        Halo[i].Desc_GrNr = Halo[Halo[i].Descendant].FoF_GrNr;
}


#ifdef SELECT_MAIN_BRANCHES

/*
 * process_tree: Initializes auxiliary tracking structures and resets * mass accumulation flags across all nodes in the given merger tree
 * before initiating recursive branch processing.
 */
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


/*
 * walk: Recursively traverses progenitor branches to compute and
 * accumulate the integrated historical mass (MassInStem) for each
 * halo within the tree.
 */
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


/*
 * modify_tree: Corrects main progenitor branches and FOF group
 * structures based on integrated mass history (MassInStem).
 * 1. Progenitor Loop: Re-evaluates mergers to ensure the main branch
 * follows the historically most massive lineage, preventing
 * unphysical lineage jumps caused by transient snapshot mass
 * fluctuations.
 * 2. FOF Group Loop: Ensures that the central subhalo of each
 * Friends-of-Friends group is the historically dominant one,
 * reassigning global halo properties and substructure pointers
 * accordingly for semi-analytic models.
 */

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
      if(ProcessHalo[nr].FirstHaloInFOFgroup == nr && HaloAux[nr].DoneFlag == 0)        /* let's check this FOF group */
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


/*
 * comp_halo_order_grnr: Comparison function used for sorting subhalos
 * hierarchically by snapshot number, file number, FOF group number,
 * and subhalo index.
 */

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

/*
 * comp_halo_order: Comparison function used for sorting subhalos
 * hierarchically by snapshot number, file number, and subhalo index.
 */

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


/*
 * tree_load_subhalo_catalogue: Loads Friends-of-Friends (FOF) group
 * and SUBFIND subhalo catalogues from disk for a given snapshot, maps
 * their properties and descendant pointers into dynamic memory, and
 * prepares the structures required for merger tree construction.
 */

void tree_load_subhalo_catalogue(int num, int lastsnap)       /*ok */
{

#ifndef HASHTABLE
    int load_id_list = 0;
    fof_subfind_load_groups(num, load_id_list); //not load id_list
    if (num != lastsnap)
        fof_subfind_load_descendants(num);

    if(Halo == NULL)
        Halo = mymalloc_movable(&Halo, "Halo", (LocNsubgroups + Nsubgroups) * sizeof(struct halo_data));
    else
        Halo = myrealloc_movable(Halo, (LocNsubgroups + Nsubgroups) * sizeof(struct halo_data));

    /* now fill the data into the halos */

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

  //load_subhalo_catalogue_single(num, nr);
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

      //Discard those that are not in the selection
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

      /* determine local extension */
      for(j = 0; j < 3; j++)
        {
          xmin[j] = 1.0e30;
          xmax[j] = -1.0e30;
        }
      for(i = 0; i < LocNsubgroups; i++)
        {
          for(j = 0; j < 3; j++)
            {
              if(xmin[j] > Halo[i].Pos[j])
                xmin[j] = Halo[i].Pos[j];

              if(xmax[j] < Halo[i].Pos[j])
                xmax[j] = Halo[i].Pos[j];
            }
        }

      MPI_Allreduce(xmin, glob_xmin, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(xmax, glob_xmax, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      len = 0;
      for(j = 0; j < 3; j++)
        if(glob_xmax[j] - glob_xmin[j] > len)
          len = glob_xmax[j] - glob_xmin[j];

      for(j = 0; j < 3; j++)
        DomainCenterTree[j] = 0.5 * (glob_xmin[j] + glob_xmax[j]);

      DomainRadius = 0.5 * len * sqrt(3);

      LocNsubgroups = 0;

    }
  // else {

  mark_hash_cells(iter, 0);

  for(ii = ThisTask; ii < NumFilestoLoad; ii += NTask)
    {
      printf("Loading %d out of %d file %d\n", ii, NumFilestoLoad, FilestoLoad[ii]);
      load_subhalo_catalogue_single(num, FilestoLoad[ii]);
    }

  //Discard those that are not in the selection
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
  // }
  //  printf("Holding [%d] %d\n", ThisTask, LocNsubgroups);
//
  HighMarkBytes_load = 0;
  report_memory_usage(&HighMarkBytes_load, "LOAD");

  myfree_movable(FilestoLoad);

  if(ThisTask == 0)
    free_hash_table();
#endif
}


/*
 * mark_hash_cells: Evaluates spatial hash cells across the simulation
 * box to determine which grid regions intersect with the local domain
 * and its ghost zones, mapping and broadcasting the required physical
 * files from disk to optimize I/O and data loading.
 */

#ifdef HASHTABLE
void mark_hash_cells(int iter, int type)
{
  int i, j, k, next, hashkey, ncellsload, flag;
  double x, y, z;
  double dx, dy, dz, rr;
  peanokey key;

  if(ThisTask == 0)
    {

      printf("Marking cells (HashTabSize = %d, BoxSize=%f, DomainFac=%f)\n", HashTabSize, All.BoxSize,
             DomainFac);

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
                if(rr < (DomainRadius + GhostRegionThickness + All.BoxSize / HashTabSize))
                  flag = 1;

              if(type == 1)
                if(ShouldWeLoadCell(iter, x, y, z, DomainFac))
                  flag = 1;


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
        if(FileMap[i] == 1)
          NumFilestoLoad++;

      printf("Loading %d cells on %d files\n", ncellsload, NumFilestoLoad);
    }

  MPI_Bcast(&HashTabSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NumFilestoLoad, 1, MPI_INT, 0, MPI_COMM_WORLD);
  FilestoLoad = mymalloc_movable(&FilestoLoad, "FilestoLoad", sizeof(int) * NumFilestoLoad);

  int n = 0;

  if(ThisTask == 0)
    for(i = 0; i < NFiles; i++)
      if(FileMap[i] == 1)
        {
          FilestoLoad[n] = i;
          n++;
        }


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
  /* get peano-hilbert number */
  int task;

  if(type == PH_DECOMPOSITION)
    {
      if(typePart == 1)         //this means it is a ghost particle
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
          if (task < 0)
              task = 0;
          if (task >= NTask)
              task = NTask-1;

          //task = (key - Ncells_per_iter * iter) / ((float) Ncells_per_iter / NTask);
        }

    }
  else if(type == RANDOM_DECOMPOSITION)
    {
      task = (int) floor((pos[0] - (int) pos[0]) * NTask);
//     printf("%f %f %f, task=%d\n",(pos[0]-(int)pos[0]), (pos[0]-(int)pos[0])*NTask, pos[0], task);
    }
  else
  Terminate("get_target_task mode =%d not recognized\n", type) if(task >= NTask)

  if(task >= NTask)  /* esto no está en el código original*/
    task = NTask - 1;
  assert(task >= 0 && task < NTask);

  return task;
}


/*
 * tree_domain_decompose: Redistributes subhalos across MPI tasks
 * based on spatial domain criteria, handling parallel communication
 * (Alltoall and Sendrecv) to migrate halos to their corresponding
 * target processors.
 */

void tree_domain_decompose(int iter, int type)
{
  int n, target, ngrp;
  long long i;

  for(n = 0; n < NTask; n++)
    {
      Send_count[n] = 0;
      Recv_count[n] = 0;
    }

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

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

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


/*
 * periodic_wrap: Wraps coordinates or distances into the central
 * periodic box domain [-0.5*BoxSize, 0.5*BoxSize] to account for
 * periodic boundary conditions in cosmological simulations.
 */

double periodic_wrap(double x)
{
  if(x > 0.5 * All.BoxSize)
    x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize)
    x += All.BoxSize;

  return x;
}


/*
 * get_ghost_halos: Computes local domain centers and radii, shares
 * them globally across MPI tasks via Allgather, and exchanges
 * boundary subhalos (ghost halos) with neighboring domains to ensure
 * seamless cross-boundary tree construction.
 */

void get_ghost_halos(void)
{
  /* first, let's determine the center's and radii of all domains */
  int i, j, ngrp;
  double len, xmin[3], xmax[3];

  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
      xmin[j] = 1.0e30;
      xmax[j] = -1.0e30;
    }

  for(i = 0; i < LocNsubgroups; i++)
    {
      if(Halo[i].SnapNum == LastSnapShotNr)
        for(j = 0; j < 3; j++)
          {
            if(xmin[j] > Halo[i].Pos[j])
              xmin[j] = Halo[i].Pos[j];

            if(xmax[j] < Halo[i].Pos[j])
              xmax[j] = Halo[i].Pos[j];
          }
    }

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax[j] - xmin[j] > len)
      len = xmax[j] - xmin[j];

  double domainCenter[3], domainRadius;

  for(j = 0; j < 3; j++)
    domainCenter[j] = 0.5 * (xmin[j] + xmax[j]);

  domainRadius = 0.5 * len * sqrt(3);

  if(LocNsubgroups == 0)
    DomainRadius = -1;

  //if(ThisTask == 0 || ThisTask == NTask - 1)
    printf("ThisTask=%d len=%f rad=%f DomainCenter=(%f,%f,%f) (%f-%f) (%f-%f) (%f-%f)\n", ThisTask, len, domainRadius,  domainCenter[0], domainCenter[1], domainCenter[2],
            xmin[0], xmax[0], xmin[1], xmax[1], xmin[2], xmax[2]);


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

          MPI_Sendrecv(&TreeNexport,
                       1, MPI_INT, target, TAG_N,
                       &TreeNimport, 1, MPI_INT, target, TAG_N, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          if(Halo == NULL)
            Halo = mymalloc_movable(&Halo, "Halo", (recv_off + TreeNimport) * sizeof(struct halo_data));
          else
            Halo = myrealloc_movable(Halo, (recv_off + TreeNimport) * sizeof(struct halo_data));

          MPI_Sendrecv(halobuf,
                       TreeNexport * sizeof(struct halo_data), MPI_BYTE, target, TAG_PDATA,
                       Halo + recv_off,
                       TreeNimport * sizeof(struct halo_data), MPI_BYTE, target, TAG_PDATA,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          recv_off += TreeNimport;
        }
    }

  LocNsubgroups = recv_off;

  printf("Task %d Number of haloes + Ghost region: %d , Exported: %d\n", ThisTask, LocNsubgroups,
         TreeNexport);

  myfree(halobuf);
  myfree(DomainRadiusTable);
  myfree(DomainCenterTable);
}


/*
 * walk_it: Recursively traverses and extracts connected halos
 * (following descendants, progenitors, and FOF group links) into a
 * sequential, dynamically resized flat array (HaloList) to build an
 * independent merger tree.
 */


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
    {
      if(Halo[Halo[i].Descendant].UsedFlag == 0)
        walk_it(Halo[i].Descendant);
    }

  p = Halo[i].FirstProgenitor;
  while(p >= 0)
    {
      if(Halo[p].UsedFlag == 0)
        walk_it(p);

      p = Halo[p].NextProgenitor;
    }

  p = Halo[i].FirstHaloInFOFgroup;
  if(Halo[p].HaloFlag == 0)
    {
      Halo[p].HaloFlag = 1;
      while(p >= 0)
        {
          if(Halo[p].UsedFlag == 0)
            walk_it(p);
          p = Halo[p].NextHaloInFOFgroup;
        }
    }
}


/*
 * check_tree_candidate: Validates whether the extracted tree belongs
 * to the local MPI task (contains native z=0 halos and wins the
 * deterministic ID ownership check), reallocates global memory
 * buffers if necessary, and ermanently stores the tree nodes into the
 * global TreeHalos and TreeLen arrays.
 */

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

          //if(target == ThisTask)
          if(HaloList[k].GhostFlag == 0)
            valid = 1;          /* this tree should be valid because it
                                   contains at least one local subhalo */

          if(kmin < 0)
            kmin = k;
          else
            {
              if(HaloList[k].SubMostBoundID < HaloList[kmin].SubMostBoundID)
                kmin = k;
            }
        }
    }

  if(valid)
    {
      target = get_target_task(iter, HaloList[kmin].Pos, 1, HaloList[kmin].GhostFlag);
      //if(target == ThisTask)
      if(HaloList[kmin].GhostFlag == 0)
        {
          /* so we identified a new tree for the local domain */

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
              /* here, fill in the relevant tree data that we want to keep */

              TreeHalos[NHalos + i].Descendant = HaloList[i].Descendant;
              TreeHalos[NHalos + i].FirstProgenitor = HaloList[i].FirstProgenitor;
              TreeHalos[NHalos + i].NextProgenitor = HaloList[i].NextProgenitor;
              TreeHalos[NHalos + i].FirstHaloInFOFgroup = HaloList[i].FirstHaloInFOFgroup;
              TreeHalos[NHalos + i].NextHaloInFOFgroup = HaloList[i].NextHaloInFOFgroup;

              /* properties of subhalo */

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

              // MODIFICACIÓN PARA LGBH: Evitar asignar memoria no declarada
              // for(j = 0; j < 3; j++)
              //  TreeHalos[NHalos + i].CM[j] = HaloList[i].CM[j];

#endif
              TreeHalos[NHalos + i].SubMostBoundID = HaloList[i].SubMostBoundID;

              /* original position in subfind output */
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


/*
 * generate_trees: Coordinates the traversal and root-finding of
 * merger trees starting from z=0 snapshots, invoking recursive tree
 * extraction (walk_it), reindexing pointers to match local array
 * offsets, and validating tree candidates for local MPI ownership.
 */

void generate_trees(int iter)
{

  NTrees = 0;
  NHalos = 0;

  MaxNTrees = 10;
  MaxNHalos = 10000;
  MaxTreeLen = 10;

  TreeLen = mymalloc_movable(&TreeLen, "TreeLen", MaxNTrees * sizeof(int));     /* gives the number of halos in each tree */

  TreeHalos = mymalloc_movable(&TreeHalos, "TreeHalos", MaxNHalos * sizeof(struct tree_data));

  HaloList = mymalloc_movable(&HaloList, "HaloList", MaxTreeLen * sizeof(struct halo_data));    /* temporarily holds the halos of current tree */

  for(int i = 0; i < LocNsubgroups; i++)
    Halo[i].UsedFlag = Halo[i].HaloFlag = 0;


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

          /* check whether this tree is local and we really want to keep it */

          check_tree_candidate(iter);

        }
    }


  //MPI_Barrier(MPI_COMM_WORLD);
  //exit(9);

}


/*
 * save_trees: Manages parallel I/O and synchronized disk writing of
 * the generated merger trees across MPI tasks, computing global
 * statistics (total trees, total halos, and maximum tree length) and
 * sequentially writing tree lengths and node structures into binary
 * output files.
 */

void save_trees(int iter)
{
  FILE *fd;
  char buf[2*MAXLEN_PATH];
  int glob_NTrees, glob_NHalos;


  if(ThisTask == 0)
    {
      /* sprintf(buf, "%s/treedata/trees_%03d", All.OutputDir, LastSnapShotNr); */
      sprintf(buf, "/cosmos_storage/home/acascales/gw_lss_phd/codes/trees_output/TNG300-mimic/trees_%03d", LastSnapShotNr);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef LITE_MERGERTREE
  /*sprintf(buf, "%s/treedata/trees_%03d/%s_%03d.%d", All.OutputDir, LastSnapShotNr, "tree", LastSnapShotNr, iter); */
  sprintf(buf, "/cosmos_storage/home/acascales/gw_lss_phd/codes/trees_output/TNG300-mimic/trees_%03d/%s_%03d.%d", LastSnapShotNr, "tree", LastSnapShotNr, iter);
#else
  sprintf(buf, "%s/treedata/trees_%03d/%s_%03d.%d", All.OutputDir, LastSnapShotNr, "lite_tree", LastSnapShotNr, iter);
#endif

/*
  int itree = ntrees = count_tree = 0;
  for (int i = 0; i < NTrees; i++)
  {
    int oldLen = TreeLen[i];
    int root = TreeHalos[itree].FirstHaloInFOFgroup + itree;
    if (ShouldWeLoadCell(TreeHalos[root].Pos[0], TreeHalos[root].Pos[1], TreeHalos[root].Pos[2], DomainFac))
    {
       memcpy(&TreeHalos[count_tree],&TreeHalos[itree],TreeLen[i]*sizeof(struct tree_data));
       TreeLen[ntrees] = TreeLen[i];
       ntrees++;
       count_tree += TreeLen[i];
    } else
      printf("not in the region\n");
    itree += oldLen;
  }
  NTrees = ntrees;
*/

  int maxlen = 0, len = 0;
  for(int i = 0; i < NTrees; i++)
    if(len < TreeLen[i])
      len = TreeLen[i];

  MPI_Reduce(&len, &maxlen, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&len, &maxlen, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&NTrees, &glob_NTrees, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&NHalos, &glob_NHalos, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(!(fd = fopen(buf, "w")))
        {
          printf("can't open file `%s'\n", buf);
          exit(1);
        }
      printf("Writing %d trees containing %d halos \n", glob_NTrees, glob_NHalos);

      printf("Longest tree has %d nodes\n", maxlen);

      // MODIFICACIÓN PARA LGBH: Evitamos volcar variables no leídas en el header
      // my_fwrite(&Niter, 1, sizeof(int), fd);
      // my_fwrite(&Ncells_per_iter, 1, sizeof(int), fd);
      my_fwrite(&glob_NTrees, 1, sizeof(int), fd);
      my_fwrite(&glob_NHalos, 1, sizeof(int), fd);
      fclose(fd);
    }

  for(int it = 0; it < NTask; it++)
    {
      if(it == ThisTask)
        {
          printf("task %d writing %d trees\n", ThisTask, NTrees);
          if(!(fd = fopen(buf, "a")))
            {
              printf("can't open file `%s'\n", buf);
              exit(1);
            }

          my_fwrite(TreeLen, NTrees, sizeof(int), fd);
          fclose(fd);
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }

  for(int it = 0; it < NTask; it++)
    {
      if(it == ThisTask)
        {
          if(!(fd = fopen(buf, "a")))
            {
              printf("can't open file `%s'\n", buf);
              exit(1);
            }
          int count_write = 0;
          for(int i = 0; i < NTrees; i++)
            {
              my_fwrite(TreeHalos + count_write, TreeLen[i], sizeof(struct tree_data), fd);
              count_write += TreeLen[i];
            }

          fclose(fd);
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }
}


/*
 * set_descendant_pointers: Connects each subhalo to its descendant in
 * the subsequent snapshot by performing an efficient binary search
 * (bisection) over sorted local subhalo catalogues, verifying valid
 * mappings, and computing global spatial separation statistics across
 * MPI tasks.
 */

void set_descendant_pointers(void)
{
  int i, j, dir;
  float minsep[3], maxsep[3], dx[3], glob_minsep[3], glob_maxsep[3];

  for(j = 0; j < 3; j++)
    {
      minsep[j] = All.BoxSize;
      maxsep[j] = -All.BoxSize;
    }

  for(i = 0; i < LocNsubgroups; i++)
    Halo[i].Descendant = -1;

  for(i = 0; i < LocNsubgroups; i++)
    {
      if((i % 100000) == 0 && ThisTask == 0)
        printf("%d|%d\n", i, LocNsubgroups);

      if(Halo[i].Desc_FileNr >= 0 && Halo[i].SnapNum < LastSnapShotNr)
        {
          int snapnum = Halo[i].SnapNum + 1;
          int filenr = Halo[i].Desc_FileNr;
          int subindex = Halo[i].Desc_SubIndex;
          int grnr = Halo[i].Desc_GrNr;

          /* let's try to find this entry with bisection */
          int left = 0;
          int right = LocNsubgroups - 1;

          do
            {
              int mid = (left + right) / 2;

              if(snapnum < Halo[mid].SnapNum)
                dir = -1;
              else if(snapnum > Halo[mid].SnapNum)
                dir = +1;
              else if(filenr < Halo[mid].FileNr)
                dir = -1;
              else if(filenr > Halo[mid].FileNr)
                dir = +1;
              else if((grnr < Halo[mid].FoF_GrNr) && (grnr != -1))
                dir = -1;
              else if((grnr > Halo[mid].FoF_GrNr) && (grnr != -1))
                dir = +1;
              else if(subindex < Halo[mid].SubIndex)
                dir = -1;
              else if(subindex > Halo[mid].SubIndex)
                dir = +1;
              else
                dir = 0;

              /*
                 printf("left=%d right=%d mid=%d dir=%d  (%d|%d|%d) (%d|%d|%d) \n",
                 left, right, mid, dir,
                 snapnum, filenr, subindex,
                 Halo[mid].SnapNum, Halo[mid].FileNr, Halo[mid].SubIndex);
               */

              if(dir == 0)
                {
                  Halo[i].Descendant = mid;
                  break;
                }

              if(left == right)
                break;

              if(dir < 0)
                {
                  if(right == mid)
                    right--;
                  else
                    right = mid;
                }

              if(dir > 0)
                {
                  if(left == mid)
                    left++;
                  else
                    left = mid;
                }
            }
          while(dir != 0);

          if(Halo[i].GhostFlag == 0)
            {
              if(Halo[i].Descendant == -1)      /* no descendant found: this may not happen for the genuine local halos */
                {
                  if(Halo[i].SnapNum != (LastSnapShotNr - 1))
                    {


                    for(int j = 0; j < LocNsubgroups; j++)
                        if ((Halo[j].SnapNum == snapnum) && (Halo[j].FileNr == filenr))
                        {
                            printf("i=%d snap=%d filenr=%d grnr=%lld subibdex=%d\n", j, Halo[j].SnapNum, Halo[j].FileNr, Halo[j].FoF_GrNr, Halo[j].SubIndex);
                            if(Halo[j].SubIndex == subindex)
                                printf("found it!!\n");
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
                      if(dx[j] > maxsep[j])
                        maxsep[j] = dx[j];
                      if(dx[j] < minsep[j])
                        minsep[j] = dx[j];
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


/*
 * set_progenitor_pointers: Establishes and links progenitor chains
 * for subhalos pointing toward their common descendants, sorting them
 * hierarchically by fusion weight (Desc_Weight) to correctly
 * designate main (FirstProgenitor) and secondary (NextProgenitor)
 * branches.
 */

void set_progenitor_pointers(void)
{
  int i, first, desc;

  for(i = 0; i < LocNsubgroups; i++)
    Halo[i].FirstProgenitor = Halo[i].NextProgenitor = -1;

  for(i = 0; i < LocNsubgroups; i++)
    {
      if((desc = Halo[i].Descendant) >= 0)
        {
          if((first = Halo[desc].FirstProgenitor) >= 0)
            {
              //if(Halo[i].Len >= Halo[first].Len)
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


/*
 * set_fof_pointers: Establishes horizontal linked lists for subhalos
 * belonging to the same Friends-of-Friends (FOF) group within a
 * snapshot by linking neighboring elements and assigning group head
 * references (FirstHaloInFOFgroup and NextHaloInFOFgroup).
 */

void set_fof_pointers(void)
{
  int i, last;

  for(i = 0; i < LocNsubgroups; i++)
    Halo[i].FirstHaloInFOFgroup = Halo[i].NextHaloInFOFgroup = -1;

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
