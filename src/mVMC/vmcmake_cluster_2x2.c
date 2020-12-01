/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */
#define N4X4MAXUPDATE 80
/*
 * TODO: Declare external global.h.
 */
#include <stdio.h>

/*
 * Runtime verifiers.
 */
#include <assert.h>
#define assert__(x) for ( ; !(x) ; assert(x) )
#ifdef _DEBUG
#define _VERIFIER
#endif
#ifdef _VERIFIER
#define MASSERT(x, msg) assert__(x) \
  fprintf(stderr, "error: " msg)
#else
#define MASSERT(x, msg)
#endif
#define MASSERT_INTERNAL(x) MASSERT(x, "Internal error.\n")

#define TBLISNAME( funcname, cblachar ) funcname##_##cblachar
#define MVMCNAME( funcname, vmcsuffix ) funcname##vmcsuffix

#define MAKEDEF_VMCMAKE( ctype, cblachar, funcname, vmcsuffix, rtype, SlaterElm_GLOB, InvM_GLOB, PfM_GLOB ) \
  void MVMCNAME(funcname, vmcsuffix)(MPI_Comm comm) \
{ \
  int nSweep; \
  int nInterval; \
  ctype logIpOld, logIpNew; \
\
  int eleProjCntNew[NProj]; \
  int newEleCfg[Nsite2]; \
  ctype pfMNew[NQPFull]; \
  rtype x, w; \
\
  int qpStart, qpEnd; \
  /*
   * TODO: Add support for QP splitting.
   *
  int rank, size; \
\
  MPI_Comm_size(comm, &size); \
  MPI_Comm_rank(comm, &rank); \
  SplitLoop(&qpStart, &qpEnd, \
            NQPFull, rank, size); \
  */ \
  qpStart = 0; \
  qpEnd = NQPFull; \
\
  StartTimer(30); \
  if (BurnFlag == 0) { \
    makeInitialSample(TmpEleIdx, TmpEleCfg, \
                      TmpEleNum, TmpEleProjCnt, \
                      qpStart, qpEnd, comm); \
  } else { \
    copyFromBurnSample(TmpEleIdx, TmpEleCfg, \
                       TmpEleNum, TmpEleProjCnt); \
  } \
\
  /*
   * TODO: Make it input parameters.
   */ \
  if (NExUpdatePath == 0) \
    NBlockUpdateSize = 4; \
  else \
    NBlockUpdateSize = 20; \
\
  /* 
   * Set one universal EleSpn and initialize update handler.
   */ \
  for (mi=0; mi<Ne;  mi++) EleSpn[mi] = 0; \
  for (mi=Ne;mi<Ne*2;mi++) EleSpn[mi] = 1; \
\
  /*
   * Meaning of EleCfg is not consistent across mVMC.
   *   Unify to fsz case here.
   */ \
  for (osi = 0; osi < Nsite; ++osi) \
    newEleCfg[osi] = TmpEleCfg[osi] + Ne*0; \
  for (osi = Nsite; osi < Nsite2; ++osi) \
    newEleCfg[osi] = TmpEleCfg[osi] + Ne*1; \
\
  /*
   * Declare and initialize (PfM, InvM) updators.
   */ \
  void *pfOrbital[NQPFull]; \
  void *pfUpdator[NQPFull]; \
\
  TBLISNAME(updated_tdi_v_init, cblachar) \
    (NQPFull, Nsite, Nsite2, Nsize, \
     SlaterElm_GLOB, Nsite2*Nsite2, \
     InvM_GLOB, Nsize*Nsize, \
     TmpEleIdx, EleSpn, \
     NBlockUpdateSize, \
     pfUpdator, \
     pfOrbital); \
  TBLISNAME(updated_tdi_v_get_pfa, cblachar) \
    (NQPFull, \
     PfM_GLOB, \
     pfUpdator); \
\
  logIpOld = MVMCNAME(CalculateLogIP, vmcsuffix) \
    (PfM_GLOB, qpStart, qpEnd, comm); \
\
  /*
   * TODO: Remake sample on infinity.
   */ \
  StopTimer(30); \
\
  nSweep = (BurnFlag==0) ? NVMCWarmUp+NVMCSample : NVMCSample+1; \
  nInterval = NVMCInterval; \
\
  for (int iSweep = 0; iSweep < nSweep; ++iSweep) { \
    /*
     * Repeat full-size sweeping for nInterval times.
     */ \
    for (int iInterval = 0; iInterval < nInterval; ++iInterval) \
      for (int iSite = 0; iSite < Nsite; ++iSite) { \
        /*
         * Select update cluster by index.
         */ \
        int clusterPos[4], \
            clusterCfg[8], \
            clusterSpn[4], \
            clusterNum[4], \
            clusterUpper = 0; \
        clusterByIdx_square2x2(iSite, \
                               clusterPos, \
                               clusterCfg, \
                               clusterSpn, \
                               clusterNum, \
                               newEleCfg, \
                               TmpEleNum); \
        /*
         * Identify target configurations and arrange update list.
         */ \
        int clusterFrom[N4X4MAXUPDATE][4], \
            clusterTo  [N4X4MAXUPDATE][4], \
            clusterNHop[N4X4MAXUPDATE], \
            nThisClusterUpdate; \
        clusterArrangeUpdate_4(clusterPos, \
                               clusterCfg, \
                               clusterSpn, \
                               clusterNum, \
                               clusterFrom, clusterTo, \
                               clusterNHop, \
                               &nThisClusterUpdate); \
        /*
         * Compute weight of each arranged (connected) update.
         */ \
        for (int iTrial = 0; iTrial < nThisClusterUpdate; ++iTrial) { \
          /*
           * Compute only weight of this exchange update and revert.
           * TODO: Assert: TmpEleIdx, newEleCfg, TmpEleNum, pfUpdator does not change state.
           */ \
          MVMCNAME(tryVMCUpdate, vmcsuffix) \
            (clusterNHop[iTrial], \
             clusterFrom[iTrial], \
             clusterTo  [iTrial], \
             TmpEleIdx, newEleCfg, TmpEleNum, EleSpn, \
             eleProjCntNew, TmpEleProjCnt, \
             &logIpNew, \
             pfMNew, pfUpdator, \
             qpStart, qpEnd, comm); \
          logIpLst[iCCWExchange] = logIpNew; \
        } \
        /*
         * Apply Suwa-Todo to determine the next configuration.
         */ \
        int selUpdate = TBLISNAME(suwaTodoSelector, cblachar)(nThisClusterUpdate, logIpLst); \
        /*
         * TODO: Accept the selected update (Change pfUpdator, TmpEle***, newEleCfg, logIpOld, etc.)
         *   and verify logIpOld with a newly calculated Pfaffian if _VERIFIER_VERBOSE in on.
         */ \
\
      } \
    \
\
    StartTimer(35); \
    /* 
     * Convert back and save electron configuration.
     * TODO: Assert in vmccal that TmpEle*** are consistent with each other.
     */ \
    for (int osi = 0; osi < Nsite2; osi++) \
      TmpEleCfg[osi] = newEleCfg[osi] % Ne; \
    if (iSweep >= nSweep - NVMCSample) { \
      int iSample = iSweep - (nSweep - NVMCSample); \
      saveEleConfig(iSample, logIpOld, \
                    TmpEleIdx, TmpEleCfg, \
                    TmpEleNum, TmpEleProjCnt); \
    } \
    StopTimer(35); \
  } \
}

#define MAKEDEF_DOVMCUPDATE( ctype, cblachar, funcname, vmcsuffix, rtype ) \
  void MVMCNAME(funcname, vmcsuffix)(int clusterNHop, \
                                     int *clusterFrom, \
                                     int *clusterTo, \
                                     int *eleIdx, \
                                     int *eleCfg, \
                                     int *eleNum, \
                                     int *eleSpn, \
                                     int eleProjCntNew, \
                                     int eleProjCntOld, \
                                     ctype *logIpNew, \
                                     ctype *pfMNew, \
                                     void *pfUpdator, \
                                     int qpStart, \
                                     int qpEnd, \
                                     MPI_Comm comm) \
{ \
  for (int iHop = 0; iHop < clusterNHop; ++iHop) { \
    /*
     * Update Gutzwiller-Jastrow counters.
     */ \
    StartTimer(65); \
    int mi = clusterFrom[iHop]; \
    int ri = eleIdx[mi] % Nsite; \
    int rj = clusterTo[iHop]; \
    int si = eleSpn[mi]; \
    /*
     * Use fsz for updating eleCfg here, but keep origin and destination spin the same.
     */ \
    updateEleConfig_fsz(mi, ri, rj, si, si, eleIdx, eleCfg, eleNum, eleSpn); \
    UpdateProjCnt(ri, rj, si, eleProjCntNew, eleProjCntOld, eleNum); \
    StopTimer(65); \
\
    /*
     * Perform Pfaffian update.
     */ \
    StartTimer(66); \
    TBLISNAME(updated_tdi_v_push, cblachar) \
      (NQPFull, \
       rj+si*Nsite, mi, \
       iHop == clusterNHop-1, \
       pfUpdator); \
    StopTimer(66); \
  } \
  /*
   * Get Proj*Pfa -> IP and write.
   */ \
  TBLISNAME(updated_tdi_v_get_pfa, cblachar) \
    (NQPFull, pfMNew, pfUpdator); \
\
  *logIpNew = MVMCNAME(CalculateLogIP, vmcsuffix) \
    (pfMNew, qpStart, qpEnd, comm); \
}

#define MAKEDEF_TRYVMCUPDATE( ctype, cblachar, funcname, vmcsuffix, rtype ) \
  void MVMCNAME(funcname, vmcsuffix)(int clusterNHop, \
                                     int *clusterFrom, \
                                     int *clusterTo, \
                                     int *eleIdx, \
                                     int *eleCfg, \
                                     int *eleNum, \
                                     int *eleSpn, \
                                     int eleProjCntNew, \
                                     int eleProjCntOld, \
                                     ctype *logIpNew, \
                                     ctype *pfMNew, \
                                     void *pfUpdator, \
                                     int qpStart, \
                                     int qpEnd, \
                                     MPI_Comm comm) \
{ \
  MVMCNAME(doVMCUpdate, vmcsuffix) \
    (clusterNHop, \
     clusterFrom, \
     clusterTo, \
     eleIdx, eleCfg, eleNum, eleSpn, \
     eleProjCntNew, eleProjCntOld, \
     logIpNew, \
     pfMNew, pfUpdator, \
     qpStart, qpEnd, comm); \
\
  /*
   * Directly pop-out updates for this trial.
   * TODO: Shortest update-route planning.
   */ \
  for (int iHop = 0; iHop < clusterNHop; ++iHop) \
    TBLISNAME(updated_tdi_v_pop, cblachar)(NQPFull, 0, pfUpdator); \
}

/*
 * Definition expansion occurs at the end.
 */

void clusterMake_n(int size,
                   int *clusterPos,
                   int *clusterCfg,
                   int *clusterSpn,
                   int *clusterNum,
                   int *eleCfg,
                   int *eleNum)
{
  MASSERT_INTERNAL(size > 0 && size <= Nsite);

  for (int i = 0; i < size; ++i) {
    // Configuration.
    int mri_up = eleCfg[clusterPos[i]+Nsite*0];
    int mri_dn = eleCfg[clusterPos[i]+Nsite*1];
    clusterCfg[i+size*0] = mri_up;
    clusterCfg[i+size*1] = mri_dn;
    MASSERT(mri_up < Nsize, "Internal error: eleCfg invalid.\n");
    MASSERT(mri_dn < Nsize, "Internal error: eleCfg invalid.\n");

    // NOTE: eleNum is not used. eleCfg instead.
    clusterNum[i] = 0;
    clusterSpn[i] = 0;
    if (mri_up > 0)
    { clusterNum[i] += 1;
      clusterSpn[i] += 1; }
    if (mri_dn > 0)
    { clusterNum[i] += 1;
      clusterSpn[i] -= 1; }
  }
}

/*
 * TODO: Include global.h.
 */
extern int NLatticeW;
extern int NLatticeL;

void clusterByIdx_square2x2(int iSite, 
                            int clusterPos[4],
                            int clusterCfg[8],
                            int clusterSpn[4],
                            int clusterNum[4],
                            int *eleCfg,
                            int *eleNum)
{
  int xSite = iSite % NLatticeW;
  int ySite = iSite / NLatticeW;

  /*
   * Index-increasing order.
   */
  clusterPos[0] = iSite;
  clusterPos[1] = (xSite+1)%NLatticeW + ySite*NLatticeW;
  clusterPos[2] = xSite + ((ySite+1)%NLatticeL)*NLatticeW;
  clusterPos[3] = (xSite+1)%NLatticeW + ((ySite+1)%NLatticeL)*NLatticeW;

  clusterMake_n(4, clusterPos, clusterCfg, clusterSpn, clusterNum);
}

void clusterArrangeUpdate_4(int clusterPos[4],
                            int clusterCfg[8],
                            int clusterSpn[4],
                            int clusterNum[4],
                            int clusterFrom[][4], int clusterTo[][4],
                            int *clusterNHop,
                            int *nThisClusterUpdate)
{
  for (int i = 0; i < 4; ++i)
    MASSERT(clusterNum[i] == 1, "Not a spin configuration. Only spin supports cluster update at the moment.");

  /** Number of upper spin.
   */
  int nUpper = 0;
  /** Any valid index with up/dn spin.
   */
  int iUp = -1, iDn = -1;
  for (int i = 0; i < 4; ++1) {
    if (clusterSpn[i] == 0) {
      nUpper += 1;
      iUp = i;
    } else
      iDn = i;
  }
  MASSERT_INTERNAL(iUp > 0);
  MASSERT_INTERNAL(iDn > 0);
  MASSERT_INTERNAL(iUp != iDn);

  /*
   * Current config itself is of course part of Suwa-Todo.
   */
  clusterNHop[0] = 0;

  switch (nUpper) {
  /*
   * All-up or all-down.
   * These clusters do not update.
   */
  case 0:
  case 4:
    *nThisClusterUpdate = 1;
    break;

  /*
   * 1 up + 3 down.
   * Index-increasing update around current up config.
   */
  case 1:
    *nThisClusterUpdate = 4;
    MASSERT_INTERNAL(clusterSpn[iUp] = 0);

    int mi = clusterCfg[iUp+0*4];
    int ri = clusterPos[iUp];
    for (int dj = 1; dj < 4; ++dj) {
      int jDn = (iUp + dj) % 4;
      int mj = clusterCfg[jDn+1*4];
      int rj = clusterPos[jDn];
      MASSERT_INTERNAL(iUp != jDn);
      MASSERT_INTERNAL(mi != mj);
      MASSERT_INTERNAL(ri != rj);
      MASSERT_INTERNAL(clusterSpn[iUp] != clusterSpn[jDn]);
      MASSERT_INTERNAL(ri < Nsite && rj < Nsite);

      /*
       * Queue update.
       */
      clusterNHop[dj] = 2;
      clusterFrom[dj][0] = mi;
      clusterTo  [dj][0] = rj;
      clusterFrom[dj][1] = mj;
      clusterTo  [dj][1] = ri;
    }
    break;
  /*
   * TODO: This can be effectively merged with case 1.
   */
  case 3:
    *nThisClusterUpdate = 4;
    MASSERT_INTERNAL(clusterSpn[iDn] = 0);

    int mi = clusterCfg[iDn+1*4];
    int ri = clusterPos[iDn];
    for (int dj = 1; dj < 4; ++dj) {
      int jUp = (iDn + dj) % 4;
      int mj = clusterCfg[jUp+0*4];
      int rj = clusterPos[jUp];
      MASSERT_INTERNAL(iDn != jUp);
      MASSERT_INTERNAL(mi != mj);
      MASSERT_INTERNAL(ri != rj);
      MASSERT_INTERNAL(clusterSpn[iDn] != clusterSpn[jUp]);
      MASSERT_INTERNAL(ri < Nsite && rj < Nsite);

      /*
       * Keep a gauge here: Odd updates are for spin up; Even for spin down.
       *   Error "Update gauge broken." instead of "Internal error" when violated.
       *   This might be lift when giving support to itinerant electrons.
       */
      clusterNHop[dj] = 2;
      clusterFrom[dj][0] = mj;
      clusterTo  [dj][0] = ri;
      clusterFrom[dj][1] = mi;
      clusterTo  [dj][1] = rj;
    }
    break;

  /*
   * 6 Connected configurations in this case.
   */
  case 2:
    *nThisClusterUpdate = 6;
    /*
     * Identify all 4 spins.
     */
    int jUp = -1, 
        jDn = -1;
    for (int j = 0; j < 4; ++j) {
      if (clusterSpn[j] == 0 && j != iUp) jUp = j;
      if (clusterSpn[j] == 1 && j != iDn) jDn = j;
    }
    MASSERT_INTERNAL(clusterSpn[iUp] == 0 && clusterSpn[jUp] == 0 && iUp != jUp);
    MASSERT_INTERNAL(clusterSpn[iDn] == 1 && clusterSpn[jDn] == 1 && iDn != jDn);

    /*
     * Exchange one side * 4.
     */
    int iiUp[2] = { jUp, iUp };
    int jjDn[2] = { jDn, iDn };
    for (int ii = 0; ii < 2; ++ii) {
      int mi = clusterCfg[iiUp[ii]+0*4];
      int ri = clusterPos[iiUp[ii]];
      for (int jj = 0; jj < 2; ++jj) {
        int mj = clusterCfg[jjDn[jj]+1*4];
        int rj = clusterPos[jjDn[jj]];
        MASSERT_INTERNAL(mi != mj);
        MASSERT_INTERNAL(ri != rj);

        clusterNHop[1+ii*2+jj] = 2;
        clusterFrom[1+ii*2+jj][0] = mi;
        clusterTo  [1+ii*2+jj][0] = rj;
        clusterFrom[1+ii*2+jj][1] = mj;
        clusterTo  [1+ii*2+jj][1] = ri;
      }
    }
    /*
     * Exchange both sides.
     * TODO: Do this by successive updating after previous step.
     */
    int mUp1 = clusterCfg[iUp+0*4];
    int mUp2 = clusterCfg[jUp+0*4];
    int mDn1 = clusterCfg[iDn+1*4];
    int mDn2 = clusterCfg[jDn+1*4];
    int rUp1 = clusterPos[iUp];
    int rUp2 = clusterPos[jUp];
    int rDn1 = clusterPos[iDn];
    int rDn2 = clusterPos[jDn];
    MASSERT_INTERNAL(mUp1 != mUp2 && mDn1 != mDn2 && mUp1 != mDn2);
    MASSERT_INTERNAL(mUp1 != mDn1 && mUp2 != mDn2 && mUp2 != mDn1);
    MASSERT_INTERNAL(rUp1 != rUp2 && rDn1 != rDn2 && rUp1 != rDn2);
    MASSERT_INTERNAL(rUp1 != rDn1 && rUp2 != rDn2 && rUp2 != rDn1);
    clusterNHop[1+4+0] = 4;
    clusterFrom[1+4+0][0] = mUp1;
    clusterTo  [1+4+0][0] = rDn1;
    clusterFrom[1+4+0][1] = mDn1;
    clusterTo  [1+4+0][1] = rUp1;
    clusterFrom[1+4+0][2] = mUp2;
    clusterTo  [1+4+0][2] = rDn2;
    clusterFrom[1+4+0][3] = mDn2;
    clusterTo  [1+4+0][3] = rUp2;
    break;

  default:
    MASSERT(false, "Forbidden branch.");
    break;
  }
}

/*
 * TODO: Maybe it's somehow good to add single-precision support?
 */
MAKEDEF_DOVMCUPDATE( double, d, doVMCUpdate, _real, _real )
MAKEDEF_DOVMCUPDATE( double complex, z, doVMCUpdate, _fcmp, _fcmp )

MAKEDEF_TRYVMCUPDATE( double, d, tryVMCUpdate, _real, _real )
MAKEDEF_TRYVMCUPDATE( double complex, z, tryVMCUpdate, _fcmp, _fcmp )

MAKEDEF_VMCMAKE( double, d, VMCMakeCluster2x2, _real, double, SlaterElm_real, InvM_real, PfM_real )
MAKEDEF_VMCMAKE( double complex, z, VMCMakeCluster2x2, _fcmp, double, SlaterElm, InvM, PfM )

/*
 * Undefine macros.
 */
#undef TBLISNAME
#undef MVMCNAME
#undef MAKEDEF_DOVMCUPDATE
#undef MAKEDEF_TRYVMCUPDATE
#undef MAKEDEF_VMCMAKE
#undef N4X4MAXUPDATE
#undef MASSERT_INTERNAL
#undef MASSERT
#undef assert__

