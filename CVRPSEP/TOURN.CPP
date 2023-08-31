#include <stdio.h>
#include <stdlib.h>
#include "memmod.h"
#include "basegrph.h"
#include "cnstrmgr.h"

/* Tournament inequalities */

void TOURN_Main(int NoOfCusts,
                int *ReadyTime,
                int *DueTime,
                int **TravelTime,
                ReachPtr RPtr, /* SupportPtr */
                double **XMatrix,
                CnstrMgrPointer CutsCMP)
{
  int i,j,k;
  int L; /* L = Level */
  int Root,Tmp;
  int LoopNr;
  char *OnPath;
  int *LIdx, *LN, *ArrTime, *OtherArrTime;
  double *LXSum;

  OnPath       = MemGetCV(NoOfCusts+1);
  LIdx         = MemGetIV(NoOfCusts+1); /* LevelIdx */
  LN           = MemGetIV(NoOfCusts+1); /* LevelNode */
  ArrTime      = MemGetIV(NoOfCusts+1); /* Arrival time */
  OtherArrTime = MemGetIV(NoOfCusts+1); /* Beginning at v_2 */
  LXSum        = MemGetDV(NoOfCusts+1); /* LevelXSum */

  for (i=0; i<=NoOfCusts; i++) OnPath[i] = 0;

  for (Root=1; Root<=NoOfCusts; Root++)
  {
    /*
    printf("Root = %d\n",Root);
    fflush(stdout);
    */

    LoopNr = 0;

    OnPath[Root-1] = 0; //Clear entrance from previous path
    OnPath[Root] = 1; //Add 1 for new root customer on path
    LN[1] = Root; //Add the first customer on the path at level 1
    LN[0] = 0; //Set level 0=0 (Depot)
    ArrTime[1] = ReadyTime[Root]; //For tournament constraints we need to set arrival time = whatever it's time window opens
    LXSum[1] = 0.0; //Our slack sum is 0 in the beginning

    L = 2; //We then set level = 2
    LIdx[L] = 1; //And start by setting the 2nd level index = 1 (We start our search from the first customer
    LN[L] = 0; //And add 0 to the current node added at level 2 (We have not added any node)

    do
    {
      /*
      printf("L = %d\n",L);
      fflush(stdout);
      */
      
      LoopNr++;
      /* if (LoopNr > 30) { printf("Stop on loopNr\n"); exit(0); } */

      /* Find the next node at this level = <Level> */
      while (LIdx[L] <= RPtr->LP[LN[L-1]].CFN)
      {
        Tmp = RPtr->LP[LN[L-1]].FAL[LIdx[L]];

        /* printf("Tmp:= %d\n",Tmp); */

        if ((Tmp > 0) && (OnPath[Tmp] == 0))
        {
          /* Extend path */
          LN[L] = Tmp;
          ArrTime[L] = ArrTime[L-1] + TravelTime[LN[L-1]][LN[L]];

          /*
          printf("ArrTime(L=%d)=%d, ReadyTime=%d\n",
                  L,ArrTime[L],ReadyTime[LN[L]]);
          */

          /* Try next at this level if ArrTime <= ReadyTime */
          if (ArrTime[L] > ReadyTime[LN[L]])
          {
            /* Update XSum */
            LXSum[L] = LXSum[L-1];
            /* printf("LXSum(%d):= %lf\n",L,LXSum[L]); */

            i = LN[L];
            for (j=1; j<=RPtr->LP[i].CBN; j++)
            {
              k = RPtr->LP[i].BAL[j];
              if (OnPath[k])
              {
                LXSum[L] += XMatrix[k][i];
                /*
                printf("LXSum(%d) += X(%d,%d) = %lf => LXSum(%d)=%lf\n",
                L,k,i,XMatrix[k][i],L,LXSum[L]);
                */
              }
            }

            if (LXSum[L] > ((1.0 * L) - 1.99))
            { /* Violation possible */
              if (ArrTime[L] > DueTime[LN[L]])
              {
                /* printf("Violated path found\n"); */
                /* Infeasible path: Violation found */
                /* Check if the path starting at level 2 is infeasible */

                OtherArrTime[2] = ReadyTime[LN[2]];
                /* printf("OtherArrTime(%d) = %d\n",2,OtherArrTime[2]); */
                for (i=3; i<=L; i++)
                {
                  OtherArrTime[i] = OtherArrTime[i-1] +
                                    TravelTime[LN[i-1]][LN[i]];
                  if (OtherArrTime[i] < ReadyTime[LN[i]])
                  OtherArrTime[i] = ReadyTime[LN[i]];
                  /* printf("OtherArrTime(%d) = %d\n",i,OtherArrTime[i]); */
                }

                /* printf("DueTime(%d)=%d\n",LN[L],DueTime[LN[L]]); */

                /* Otherwise generate the cut */
                if (OtherArrTime[L] <= DueTime[LN[L]])
                {
                  /* printf("Not dominated\n"); */
                  /* Add cut */
                  CMGR_AddCnstr(CutsCMP,
                                CMGR_CT_TOURN,0,
                                L,
                                LN,
                                (1.0 * L) - 2.0);
                }

                /*
                else
                printf("Dominated\n");
                */

                /* Continue with next node at this level */
              } /* ArrTime > DueTime */
              else
              {
                /* Violation still possible, but the path is feasible */
                /* Continue by branching from the new node */
                OnPath[LN[L]] = 1;
                /* printf("OnPath(%d) := 1\n",LN[L]); */

                L++; /* Branch to next level */
                LIdx[L] = 0; /* => begin with #1 in next pass */
                LN[L] = 0;
              }
            } /* if Violation possible */
          } /* if ArrTime > ReadyTime */
        } /* if Tmp>0 and not on path */

        (LIdx[L])++;
      } /* while LIdx[L] <= CFN[L-1] */

      /* Backtrack */

      L--;
      OnPath[LN[L]] = 0;
      /* printf("OnPath(%d) := 0\n",LN[L]); */

      (LIdx[L])++;
    } while (L > 1);
  } /* Root */

  MemFree(OnPath);
  MemFree(LIdx);
  MemFree(LN);
  MemFree(ArrTime);
  MemFree(OtherArrTime);
  MemFree(LXSum);
}

