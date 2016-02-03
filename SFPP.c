#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>



# define NUM 10000		/* Number of simulations */
# define MAX 5000    		/* Maximal depth */
# define DUB 10000 		/* Set to double the depth */



/* The algorithms maintains the situation at the boundary of the explored region. For every integer x PtimeL(mod(x)) and PtimeR(mod(x)) are the passage times to (x-t,t) and (x+t,t) respectively, for largest value of t (node at the boundary). NodeL(mod(x)) and NodeR(mod(x)) keep track whether the above nodes belong to the tree rooted at zero (value 0) or not (value 1).  */
int NodeR[DUB], NodeL[DUB], DepthL[DUB], DepthR[DUB];
double PtimeR[DUB], PtimeL[DUB];

/* These keep track of some statistics of the origin tree. */
int Heights[NUM];
int LeftDs[NUM];
int RightDs[NUM];
int Widhts[NUM];



/* Used for encoding positive and negative integer to positive array indices. */
int mod(int x);




int main()
{
  int level, left_edge, right_edge, width, max_width, left_consider, right_consider, low_left, low_right, running_node, left_max_disp, right_max_disp, running_left, running_right, cand_left, cand_right, change_left;
  int i, j;
  double left_cand, right_cand, running_value;
  FILE *fp;

  srand((unsigned)time(NULL));

  
  for (i = 0; i<NUM; i++)
    {
      level = 0;
      NodeR[0] = 0;
      NodeL[0] = 0;
      DepthR[0] = 0;
      DepthL[0] = 0;
      PtimeR[0] = 0.0;
      PtimeL[0] = 0.0;
      left_edge = 0;
      right_edge = 0;
      width = 1;
      
      left_max_disp = 0;
      right_max_disp = 0;
      max_width = 1;
      
      
      while ((width > 0) && (level < MAX - 1)) /* Running a single simulation */
	{
	  level ++;
	  NodeR[mod(level)] = 1;
	  NodeL[mod(level)] = 1;
	  DepthR[mod(level)] = 0;
	  DepthL[mod(level)] = 0;
	  PtimeR[mod(level)] = 0.0;
	  PtimeL[mod(level)] = 0.0;
	  NodeR[mod((-1)*level)] = 1;
	  NodeL[mod((-1)*level)] = 1;
	  DepthR[mod((-1)*level)] = 0;
	  DepthL[mod((-1)*level)] = 0;
	  PtimeR[mod((-1)*level)] = 0.0;
	  PtimeL[mod((-1)*level)] = 0.0;


	  
	  /* Computing the passage time and the root node for the rightmost candidate in the new layer which might involve doing the same computation for the whole right external boundary bottom-to-top. */
	  right_consider = left_edge;
	  left_consider = right_consider - level;
	  low_left = DepthR[mod(left_consider)];
	  
	  for (j = low_left + 1; j <= level; j++)    
	    {
	      running_right = left_consider + j;
	      right_cand = PtimeL[mod(running_right)] - log ((float)rand()/(float)RAND_MAX);
	      left_cand = PtimeR[mod(left_consider)] - log ((float)rand()/(float)RAND_MAX);
	      if (left_cand < right_cand)
		{
		  running_value = left_cand;
		  running_node = NodeR[mod(left_consider)];
		}
	      else
		{
		  running_value = right_cand;
		  running_node = NodeL[mod(running_right)];
		}
	      DepthL[mod(running_right)] ++;
	      PtimeL[mod(running_right)] = running_value;
	      NodeL[mod(running_right)] = running_node;
	      DepthR[mod(left_consider)] ++;
	      PtimeR[mod(left_consider)] = running_value;
	      NodeR[mod(left_consider)] = running_node;
	    }
	  change_left = running_node;
		
	      
	      



	  
	  /* Computing the passage times and root nodes for the candidates in the new layer (expect the leftmost and the rightmost) */	  
	  for (j = left_edge + 1; j <= right_edge; j++) 
	    {
	      left_consider ++;
	      right_consider ++;
	      left_cand =  PtimeR[mod(left_consider)] - log ((float)rand()/(float)RAND_MAX);
	      right_cand = PtimeL[mod(right_consider)] - log ((float)rand()/(float)RAND_MAX);
	      DepthR[mod(left_consider)] ++;
	      DepthL[mod(right_consider)] ++;	      
	      if (left_cand < right_cand)
		{
		  NodeL[mod(right_consider)] = NodeR[mod(left_consider)];
		  PtimeL[mod(right_consider)] = left_cand;
		  PtimeR[mod(left_consider)] = left_cand;
		}
	      else
		{
		  NodeR[mod(left_consider)] = NodeL[mod(right_consider)];
		  PtimeL[mod(right_consider)] = right_cand;
		  PtimeR[mod(left_consider)] = right_cand;
		}
	    }


	  /* Computing the passage time and the root node for the rightmost candidate in the new layer which might involve doing the same computation for the whole right external boundary bottom-to-top. */
	  left_consider ++;
	  right_consider ++;
	  low_right = DepthL[mod(right_consider)];

	  for (j = low_right + 1; j <= level; j++)   
	    {
	      running_left = right_consider - j;
	      right_cand = PtimeL[mod(right_consider)] - log ((float)rand()/(float)RAND_MAX);
	      left_cand = PtimeR[mod(running_left)] - log ((float)rand()/(float)RAND_MAX);
	      if (left_cand < right_cand)
		{
		  running_value = left_cand;
		  running_node = NodeR[mod(running_left)];
		}
	      else
		{
		  running_value = right_cand;
		  running_node = NodeL[mod(right_consider)];
		}
	      DepthL[mod(right_consider)] ++;
	      PtimeL[mod(right_consider)] = running_value;
	      NodeL[mod(right_consider)] = running_node;
	      DepthR[mod(running_left)] ++;
	      PtimeR[mod(running_left)] = running_value;
	      NodeR[mod(running_left)] = running_node;
	    }
	  
	  if (running_node == 0) right_edge ++;
	  if (change_left == 1) left_edge ++;



	  /* Updating the statistics: width captures the maximal number of nodes in a level; left_max_disp and right_max_disp capture min_i(l(i)) and max_i(r(i)) where l(i) and r(i) are the minimal and the maximal value of the projection of i-th level nodes */
	  width = right_edge - left_edge + 1;
	  if (width > 0)
	    {
	      if (width > max_width) max_width = width;
	      cand_left = 2*left_edge - level;
	      cand_right = 2*right_edge - level;
	      if (cand_left < left_max_disp)
		{
		  left_max_disp = cand_left;
		}
	      if (cand_right > right_max_disp)
		{
		  right_max_disp = cand_right;
		}
	    }
	}

      
      Heights[i] = level - 1;
      LeftDs[i] = left_max_disp;
      RightDs[i] = right_max_disp;
      Widhts[i] = max_width;

    }


  fp = fopen("./simulation_results", "w");
  for (i=0; i< NUM; i++)
    {
      fprintf(fp, "%d %d %d %d \n", Heights[i], LeftDs[i], RightDs[i], Widhts[i]);
    }
  
  fclose(fp);


      
  return 0;
}





int mod(int x)
{
  if (x == 0)
    {
      return 0;
    }
  if (x > 0)
    {
      return 2*x-1;
    }
  if (x < 0)
    {
      return (-2)*x;
    }
}

