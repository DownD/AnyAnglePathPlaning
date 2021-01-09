/*
This code was developed by Tansel Uras (turas@usc.edu) at USC.
The code is hosted at 'http://idm-lab.org/anyangle'.
If you use this code in your research, please  cite our SoCS paper:

T. Uras and S. Koenig,  2015. An Empirical Comparison of Any-Angle Path-Planning Algorithms. In: Proceedings of the 8th Annual Symposium on Combinatorial
Search. Code available at: http://idm-lab.org/anyangle

Bibtex:
@inproceedings{uras:15,
  author = "T. Uras and S. Koenig",
  title = "An Empirical Comparison of Any-Angle Path-Planning Algorithms",
  booktitle = {Proceedings of the 8th Annual Symposium on Combinatorial Search},
  year = "2015",
  note = "Code available at: http://idm-lab.org/anyangle",
}
*/

#include <stdio.h>
#include <stdint.h>
#include <numeric>
#include <algorithm>
#include "ScenarioLoader.h"
#include "Timer.h"

#include "ThetaStar.h"
#include "BlockAStar.h"
#include "FieldAStar.h"
#include "SubgoalAA.h"
#include "ANYA.h"

// These should be defined in the makefile

#define EXPERIMENT_A_EUC
#define EXPERIMENT_A_OCT
#define EXPERIMENT_T
#define EXPERIMENT_L
#define EXPERIMENT_F
#define EXPERIMENT_B
#define EXPERIMENT_SUB_1_A
#define EXPERIMENT_SUB_1_T
#define EXPERIMENT_SUB_2_A
#define EXPERIMENT_SUB_2_T
#define EXPERIMENT_SUB_10000_A
#define EXPERIMENT_SUB_10000_T
#define EXPERIMENT_ANYA


void LoadMap(const char *fname, std::vector<bool> &map, int &w, int &h);

#ifndef COMPILE_STATIC_LIB

int main(int argc, char **argv)
{
#ifndef ANY_ANGLE_STATISTICS
	std::cout<<"ANY_ANGLE_STATISTICS undefined. No point in running benchmarks!"<<std::endl;
	return 1;
#endif

	// Process the arguments
	char mapname[255];
	char scenname[255];
	std::vector<bool> mapData;
	int width, height;
	
	sprintf(mapname, "%s.map", argv[1]);
	sprintf(scenname, "%s.map.scen", argv[1]);
	LoadMap(mapname, mapData, width, height);	// Read the map
	ScenarioLoader scen(scenname);				// Read the scenario file
	
	Timer t;
	t.StartTimer();
	double totalTime = 0.0;
	
	// Prepare the algorithms
	std::vector<AnyAngleAlgorithm*> algorithms;
	
#ifdef EXPERIMENT_A_EUC
	// A* with Euclidean distance heuristic
	algorithms.push_back(new ThetaStar(mapData, width, height, A_STAR_EUC));
#endif
	
#ifdef EXPERIMENT_A_OCT
	// A* with Octile distance heuristic
	algorithms.push_back(new ThetaStar(mapData, width, height, A_STAR_OCT));
#endif
	
#ifdef EXPERIMENT_T
	// Theta*
	algorithms.push_back(new ThetaStar(mapData, width, height, THETA_STAR));
#endif
	
#ifdef EXPERIMENT_L
	// Lazy Theta*
	algorithms.push_back(new ThetaStar(mapData, width, height, LAZY_THETA_STAR));
#endif
	
#ifdef EXPERIMENT_F
	// Field A*
	algorithms.push_back(new FieldAStar(mapData, width, height));
#endif
	
#ifdef EXPERIMENT_B
	// Block A*
	algorithms.push_back(new BlockAStar(mapData, width, height));
#endif
	
#ifdef EXPERIMENT_SUB_1_A
	// Simple subgoal graph with A*
	algorithms.push_back(new SubgoalAA(mapData, width, height, 1, SUBGOAL_A_STAR));
#endif
	
#ifdef EXPERIMENT_SUB_1_T
	// Simple subgoal graph with Theta*
	algorithms.push_back(new SubgoalAA(mapData, width, height, 1, SUBGOAL_THETA_STAR));
#endif
	
#ifdef EXPERIMENT_SUB_2_A
	// Two-level subgoal graph with A*
	algorithms.push_back(new SubgoalAA(mapData, width, height, 2, SUBGOAL_A_STAR));
#endif
	
#ifdef EXPERIMENT_SUB_2_T
	// Two-level subgoal graph with Theta*
	algorithms.push_back(new SubgoalAA(mapData, width, height, 2, SUBGOAL_THETA_STAR));
#endif
	
#ifdef EXPERIMENT_SUB_10000_A
	// N-Level subgoal graph with A*
	algorithms.push_back(new SubgoalAA(mapData, width, height, 10000, SUBGOAL_A_STAR));
#endif
		
#ifdef EXPERIMENT_SUB_10000_T
	// N-Level subgoal graph with Theta*
	algorithms.push_back(new SubgoalAA(mapData, width, height, 10000, SUBGOAL_THETA_STAR));
#endif
	
#ifdef EXPERIMENT_ANYA
	// ANYA
	algorithms.push_back(new ANYA(mapData, width, height));
#endif

#ifdef ANY_ANGLE_STATISTICS
	for (int i = 0; i < algorithms.size(); i++)
		algorithms[i]-> SetStatisticsFiles(argv[1]);
#endif

	// Starting the experiment
    //std::cout<<"MAP: "<<mapname<<"\tAlgorithm(s):";
    //for (int i = 0; i < algorithms.size(); i++)
    //	std::cout<<" "<<algorithms[i]->GetName();
    //std::cout<<std::endl;	
	
	// Solve the instances
	for (int x = 0; x < scen.GetNumExperiments(); x++)
    {
			xyLoc s, g;
			s.x = scen.GetNthExperiment(x).GetStartX();
			s.y = scen.GetNthExperiment(x).GetStartY();
			g.x = scen.GetNthExperiment(x).GetGoalX();
			g.y = scen.GetNthExperiment(x).GetGoalY();

			for (int i = 0; i < algorithms.size(); i++)
				algorithms[i]-> FindPath(s,g);
    }
	
	totalTime += t.EndTimer();
    
    // Print statistics
#ifdef ANY_ANGLE_STATISTICS
	for (int i = 0; i < algorithms.size(); i++)
		algorithms[i]->PrintStatistics(); 
#endif
    
    // Mark the end of the experiment
	std::cout<<"Total time spent: "<<totalTime*1000<<"ms."<<std::endl<<std::endl;
	
	// Clean up
	for (int i = 0; i < algorithms.size(); i++)
		delete algorithms[i];

	return 0;
}
#endif
void LoadMap(const char *fname, std::vector<bool> &map, int &width, int &height)
{
	FILE *f;
	f = fopen(fname, "r");
	if (f)
    {
		fscanf(f, "type octile\nheight %d\nwidth %d\nmap\n", &height, &width);
		map.resize(height*width);
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				char c;
				do {
					fscanf(f, "%c", &c);
				} while (isspace(c));
				map[y*width+x] = (c == '.' || c == 'G' || c == 'S');
				//printf("%c", c);
			}
			//printf("\n");
		}
		fclose(f);
    }
}
