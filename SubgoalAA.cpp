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

#include "SubgoalAA.h"

#ifdef ANY_ANGLE_RUNNING_IN_HOG
SubgoalAA::SubgoalAA(MapEnvironment *env, int level, int search_method)
	:AnyAngleAlgorithm::AnyAngleAlgorithm(env)
{
	SubgoalAAInitialize(level, search_method);
}
#endif

SubgoalAA::SubgoalAA(std::vector<bool> &bits, int _width, int _height, int level, int search_method)
	:AnyAngleAlgorithm::AnyAngleAlgorithm(bits, _width, _height)
{
	SubgoalAAInitialize(level, search_method);
}
void SubgoalAA::SubgoalAAInitialize(int level, int search_method)
{
    max_graph_level_ = level;
    search_method_ = search_method;
	CreateGraph(level);

	// Set up statistics
#ifdef ANY_ANGLE_STATISTICS
	SetupStatistics();
#endif
}
SubgoalAA::~SubgoalAA()
{}

// CONSTRUCTION

void SubgoalAA::CreateGraph(int level)
{
	corners_.resize(corner_locations_.size());
	SetDirections();

#ifdef ANY_ANGLE_STATISTICS
	Timer t;
	t.StartTimer();
#endif

	IdentifySubgoals();
	ComputeClearances();
	LinkSubgoals();

#ifdef ANY_ANGLE_STATISTICS
	t.EndTimer();
	ssg_construction_time_ = t.GetElapsedTime();
	t.StartTimer();
#endif

	graph_level_ = 1;

	// Reset the generated values
	ResetSearch();
	open_list_.reserve(10000);

	while(graph_level_ < max_graph_level_ && AddLevel());

#ifdef ANY_ANGLE_STATISTICS
	t.EndTimer();
	partitioning_time_ = t.GetElapsedTime();
#endif

	ResetSearch();

#ifdef ANY_ANGLE_VERBOSE
	printf("Graph level: %u\nGlobal subgoals: %u\n", graph_level_, num_global_subgoals_);
#endif

	FinalizeGraph();
}

const std::string SubgoalAA::GetName() const
{
    std::stringstream s;
    s << "S_";

    s<<max_graph_level_;

    switch (search_method_) {

    case SUBGOAL_A_STAR:
        s << "_A";
        break;

    case SUBGOAL_THETA_STAR:
        s << "_T";
        break;

    default:
        break;
    }

    return s.str();
}

// Functions to create the simple subgoal graph
void SubgoalAA::SetDirections()
{
	delta_map_loc_[DIR_N] = -(int)width_ + 1;// North	(We have 'width - 1' number of corners in each row, so going north subtracts width - 1 from its current location)
	delta_map_loc_[DIR_E] = 1;			// East
	delta_map_loc_[DIR_S] = width_ - 1;	// South
	delta_map_loc_[DIR_W] = -1;		// West
	delta_map_loc_[DIR_NE] = delta_map_loc_[DIR_N] + delta_map_loc_[DIR_E];	// North-East
	delta_map_loc_[DIR_SE] = delta_map_loc_[DIR_S] + delta_map_loc_[DIR_E];	// South-East
	delta_map_loc_[DIR_SW] = delta_map_loc_[DIR_S] + delta_map_loc_[DIR_W];	// South-West
	delta_map_loc_[DIR_NW] = delta_map_loc_[DIR_N] + delta_map_loc_[DIR_W];	// North-West

	// Create the extra copies
	for (direction d = 0; d < 8; d++)
	{
		delta_map_loc_[d + 8] = delta_map_loc_[d];
		delta_map_loc_[d + 16] = delta_map_loc_[d];
	}
}
void SubgoalAA::IdentifySubgoals()
{
	for (cornerId i = 0; i < corner_locations_.size(); i++)
	{
		if (IsConvexCorner(i))
		{
			Subgoal s;
			s.loc = ToXYLoc(i);
			s.lvl = 1;
			s.has_extra_edge = false;
			subgoals_.push_back(s);

			corners_[i].sg_id = subgoals_.size()-1;
		}
		else
			corners_[i].sg_id = NON_SUBGOAL;
	}

	num_subgoals_ = subgoals_.size();
	num_global_subgoals_ = num_subgoals_;

#ifdef ANY_ANGLE_VERBOSE
	printf("%u corners, %u subgoals\n", corners_.size(), subgoals_.size());
#endif

	// Preallocate the extra subgoals for start and goal
	Subgoal s;
	s.lvl = 0;
	s.has_extra_edge = false;
	subgoals_.push_back(s);
	subgoals_.push_back(s);

}
void SubgoalAA::ComputeClearances()
{
	// North clearances
	for (int x = 0; x < (int)width_-1; x++)
	{
		int clearance = 0;
		for (int y = 0; y < (int)height_-1; y++)
		{
			cornerId id = ToCornerId(x, y);
			if (!(IsTraversable(NorthWestCell(id)) || IsTraversable(NorthEastCell(id))))
			{
				clearance = 0;
			}
			else
			{
				clearance ++;
			}

			SetClearance(id, DIR_N, clearance);

			if (IsSubgoal(id))
				clearance = 0;
		}
	}

	// South clearances
	for (int x = 0; x < (int)width_-1; x++)
	{
		int clearance = 0;
		for (int y = (int)height_-2; y >= 0; y--)
		{
			cornerId id = ToCornerId(x, y);
			if (!(IsTraversable(SouthWestCell(id)) || IsTraversable(SouthEastCell(id))))
			{
				clearance = 0;
			}
			else
			{
				clearance ++;
			}

			SetClearance(id, DIR_S, clearance);

			if (IsSubgoal(id))
				clearance = 0;
		}
	}

	// West clearances
	for (int y = 0; y < (int)height_-1; y++)
	{
		int clearance = 0;
		for (int x = 0; x < (int)width_-1; x++)
		{
			cornerId id = ToCornerId(x, y);
			if (!(IsTraversable(NorthWestCell(id)) || IsTraversable(SouthWestCell(id))))
			{
				clearance = 0;
			}
			else
			{
				clearance ++;
			}

			SetClearance(id, DIR_W, clearance);

			if (IsSubgoal(id))
				clearance = 0;
		}
	}

	// East clearances
	for (int y = 0; y < (int)height_-1; y++)
	{
		int clearance = 0;
		for (int x = (int)width_-2; x >= 0; x--)
		{
			cornerId id = ToCornerId(x, y);
			if (!(IsTraversable(NorthEastCell(id)) || IsTraversable(SouthEastCell(id))))
			{
				clearance = 0;
			}
			else
			{
				clearance ++;
			}

			SetClearance(id, DIR_E, clearance);

			if (IsSubgoal(id))
				clearance = 0;
		}
	}

	// Compute diagonal clearances
	for (int y = 0; y < (int)height_-1; y++)
	{
		for (int x = 0; x < (int)width_-1; x++)
		{
			cornerId id = ToCornerId(x, y);
			int clearance;

			// Northwest clearance;
			direction d = DIR_NW;

			if (!(IsTraversable(NorthWestCell(id))))
			{
				clearance = 0;
			}
			else
			{
				cornerId id2 = ToCornerId(x-1, y-1);
				if (IsSubgoal(id2))
					clearance = 1;
				else
				{
					clearance = GetClearance(id2, d);
					clearance++;
				}
			}

			SetClearance(id, d, clearance);

			// Northeast clearance;
			d = DIR_NE;

			if (!(IsTraversable(NorthEastCell(id))))
			{
				clearance = 0;
			}
			else
			{
				cornerId id2 = ToCornerId(x+1, y-1);
				if (IsSubgoal(id2))
					clearance = 1;
				else
				{
					clearance = GetClearance(id2, d);
					clearance++;
				}
			}

			SetClearance(id, d, clearance);
		}
	}

	for (int y = (int)height_-2; y >= 0; y--)
	{
		for (int x = 0; x < (int)width_-1; x++)
		{
			cornerId id = ToCornerId(x, y);
			int clearance;

			// Southwest clearance;
			direction d = DIR_SW;

			if (!(IsTraversable(SouthWestCell(id))))
			{
				clearance = 0;
			}
			else
			{
				cornerId id2 = ToCornerId(x-1, y+1);
				if (IsSubgoal(id2))
					clearance = 1;
				else
				{
					clearance = GetClearance(id2, d);
					clearance++;
				}
			}

			SetClearance(id, d, clearance);

			// Southeast clearance;
			d = DIR_SE;

			if (!(IsTraversable(SouthEastCell(id))))
			{
				clearance = 0;
			}
			else
			{
				cornerId id2 = ToCornerId(x+1, y+1);
				if (IsSubgoal(id2))
					clearance = 1;
				else
				{
					clearance = GetClearance(id2, d);
					clearance++;
				}
			}

			SetClearance(id, d, clearance);
		}
	}
}
void SubgoalAA::LinkSubgoals()
{
	std::vector<subgoalId> directHReachableNeigbors;

	int nEdges = 0;
	for (subgoalId s = 0; s < num_subgoals_; s++)
	{
		GetDirectHReachableSubgoals(subgoals_[s].loc, subgoals_[s].neighbors);
		nEdges += subgoals_[s].neighbors.size();

		subgoals_[s].neighbor_dist.clear();

		// Add the edge lengths as octile distance first, for partitioning. Finalize graph will convert them to euclidean distances after partitioning
		for (unsigned int i = 0; i < subgoals_[s].neighbors.size(); i++)
		{
			subgoals_[s].neighbor_dist.push_back(OctileDistanceSG(s, subgoals_[s].neighbors[i]));
		}
	}
}
void SubgoalAA::GetDirectHReachableSubgoals(xyLoc & from, std::vector<subgoalId> & subgoals)
{
	cornerId origin = ToCornerId(from);
	subgoals.clear();
	int clearance[9];

	// Get clearances
	for (direction d = 0; d < 8; d++)
	{
		clearance[d] = GetClearance(origin, d);

		cornerId c = origin + delta_map_loc_[d] * clearance[d];
		if (IsSubgoal(c) && clearance[d] != 0)
		{
			subgoals.push_back(corners_[c].sg_id);
			clearance[d]--;
		}
	}

	clearance[8] = clearance[0];	// Report North twice, to avoid using the modulus operation

	/* Now, explore the 8 areas that the cardinal and diagonal separators create.
	 * Each area is at most as large as a parallelogram whose dimensions are the
	 * clearances of its associated cardinal and diagonal directions.
	 */

	for (direction d = 1; d <= 7; d+=2)
	{
		for (direction c = d-1; c <= d+1; c+=2)
		{

			//int maxExt = GetClearance(from,c);
			int maxExt = clearance[c];
			cornerId loc = origin;

			for (int i = 1; i <= clearance[d]; i++)	// Go along the diagonal, from each corner on the diagonal ..
			{
				loc +=delta_map_loc_[d];
				int cardClearance = GetClearance(loc,c % 8);

				if (cardClearance <= maxExt)
				{
					maxExt = cardClearance;
					cornerId loc2 = loc + delta_map_loc_[c]*cardClearance;
					if (IsSubgoal(loc2))
					{
						subgoals.push_back(corners_[loc2].sg_id);
						maxExt--;
					}
				}
			}
		}
	}
}
void SubgoalAA::FinalizeGraph()
{
	// Set the number of buckets that will be used when creating descending paths to the goal
	buckets_.resize(graph_level_+1);	// From level 0 (for goal) to level graphLevel

	for (subgoalId s = 0; s < num_subgoals_; s++)
	{
		for (unsigned int i = 0; i < subgoals_[s].neighbors.size(); i++)
		{
			subgoals_[s].neighbor_dist[i] = EuclideanDistanceSG(s, subgoals_[s].neighbors[i]);
		}
		for (unsigned int i = 0; i < subgoals_[s].local_neighbors.size(); i++)
		{
			subgoals_[s].local_neighbor_dist[i] = EuclideanDistanceSG(s, subgoals_[s].local_neighbors[i]);
		}
	}
}

// Functions to create n-level subgoal graphs from simple subgoal graphs
bool SubgoalAA::AddLevel()
{

	std::vector<Subgoal> sgCopy = subgoals_;
	int nGlobalSubgoalsCopy = num_global_subgoals_;

	for (subgoalId s = 0; s < num_subgoals_; s++)
	{
		if (IsGlobalSubgoal(s))
		{
			subgoals_[s].lvl++;
		}
	}

	graph_level_++;
	std::vector<subgoalId> from, to;

	for (subgoalId s = 0; s < num_subgoals_; s++)
	{
		from.clear();
		to.clear();

		if (IsGlobalSubgoal(s) && CanDemoteSubgoal(s, from, to))
			DemoteSubgoal(s, from, to);

		/*
		for (subgoalId s = 0; s < nSubgoals; s++) {
			for (int i = 0; i < sg[s].neighbors.size(); i++) {
				for (int j = i+1; j < sg[s].neighbors.size(); j++)	{
					assert(sg[s].neighbors[i] != sg[s].neighbors[j]);
				}
			}
		}
		*/
	}

	// Highest level did not change or there are no global subgoals remaining
	if (num_global_subgoals_ == nGlobalSubgoalsCopy || num_global_subgoals_ == 0)
	{
#ifdef ANY_ANGLE_VERBOSE
		if (num_global_subgoals_ == 0)
			std::cout<<"No global subgoals remaining, going back a level"<<std::endl;
		else
			std::cout<<"Global subgoal graph cannot be partitioned any more"<<std::endl;
#endif
		num_global_subgoals_ = nGlobalSubgoalsCopy;
		subgoals_ = sgCopy;
		graph_level_ --;
		return false;
	}

	else
	{
		for (subgoalId s = 0; s < num_subgoals_; s++)
		{
			if (IsGlobalSubgoal(s))
			{
				subgoals_[s].local_neighbors.clear();
				subgoals_[s].local_neighbor_dist.clear();
			}
		}
		return true;
	}
}
bool SubgoalAA::CanDemoteSubgoal(subgoalId s, std::vector<subgoalId> & from, std::vector<subgoalId> & to)
{
	// Collect the neighbors and local neighbors of s into a new vector
	std::vector<subgoalId> neighbors = subgoals_[s].neighbors;

	for (unsigned int i = 0; i < subgoals_[s].local_neighbors.size(); i++)
		neighbors.push_back(subgoals_[s].local_neighbors[i]);

	for (unsigned int i = 0; i+1 < neighbors.size(); i++)	// Try to find a pair of subgoals that needs this one
		for (unsigned int j = i+1; j < neighbors.size(); j++)
		{
			if(IsNecessaryToConnect(s, neighbors[i], neighbors[j], from, to))
				return false;
		}

	return true;
}
void SubgoalAA::DemoteSubgoal(subgoalId s, std::vector<subgoalId> & from, std::vector<subgoalId> & to)
{

	// Since this is now a local subgoal, make all incoming edges local edges
	for (unsigned int i = 0; i < subgoals_[s].neighbors.size(); i++)
	{
		MakeEdgeLocal(subgoals_[s].neighbors[i], s);
	}

	// Add the necessary connections between surrounding subgoals
	for (unsigned int i = 0; i < from.size(); i++)
	{
		cost c = OctileDistanceSG(from[i], to[i]);
		if (CostOtherPath(s, from[i], to[i], c+1) > c + EPSILON)
		{

			if (!EdgeExists(from[i], to[i]))
			{
				AddOctileEdge(from[i], to[i]);

				if (!IsGlobalSubgoal(to[i]))
					MakeEdgeLocal(from[i], to[i]);
			}

			if (!EdgeExists(to[i], from[i]))
			{
				AddOctileEdge(to[i], from[i]);

				if (!IsGlobalSubgoal(from[i]))
					MakeEdgeLocal(to[i], from[i]);
			}
		}
	}

	DecrementLevel(s);
	num_global_subgoals_--;
}
cost SubgoalAA::CostOtherPath(subgoalId sg0, subgoalId sg1, subgoalId sg2, cost limit)
{
	assert(sg0 != sg1 && sg1 != sg2 && sg0 != sg2);

	// If sg2 is already pruned, make sure it is reachable by adding the relevant global-to-local edges
	if(!IsGlobalSubgoal(sg2))
	{
		for (unsigned int i = 0; i < subgoals_[sg2].neighbors.size(); i++)
		{
			subgoalId tempSg = subgoals_[sg2].neighbors[i];
			AddOctileEdge(tempSg, sg2);
		}
	}

	// Temporarily remove all the outgoing edges from sg, essentially removing it from the search
	std::vector<subgoalId> edges = subgoals_[sg0].neighbors;
	subgoals_[sg0].neighbors.clear();

	// Do the search and get the cost. Since we are only checking if we can do better
	// than costThrough, add a search limit to stop the search if we exceed costThrough

	cost alternateCost = SubgoalAStarSearch(sg1,sg2, limit);

	// Add back the outgoing edges from sg
	subgoals_[sg0].neighbors = edges;

	// If sg2 is already pruned, make sure to remove the global-to-local edges we added earlier
	if(!IsGlobalSubgoal(sg2))
	{
		for (unsigned int i = 0; i < subgoals_[sg2].neighbors.size(); i++)
		{
			subgoalId tempSg = subgoals_[sg2].neighbors[i];
			subgoals_[tempSg].neighbors.pop_back();
			subgoals_[tempSg].neighbor_dist.pop_back();
		}
	}

	return alternateCost;
}
bool SubgoalAA::IsNecessaryToConnect(subgoalId sg0, subgoalId sg1, subgoalId sg2, std::vector<subgoalId> & from, std::vector<subgoalId> & to)
{
	cost costThrough = OctileDistanceSG(sg1, sg0) + OctileDistanceSG(sg0, sg2);

	// If the length of path s1-s0-s2 is equal to the octile distance between s1 and s2, then s1 and s2 must be h-reachable (but not necesarily have LOS, so we skip this for now
//	if (costThrough < OctileDistanceSG(sg1, sg2) + EPSILON)
//		return false;

	if (LineOfSight(subgoals_[sg1].loc, subgoals_[sg2].loc))
	{
		// LOS check was successful for this pair, so we know we can add an edge between them, if necessary
		from.push_back(sg1);
		to.push_back(sg2);
		return false;
	}

	// Remove sg from the graph, do an A* search from sg1 to sg2. If the path found is larger then h(sg1,sg) + h(sg,sg2), then sg is necessary to connect
	if (CostOtherPath(sg0, sg1,sg2, costThrough+1) <= costThrough + EPSILON)
		return false;

	return true;
}

// SEARCH
cost SubgoalAA::FindXYLocPath(xyLoc from, xyLoc to, std::vector<xyLoc> &path)
{
    switch (search_method_) {

    case SUBGOAL_A_STAR:
        return FindSubgoalAnyAnglePath(from, to, path, false);

    case SUBGOAL_THETA_STAR:
        return FindSubgoalAnyAnglePath(from, to, path, true);

    default:
        return INFINITE_COST;
    }
}

// Functions to answer queries by setting up the graph for search and calling relevant search functions
cost SubgoalAA::FindSubgoalAnyAnglePath(xyLoc start, xyLoc goal, std::vector<xyLoc> & thePath, bool theta)
{
	thePath.clear();

	// First, check if there is a direct path between start and goal
	if (LineOfSight(start, goal))
	{
		thePath.push_back(start);
		thePath.push_back(goal);
		return EuclideanDistance(start,goal);
	}


	// If not, find a high level path between them
	subgoalId sgStart;
	subgoalId sgGoal;
	std::vector<subgoalId> startDHR;
	std::vector<subgoalId> goalDHR;
	cost pathCost;

	// Add the relevant edges to the graph
	ConnectStartAndGoalToGraph(start, goal, sgStart, sgGoal, startDHR, goalDHR);

	if (theta)
		pathCost = SubgoalThetaStarSearch(sgStart, sgGoal);
	else
		pathCost = SubgoalAStarSearch(sgStart, sgGoal);

	ExtractPath(sgStart, sgGoal, thePath);

	// Restore the original graph
	for (level l = 0; l <= graph_level_; l++)
	{
		while(!buckets_[l].empty())
		{
			subgoalId s = buckets_[l].back();
			buckets_[l].pop_back();
			subgoals_[s].has_extra_edge = false;
			subgoals_[s].neighbors.pop_back();
			subgoals_[s].neighbor_dist.pop_back();
		}
	}

	return pathCost;
}
void SubgoalAA::ConnectStartAndGoalToGraph(
		xyLoc & start, xyLoc & goal, subgoalId & sgStart, subgoalId & sgGoal,
		std::vector<subgoalId> & startDHR, std::vector<subgoalId> & goalDHR)
{


	cornerId startCorner = ToCornerId(start);
	cornerId goalCorner = ToCornerId(goal);

	if (IsSubgoal(startCorner))	// If the start location is already a subgoal, use it
	{
		sgStart = corners_[startCorner].sg_id;
	}

	else	// Create a new subgoal for it, with id nSubgoals (at the end of all the actual subgoals)
	{
		sgStart = num_subgoals_;
		subgoals_[sgStart].loc = start;
		GetDirectHReachableSubgoals(start, startDHR);

		subgoals_[sgStart].neighbors.clear();
		subgoals_[sgStart].neighbor_dist.clear();

		for (unsigned int i = 0; i < startDHR.size(); i++)
			AddEuclideanEdge(sgStart, startDHR[i]);
	}

	if (IsSubgoal(goalCorner))	// If the goal location is already a subgoal, use it
	{
		sgGoal = corners_[goalCorner].sg_id;
	}

	else	// Create a new subgoal for it, with id nSubgoals + 1 (at the end of all the actual subgoals and start)
	{
		sgGoal = num_subgoals_+1;
		subgoals_[sgGoal].loc = goal;
		GetDirectHReachableSubgoals(goal, goalDHR);

		subgoals_[sgGoal].neighbors.clear();
		subgoals_[sgGoal].neighbor_dist.clear();

		for (unsigned int i = 0; i < goalDHR.size(); i++)
		{
			AddEuclideanEdge(sgGoal, goalDHR[i]);
		}
	}

	// As a hack for the next step, increase its number of neighbors by 1 (the extra neighbor will never be checked)
	// We do this because the next step avoids checking the extra neighbor by assuming the number of neigbors is -1

	AddEdge(sgGoal, sgGoal, 0);
	subgoals_[sgGoal].has_extra_edge = true;
	subgoals_[sgGoal].g_val = 0;

	// Add new edges to the graph, back from the goal, so that the search is optimal
	buckets_[subgoals_[sgGoal].lvl].push_back(sgGoal);

	/* Bucket initially contains only the goal
	 * In the first iteration, only edges to higher level subgoals are considered
	 * When a higher level successor of a subgoal in a bucket is found, it is also added to the bucket, if it is not already included
	 * If it is already in the bucket, its g-value and extra edge is updated if the new g-value is lower
	 * A subgoal that is put in a bucket is not removed until the end of the search (that's how we keep track of subgoals with extra edges)
	 *
	 * In the second iteration, we iterate over all the subgoals already placed in the buckets in the first iteration
	 * This time, we expand its local successors
	 * We add these new local-successors to the bucket, but we do not expand them again
	 */

	// First iteration
	for (level l = subgoals_[sgGoal].lvl; l < graph_level_; l++)	// Exclude the graphLevel subgoals, they are only for keeping track of the extra edges so that we can delete them later
	{
		for (unsigned int i = 0; i < buckets_[l].size(); i++)
		{
			subgoalId sg1 = buckets_[l][i];

			for (unsigned int j = 0; j < subgoals_[sg1].neighbors.size()-1; j++)	// -1 because we do not want to check the extra edge
			{														// Note that we are always iterating the neighbors of a subgoal with an extra edge (including the goal, with the hack above)
				subgoalId sg2 = subgoals_[sg1].neighbors[j];
				cost edgeCost = subgoals_[sg1].neighbor_dist[j];

				if(!subgoals_[sg2].has_extra_edge)	// Has no other extra edge, just make sg2->sg1 the extra edge
				{
					AddEdge(sg2, sg1, edgeCost);
					subgoals_[sg2].has_extra_edge = true;
					subgoals_[sg2].g_val = subgoals_[sg1].g_val + edgeCost;

					buckets_[subgoals_[sg2].lvl].push_back(sg2);	// Add it to the relevant bucket
				}

				else	// Already has an extra edge, make sg2->sg the extra edge only if it yields a shorter path
				{
					cost newCost = subgoals_[sg1].g_val + edgeCost;
					if (newCost + EPSILON < subgoals_[sg2].g_val)
					{
						subgoals_[sg2].neighbors.back() = sg1;
						subgoals_[sg2].neighbor_dist.back() = edgeCost;
						subgoals_[sg2].g_val = newCost;
					}
				}
			}
		}
	}

//*
	// Second iteration
	for (level l = subgoals_[sgGoal].lvl<subgoals_[sgStart].lvl?subgoals_[sgStart].lvl:subgoals_[sgGoal].lvl; l < graph_level_; l++)	// We skip the buckets with level < start's level,
	{																										// since lower level local edges are irrelevant
		int maxSize = buckets_[l].size();
		for (int i = 0; i < maxSize; i++)
		{
			subgoalId sg1 = buckets_[l][i];

			for (unsigned int j = 0; j < subgoals_[sg1].local_neighbors.size(); j++)	// Iterate over its local edges
			{
				subgoalId sg2 = subgoals_[sg1].local_neighbors[j];
				cost edgeCost = subgoals_[sg1].local_neighbor_dist[j];

				if(!subgoals_[sg2].has_extra_edge)	// Has no other extra edge, just make sg2->sg1 the extra edge
				{
					AddEdge(sg2, sg1, edgeCost);
					subgoals_[sg2].has_extra_edge = true;
					subgoals_[sg2].g_val = subgoals_[sg1].g_val + edgeCost;

					buckets_[subgoals_[sg2].lvl].push_back(sg2);	// Add it to the relevant bucket
				}

				else	// Already has an extra edge, make sg2->sg the extra edge only if it yields a shorter path
				{
					cost newCost = subgoals_[sg1].g_val + edgeCost;
					if (newCost < subgoals_[sg2].g_val)
					{
						subgoals_[sg2].neighbors.back() = sg1;
						subgoals_[sg2].neighbor_dist.back() = edgeCost;
						subgoals_[sg2].g_val = newCost;
					}
				}
			}
		}
	}
//*/
}

// Functions to search the graph
cost SubgoalAA::SubgoalAStarSearch(subgoalId start, subgoalId goal, cost limit)
{
	if (search_ >= MAX_SEARCH)
		ResetSearch();

	// Initialize search
	search_++;

	open_list_.clear();
	GenerateState(start,goal);
	GenerateState(goal,goal);

	subgoals_[start].g_val = 0;
	AddToOpen(start);

	subgoals_[goal].g_val = limit;

	while(!open_list_.empty() && subgoals_[goal].g_val > GetMin().f_val + EPSILON)
	{
		// Expand the state with the minimum f-value
		subgoalId curr = GetMin().id;
		PopMin();

		#ifdef ANY_ANGLE_STATISTICS
			num_expansions_++;
		#endif

		for (unsigned int i = 0; i < subgoals_[curr].neighbors.size(); i++)
		{
			subgoalId succ = subgoals_[curr].neighbors[i];
			GenerateState(succ, goal);

			if (subgoals_[succ].list != CLOSED_LIST)
			{
				cost newGVal = subgoals_[curr].g_val + subgoals_[curr].neighbor_dist[i];

				if (newGVal + EPSILON < subgoals_[succ].g_val)
				{
					subgoals_[succ].g_val = newGVal;
					subgoals_[succ].parent = curr;
					AddToOpen(succ);
				}
			}
		}
	}

#ifdef ANY_ANGLE_ASSERTIONS
    // Make sure that the heap property is satisfied for the open list
    for (int i = 0; i < open_list_.size(); i++) {
        int child1 = (i << 1) + 1;
        int child2 = (i << 1) + 2;

        if (child2 < open_list_.size()) {
            assert(open_list_[i].f_val < open_list_[child1].f_val + 100*EPSILON);
            assert(open_list_[i].f_val < open_list_[child2].f_val + 100*EPSILON);
        }
    }
#endif

	return subgoals_[goal].g_val;
}
cost SubgoalAA::SubgoalThetaStarSearch(subgoalId start, subgoalId goal, cost limit)
{
	if (search_ >= MAX_SEARCH)
		ResetSearch();

	search_++;

	open_list_.clear();
	GenerateState(start,goal);
	GenerateState(goal,goal);

	subgoals_[start].g_val = 0;
	subgoals_[start].parent = start;	// Just a trick: when 'start' is expanded for the first time,
									// the grandparents of 'start's successors will be 'start' as well
	AddToOpen(start);

	subgoals_[goal].g_val = limit;

	//while(!theHeap.empty() && sg[goal].gVal > GetMin().fVal + EPSILON)
	while(!open_list_.empty() && subgoals_[goal].list != CLOSED_LIST)
	{
		subgoalId curr = GetMin().id;
		PopMin();

		#ifdef ANY_ANGLE_STATISTICS
			num_expansions_++;
		#endif

		for (unsigned int i = 0; i < subgoals_[curr].neighbors.size(); i++)
		{
			subgoalId succ = subgoals_[curr].neighbors[i];
			GenerateState(succ, goal);

			/* TODO: If we allow re-expansions, the quality of the path seems to be better, but still not good enough. What if:
			 * instead of the g and h value of the current node we store the following:
			 * - g_p = g-value of its parent
			 * - h_p = h-value of its parent
			 * - g = g-value of the node itself
			 * The node's f-value will be f_p = g_p + h_p, but we choose its parent as the one that minimizes its own g-value
			 */

			if (subgoals_[succ].list != CLOSED_LIST)
			{
				cost newGVal;

//*
				subgoalId newParent;
				subgoalId grandParent = subgoals_[curr].parent;
				// First, check if there is LOS from the grandparent
				if (LineOfSight(subgoals_[grandParent].loc, subgoals_[succ].loc))
				{
					newParent = grandParent;
					newGVal = subgoals_[newParent].g_val + EuclideanDistanceSG(newParent,succ);
				}
				else
				{
					newParent = curr;
					newGVal = subgoals_[curr].g_val + subgoals_[curr].neighbor_dist[i];
				}

/*/
				//  New idea: Check if there is line of sight with grand parent. If there is, check if there is line of sight with grand grand parent. If there is, check grand grand grand parent and so on..

				subgoalId newParent = curr;
				subgoalId potentialParent = sg[newParent].parent;
				while(newParent != potentialParent && LineOfSight(sg[potentialParent].loc, sg[succ].loc))
				{
					newParent = potentialParent;
					potentialParent = sg[potentialParent].parent;
				}
				newGVal = sg[newParent].gVal + EuclideanDistanceSG(newParent,succ);
//*/
				if (newGVal + EPSILON < subgoals_[succ].g_val)
				{
					subgoals_[succ].g_val = newGVal;
					subgoals_[succ].parent = newParent;
					AddToOpen(succ);
				}
			}
		}
	}

#ifdef ANY_ANGLE_ASSERTIONS
    // Make sure that the heap property is satisfied for the open list
    for (int i = 0; i < open_list_.size(); i++) {
        int child1 = (i << 1) + 1;
        int child2 = (i << 1) + 2;

        if (child2 < open_list_.size()) {
            assert(open_list_[i].f_val < open_list_[child1].f_val + 100*EPSILON);
            assert(open_list_[i].f_val < open_list_[child2].f_val + 100*EPSILON);
        }
    }
#endif

	return subgoals_[goal].g_val;
}
cost SubgoalAA::ExtractPath(subgoalId start, subgoalId goal, std::vector<xyLoc> & thePath)
{
	thePath.clear();
	// Extract path
	if (subgoals_[goal].g_val < INFINITE_COST)
	{
		subgoalId curr = goal;
		while (curr != start)
		{
			thePath.push_back(subgoals_[curr].loc);
			curr = subgoals_[curr].parent;
		}

		thePath.push_back(subgoals_[curr].loc);
		std::reverse(thePath.begin(), thePath.end());
	}

	return subgoals_[goal].g_val;
}

// Search utility functions
void SubgoalAA::ResetSearch()
{
	// Set last search and generated values to 0, so that when the search is incremented, all states are un-generated
	search_ = 0;
	for (subgoalId i = 0; i < subgoals_.size(); i++)	// Includes the additional start and goal
		subgoals_[i].generated = 0;
}
void SubgoalAA::GenerateState(subgoalId s, subgoalId goal)
{	if (subgoals_[s].generated < search_)
	{
		#ifdef ANY_ANGLE_STATISTICS
			num_generated_++;
		#endif
		subgoals_[s].generated = search_;
		subgoals_[s].h_val = EuclideanDistanceSG(s, goal);
		subgoals_[s].g_val = INFINITE_COST;
		subgoals_[s].list = NO_LIST;
	}
}


// Heap operations.
void SubgoalAA::AddToOpen(const subgoalId id)
{
    // If it is already in the open list, update its location.
	if (subgoals_[id].list == OPEN_LIST) {
		int i = subgoals_[id].heap_index;
		open_list_[i].g_val = subgoals_[id].g_val;
		open_list_[i].f_val = subgoals_[id].g_val + subgoals_[id].h_val;
		PercolateUp(i);
	}
	// Otherwise, add it to the open list.
	else {
		subgoals_[id].list = OPEN_LIST;
		subgoals_[id].heap_index = open_list_.size();
		open_list_.push_back(HeapElement(id, subgoals_[id].g_val, subgoals_[id].g_val + subgoals_[id].h_val));
		PercolateUp(open_list_.size()-1);
	}
}
void SubgoalAA::PopMin()
{
	subgoals_[open_list_[0].id].list = CLOSED_LIST;
	open_list_[0] = open_list_.back();
	subgoals_[open_list_[0].id].heap_index = 0;
	open_list_.pop_back();
	PercolateDown(0);
}
void SubgoalAA::PercolateUp(int index)
{
	HeapElement elem = open_list_[index];
	int parent = (index-1) >> 1;

	while(index > 0 && elem < open_list_[parent])
	{
		open_list_[index] = open_list_[parent];
	    subgoals_[open_list_[index].id].heap_index = index;

		index = parent;
		parent = (index-1) >> 1;

#ifdef ANY_ANGLE_STATISTICS
        num_percolations_++;
#endif
	}

	open_list_[index] = elem;
	subgoals_[elem.id].heap_index = index;
}
void SubgoalAA::PercolateDown(int index)
{
	HeapElement elem = open_list_[index];
	int maxSize = open_list_.size();

	while(true)
	{
		int child1 = (index << 1) + 1;
		if (child1 >= maxSize)
			break;

		int child2 = child1+1;

		// If the first child has smaller key
		if (child2 == maxSize || open_list_[child1] < open_list_[child2])
		{
			if (open_list_[child1] < elem)
			{
				open_list_[index] = open_list_[child1];
				subgoals_[open_list_[index].id].heap_index = index;
				index = child1;

#ifdef ANY_ANGLE_STATISTICS
				num_percolations_++;
#endif
			}
			else
				break;
		}

		else if (open_list_[child2] < elem)
		{
			open_list_[index] = open_list_[child2];
			subgoals_[open_list_[index].id].heap_index = index;
			index = child2;

#ifdef ANY_ANGLE_STATISTICS
			num_percolations_++;
#endif
		}
		else
			break;
	}

	open_list_[index] = elem;
	subgoals_[elem.id].heap_index = index;
}

// Edge functions
void SubgoalAA::AddEuclideanEdge(subgoalId sg1, subgoalId sg2)
{
	subgoals_[sg1].neighbors.push_back(sg2);
	subgoals_[sg1].neighbor_dist.push_back(EuclideanDistanceSG(sg1, sg2));
}

void SubgoalAA::AddOctileEdge(subgoalId sg1, subgoalId sg2)
{
	subgoals_[sg1].neighbors.push_back(sg2);
	subgoals_[sg1].neighbor_dist.push_back(OctileDistanceSG(sg1, sg2));
}

void SubgoalAA::AddEdge(subgoalId sg1, subgoalId sg2, cost c)
{
	subgoals_[sg1].neighbors.push_back(sg2);
	subgoals_[sg1].neighbor_dist.push_back(c);
}

void SubgoalAA::MakeEdgeLocal(subgoalId sg1, subgoalId sg2)
{
	for (int i = subgoals_[sg1].neighbors.size()-1; i >= 0; i--)
	{
		if (subgoals_[sg1].neighbors[i] == sg2)
		{
			// Copy the information to the local neighbors vector
			subgoals_[sg1].local_neighbors.push_back(sg2);
			subgoals_[sg1].local_neighbor_dist.push_back(subgoals_[sg1].neighbor_dist[i]);

			subgoals_[sg1].neighbors[i] = subgoals_[sg1].neighbors.back();
			subgoals_[sg1].neighbor_dist[i] = subgoals_[sg1].neighbor_dist.back();

			subgoals_[sg1].neighbors.pop_back();
			subgoals_[sg1].neighbor_dist.pop_back();

			return;
		}
	}
}

bool SubgoalAA::EdgeExists(subgoalId sg1, subgoalId sg2)
{
	for (int i = subgoals_[sg1].neighbors.size()-1; i >= 0; i--)
	{
		if (subgoals_[sg1].neighbors[i] == sg2)
			return true;
	}

	for (int i = subgoals_[sg1].local_neighbors.size()-1; i >= 0; i--)
	{
		if (subgoals_[sg1].local_neighbors[i] == sg2)
			return true;
	}

	return false;
}

// STATISTICS
#ifdef ANY_ANGLE_STATISTICS
void SubgoalAA::PrintAdditionalStatistics(AnyAngleStatistics* stat)
{
	if(stat == NULL)
		return;

	int num_edges = 0;
	int num_global_edges = 0;
	int num_ascending_edges = 0;
	int num_local_edges = 0;

	for (subgoalId s = 0; s < num_subgoals_; s++)
	{
		num_edges += subgoals_[s].neighbors.size() + subgoals_[s].local_neighbors.size();

		if (IsGlobalSubgoal(s))
		{
			num_global_edges += subgoals_[s].neighbors.size();
		}
		else
		{
			num_ascending_edges += subgoals_[s].neighbors.size();
			num_local_edges += subgoals_[s].local_neighbors.size();
		}
	}

	stat->ReportIntStatistic("Graph level", graph_level_);
	stat->ReportDoubleStatistic("Time to construct SSG (ms)", ssg_construction_time_ * 1000.0);
	stat->ReportDoubleStatistic("Time to partition the SSG (ms)", partitioning_time_ * 1000.0);

	stat->AddRemark("");	// endline

	stat->ReportIntStatistic("Number of corners", corners_.size());
	stat->ReportIntStatistic("Number of subgoals", num_subgoals_);
	stat->ReportIntStatistic("Number of global subgoals", num_global_subgoals_);

	stat->ReportIntStatistic("Total number of directed edges", num_edges);
	stat->ReportIntStatistic("Total number of directed global edges", num_global_edges);
	stat->ReportIntStatistic("Total number of directed ascending edges", num_ascending_edges);
	stat->ReportIntStatistic("Total number of directed local edges", num_local_edges);

	stat->ReportDoubleStatistic("Branching factor of the subgoal graph", num_edges / (double)num_subgoals_);
	stat->ReportDoubleStatistic("Branching factor of the global subgoal graph", num_global_edges / (double)num_global_subgoals_);

	stat->AddRemark("");	// endline
}
#endif

// DISPLAY
#ifdef ANY_ANGLE_RUNNING_IN_HOG
void SubgoalAA::ProcessDebugLoc(const xyLoc loc) {

    if (IsSubgoal(loc)) {
       std::cout<<"Debug loc: "<<loc<<std::endl;
       const Subgoal* s = &subgoals_[ToSubgoalId(loc)];
       for (unsigned int i = 0; i < s->neighbors.size(); i++) {
           xyLoc loc2 = subgoals_[s->neighbors[i]].loc;
           std::cout<<loc2<<"\t"<<s->neighbor_dist[i]<<"\t"<<EuclideanDistance(loc, loc2)<<std::endl;
       }
    }
    else
        std::cout<<"Not a subgoal"<<std::endl;
/*
    printf("Debug point: (%d, %d)\n", loc.x, loc.y);
    cornerId id = ToCornerId(loc);
    printf("Clearances: N(%u), NE(%u), E(%u), SE(%u), S(%u), SW(%u), W(%u), NW(%u)\n",
            GetClearance(id, DIR_N),  GetClearance(id, DIR_NE),  GetClearance(id, DIR_E),  GetClearance(id, DIR_SE),
            GetClearance(id, DIR_S),  GetClearance(id, DIR_SW),  GetClearance(id, DIR_W),  GetClearance(id, DIR_NW));
*/
}

void SubgoalAA::VisualizeAlgorithm(const MapEnvironment *env)
{
	bool display_subgoals = 0;
	bool display_edges = 0;
	bool display_global_only = 0;
    bool display_search_tree = 0;
    bool display_debug_loc = 1;

    if (display_debug_loc) {
        env->SetColor(1,0,0);
        DrawPoint(env, debug_loc_);

        std::vector<subgoalId> dhr_subgoals;
        GetDirectHReachableSubgoals(debug_loc_, dhr_subgoals);

        for (unsigned int i = 0; i < dhr_subgoals.size(); i++)
        {
            DrawLine(env, debug_loc_, subgoals_[dhr_subgoals[i]].loc);
            DrawPoint(env, subgoals_[dhr_subgoals[i]].loc);
        }
    }

	env->SetColor(0,0,1);

	for (subgoalId i = 0; i < num_subgoals_; i++) //  && i < 1000
	{
		if (display_subgoals && (!display_global_only || IsGlobalSubgoal(i)))
			DrawPoint(env, subgoals_[i].loc);

		for (unsigned int j = 0; j < subgoals_[i].neighbors.size(); j++)
		{
			subgoalId i2 = subgoals_[i].neighbors[j];
			if (i > i2)	// Do not draw the same edge twice, from different corners
				if (display_edges && (!display_global_only || (IsGlobalSubgoal(i) && IsGlobalSubgoal(i2))))
					DrawLine(env, subgoals_[i].loc, subgoals_[i2].loc);
		}

		if (display_search_tree)
		{

			if (search_ != 0 && subgoals_[i].generated == search_)
			{
				subgoalId parent = subgoals_[i].parent;
				DrawLine(env, subgoals_[i].loc, subgoals_[parent].loc);
			}
		}
	}
}
#endif

