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

#include "BlockAStar.h"

#ifdef ANY_ANGLE_RUNNING_IN_HOG
BlockAStar::BlockAStar(MapEnvironment *env, int _blockSize)
	:AnyAngleAlgorithm::AnyAngleAlgorithm(env)
{
	block_size_ = _blockSize;
	BlockAStarInitialize();
}
#endif
BlockAStar::BlockAStar(std::vector<bool> &bits, int _width, int _height, int _blockSize)
	:AnyAngleAlgorithm::AnyAngleAlgorithm(bits, _width, _height)
{
	block_size_ = _blockSize;
	BlockAStarInitialize();
}
BlockAStar::~BlockAStar()
{
	for (unsigned int b = 0; b < blocks_.size(); b++){
		blocks_[b].clear();
	}
}
void BlockAStar::Block::clear()
{
	if (dist != NULL)
	{
		for (int i = 0; i < visibility_corners.size(); i++)
			delete [] dist[i];
		delete [] dist;
	}
	if (next != NULL)
	{
		for (unsigned int i = 0; i < visibility_corners.size(); i++)
			delete [] next[i];
		delete [] next;
	}

	if (last_g_val != NULL)
		delete [] last_g_val;
}
int BlockAStar::Block::getCornerIndex(cornerId c)
{
	for (unsigned int i = 0; i < visibility_corners.size(); i++)
		if (visibility_corners[i] == c)
			return i;
	return -1;
}

void BlockAStar::BlockAStarInitialize()
{
	corners_.resize(corner_locations_.size());

	block_cell_width_ = block_size_ - 1;

	// -2 to remove padding. (-1 mod blockCellWidth) + 1 to determine the number of blocks necessary to cover the whole row (or column) of a map
	map_block_width_ = ((width_ - 3)/block_cell_width_) + 1;
	map_block_height_ = ((height_ - 3)/block_cell_width_) + 1;

	for(unsigned int x = 0; x < map_block_width_; x++)
		for (unsigned int y = 0; y < map_block_height_; y++)
			InitializeBlock(xyLoc(x*block_cell_width_, y*block_cell_width_));	// Initialize the block specified by the top left corner and add it to the blocks list

	ResetSearch();

#ifdef ANY_ANGLE_STATISTICS
    SetupStatistics();
#endif
}
void BlockAStar::InitializeBlock(xyLoc l)
{
	blockId id = blocks_.size();
	blocks_.push_back(Block(l));

	// Don't go over the map boundaries
	int xMax = width_-2 < l.x + block_cell_width_ ? width_-2 : l.x + block_cell_width_;
	int yMax = height_-2 < l.y + block_cell_width_ ? height_-2 : l.y + block_cell_width_;

	std::vector<cornerId> convexCorners;	// Gathered separately from the boundary corners, because we want to append these after the boundary corners

	// Go over all the corners in the block and identify the visibility and boundary corners
	for (int x = l.x; x <= xMax; x++)
		for (int y = l.y; y <= yMax; y++)
		{
			cornerId c = corner_ids_[x][y];

			if (c != -1)	// If the corner is a valid corner
			{
				corners_[c].belongs_to.push_back(id);		// Link corner to block
				corners_[c].index_in_block.push_back(-1);		// It might not have an index if it is not a corner or visibility cell

				if (OnBlockEdge(c))
				{
					blocks_[id].visibility_corners.push_back(c);
				}

				else if (IsConvexCorner(c))
					convexCorners.push_back(c);

			}
		}

	// We want the vCorners array to first contain all the boundary corners, then the internal convex corners
	blocks_[id].num_boundary = blocks_[id].visibility_corners.size();
	for (unsigned int i = 0; i < convexCorners.size(); i++)
	{
		blocks_[id].visibility_corners.push_back(convexCorners[i]);
	}

	for (unsigned int i = 0; i < blocks_[id].visibility_corners.size(); i++)
		corners_[blocks_[id].visibility_corners[i]].index_in_block.back() = i;

	// Allocate space for the matrices
	unsigned int nB = blocks_[id].num_boundary;
	unsigned int nV = blocks_[id].visibility_corners.size();

	if (nB == 0)	// No corners
		return;

	cost** dist = new cost*[nV];	// nB x nB	- but use nV x nV for now, for convenience
	for (unsigned int i = 0; i < nV; i++)
		dist[i] = new cost[nV];

	char** next = new char*[nV];	// nV x nB	- but use nV x nV for now, for convenience
	for (unsigned int i = 0; i < nV; i++)
		next[i] = new char[nV];

	bool** vis = new bool*[nV];		// nV x nV
	for (unsigned int i = 0; i < nV; i++)
		vis[i] = new bool[nV];

	// We have to calculate pairwise distances for all visibility cells, not just the boundary corners extDist = extended distance matrix
	cost** extDist = new cost*[nV];	// nV x nV
	for (unsigned int i = 0; i < nV; i++)
		extDist[i] = new cost[nV];

	// Fill the visibility matrix and initialize the extDist matrix
	for (unsigned int i = 0; i < nV; i++)
	{
		vis[i][i] = true;	// A corner is visible from itself
		extDist[i][i] = 0;	// A corner is trivially reachable from itself
		for (unsigned int j = 0; j < i; j++)	// Go over potential neighbors
		{
			vis[i][j] = LineOfSight(blocks_[id].visibility_corners[i], blocks_[id].visibility_corners[j]);	// Determine if it is a neighbor

			if(vis[i][j])	// Set its cost
				extDist[i][j] = EuclideanDistance(blocks_[id].visibility_corners[i], blocks_[id].visibility_corners[j]);
			else
				extDist[i][j] = INFINITE_COST;

			vis[j][i] = vis[i][j];	// Copy the symmetric information
			extDist[j][i] = extDist[i][j];
		}
	}

	// Apply Floyd-Warshall on the extDist matrix
	for (unsigned int k = 0; k < nV; k++)
		for (unsigned int i = 0; i< nV; i++)
			for (unsigned int j = 0; j < nV; j++)
				if (extDist[i][k] + extDist[k][j] < extDist [i][j])
					extDist[i][j] = extDist[i][k] + extDist [k][j];

	// Fill the next pointers
	for (unsigned int i = 0; i< nV; i++)	// For each pair i,j
		for (unsigned int j = 0; j < nV; j++)
		{
			if(extDist[i][j] == INFINITE_COST)	// if j is not reachable from i
				next[i][j] = -1;	// Indicates non-reachability
			else if (i == j)
				next[i][j] = i;
			else
			{
				//cost best = INFINITE_COST;
				for (unsigned int k = 0; k < nV; k++)
					if (vis[i][k] && i != k)	// Go over i's neighbors k
						if (extDist[i][k] + extDist[k][j] - EPSILON <= extDist[i][j])	// Find the best one
						{
							//best = extDist[i][k] + extDist[k][j];
							next[i][j] = k;
						}
			}
		}

	// Copy the relevant part of extDist to dist
	for (unsigned int i = 0; i < nV; i++)		// used to be nB x nB, extended to nV x nV for convenience
		for (unsigned int j = 0; j < nV; j++)
			dist[i][j] = extDist[i][j];

	// Delete the extended dist matrix to avoid memory leaks
	for (unsigned int i = 0; i < nV; i++)
		delete [] extDist[i];
	delete [] extDist;

	// Delete the vis matrix
	for (unsigned int i = 0; i < nV; i++)
		delete [] vis[i];
	delete [] vis;

	blocks_[id].last_g_val = new cost[nB];

	// Pass along the information to the block
	blocks_[id].dist = dist;
	blocks_[id].next = next;
}

void BlockAStar::GetBoundaryCorners(cornerId c, std::vector<cornerId> & bCorners, std::vector<cost> & bCosts, std::vector<cornerId> & bNext)
{
	// Given a corner that belongs in the interior of a block, compute the shortest paths to the blocks boundary corners

	// c is in the interior of a block, get the block that it belongs to
	Block* b = &blocks_[corners_[c].belongs_to[0]];

	// Set the g-values of all vCorners in the block to infinity
	for (unsigned int i = 0; i < b->visibility_corners.size(); i++)
		corners_[b->visibility_corners[i]].g_val = INFINITE_COST;

	// Find which vCorners are visible from c
	std::vector<int> vCornerIndices;
	if (IsConvexCorner(c))	// If c itself is a visibility corner in the block
		vCornerIndices.push_back(b->getCornerIndex(c));

	else	// Find which vCorners are visible from c
		for (unsigned int i = 0; i < b->visibility_corners.size(); i++)
			if (LineOfSight(c, b->visibility_corners[i]))
				vCornerIndices.push_back(i);

	// We have determined the visible corners, now use the block's data to figure out the rest
	for (unsigned int i = 0; i < vCornerIndices.size(); i++)
	{
		int ind = vCornerIndices[i];
		cornerId vCorner = b->visibility_corners[ind];

		cost baseGVal = EuclideanDistance(c, vCorner);

		for (unsigned int j = 0; j < b->num_boundary; j++)
		{
			cornerId bCorner = b->visibility_corners[j];
			cost newGVal = baseGVal + b->dist[ind][j];
			if (newGVal < corners_[bCorner].g_val)
			{
				corners_[bCorner].g_val = newGVal;
				corners_[bCorner].parent = vCorner;
			}
		}
	}

	for (unsigned int i = 0; i < b->num_boundary; i++)
	{
		cornerId bCorner = b->visibility_corners[i];
		if (corners_[bCorner].g_val < INFINITE_COST)
		{
			bCorners.push_back(bCorner);
			bCosts.push_back(corners_[bCorner].g_val);
			bNext.push_back(corners_[bCorner].parent);
		}
	}
}

void BlockAStar::GenerateBlock(blockId b, cornerId goal)
{
	if (blocks_[b].generated < search_)
	{
		#ifdef ANY_ANGLE_STATISTICS
			num_generated_++;
		#endif

		blocks_[b].generated = search_;
		blocks_[b].g_val = INFINITE_COST;
		blocks_[b].h_val = INFINITE_COST;
		blocks_[b].list = NO_LIST;
		blocks_[b].is_goal_block = false;

		for(unsigned int i = 0; i < blocks_[b].num_boundary; i++)
		{
			blocks_[b].last_g_val[i] = INFINITE_COST;
			GenerateCorner(blocks_[b].visibility_corners[i], goal);
		}
	}
}
void BlockAStar::GenerateCorner(cornerId c, cornerId goal)
{
	if (corners_[c].generated < search_)
	{
		corners_[c].generated = search_;
		corners_[c].h_val = HeuristicDistance(c, goal);
		corners_[c].g_val = INFINITE_COST;
	}
}
blockId BlockAStar::GetCommonBlock(cornerId c1, cornerId c2)
{
	for(unsigned int i = 0; i < corners_[c1].belongs_to.size(); i++)
		for(unsigned int j = 0; j < corners_[c2].belongs_to.size(); j++)
			if(corners_[c1].belongs_to[i] == corners_[c2].belongs_to[j])
				return corners_[c1].belongs_to[i];

	return -1;
}
cost BlockAStar::FindPathWithinBlock(blockId b, cornerId start, cornerId goal, std::vector<xyLoc> & thePath)
{
	if (LineOfSight(start, goal))
	{
		thePath.push_back(corner_locations_[start]);
		thePath.push_back(corner_locations_[goal]);
		return EuclideanDistance(start, goal);
	}

	Block* bp = &blocks_[b];

	// Identify the visible vertices from the start or the goal
	std::vector<cornerId> vStart, vGoal;

	if (OnBlockEdge(start) || IsConvexCorner(start))
			vStart.push_back(start);
	else
		for (unsigned int i = 0; i < bp->visibility_corners.size(); i++)
			if (LineOfSight(start,bp->visibility_corners[i]))
				vStart.push_back(bp->visibility_corners[i]);

	if (OnBlockEdge(goal) || IsConvexCorner(goal))
			vGoal.push_back(goal);
	else
		for (unsigned int i = 0; i < bp->visibility_corners.size(); i++)
			if (LineOfSight(goal,bp->visibility_corners[i]))
				vGoal.push_back(bp->visibility_corners[i]);

	// Find the pair of vCorners between start and goal that would yield a shortest path
	cost bestCost = INFINITE_COST;
	int bestStartInd, bestGoalInd;

	for (unsigned int i = 0; i < vStart.size(); i++)
	{
		int startInd = bp->getCornerIndex(vStart[i]);
		cost baseGVal = EuclideanDistance(start, vStart[i]);
		for (unsigned int j = 0; j < vGoal.size(); j++)
		{
			int goalInd = bp->getCornerIndex(vGoal[j]);
			cost newGVal = baseGVal + bp->dist[startInd][goalInd] +  EuclideanDistance(vGoal[j], goal);
			if (newGVal < bestCost)
			{
				bestCost = newGVal;
				bestStartInd = startInd;
				bestGoalInd = goalInd;
			}
		}
	}

	if (!(bestCost < INFINITE_COST))
		return INFINITE_COST;

	// Find the path
	thePath.push_back(corner_locations_[start]);
	while (bestStartInd != bestGoalInd)
	{
		thePath.push_back(corner_locations_[bp->visibility_corners[bestStartInd]]);
		bestStartInd = bp->next[bestStartInd][bestGoalInd];
	}
	thePath.push_back(corner_locations_[bp->visibility_corners[bestStartInd]]);
	thePath.push_back(corner_locations_[goal]);

	return bestCost;
}

cost BlockAStar::FindXYLocPath(xyLoc from, xyLoc to, std::vector<xyLoc> &thePath)
{
	thePath.clear();

#ifdef BLOCK_A_STAR_DEBUG
	printf("Searching from (%d, %d) to (%d, %d)\n", from_.x, from_.y, to_.x, to_.y);
#endif

	cornerId start = corner_ids_[from.x][from.y];
	cornerId goal = corner_ids_[to.x][to.y];

	// If start = goal, return the trivial path
	if (start == goal){
		thePath.push_back(from);
		return 0;
	}

	// If start and goal are in the same block, try to find a path between them
//*
	blockId common = GetCommonBlock(start,goal);
	if (common != -1){
		cost c = FindPathWithinBlock(common, start, goal, thePath);	// TODO: still try to find a better path throught other blocks?
		if (c < INFINITE_COST)
			return c;

		// If no path exists, clear the path for the actual search
		thePath.clear();
	}
//*/
	// Initialize start and goal
	std::vector<cornerId> startBoundaryCorners, goalBoundaryCorners;	// List of boundary corners reachable from the start
	std::vector<cost> startBoundaryCosts, goalBoundaryCosts;
	std::vector<cornerId> startBoundaryNext, goalBoundaryNext;

	if (!OnBlockEdge(goal))
		GetBoundaryCorners(goal, goalBoundaryCorners, goalBoundaryCosts, goalBoundaryNext);

	if (!OnBlockEdge(start))
		GetBoundaryCorners(start, startBoundaryCorners, startBoundaryCosts, startBoundaryNext);

	// Initialize search
	if (search_ >= MAX_SEARCH)
		ResetSearch();

	search_++;
	open_list_.clear();

	GenerateCorner(start,goal);
	corners_[start].g_val = 0;
	corners_[start].parent = start;
	if (OnBlockEdge(start))
	{
		for (unsigned int i = 0; i < corners_[start].belongs_to.size(); i++)
		{
			blockId startBlock = corners_[start].belongs_to[i];
			GenerateBlock(startBlock, goal);
			blocks_[startBlock].g_val = corners_[start].g_val;
			blocks_[startBlock].h_val = corners_[start].h_val;
			AddToOpen(startBlock);
		}
	}
	else
	{
		blockId startBlock = corners_[start].belongs_to[0];

		GenerateBlock(startBlock, goal);
		for (unsigned int i = 0; i < startBoundaryCorners.size(); i++)
		{
			cornerId c = startBoundaryCorners[i];
			corners_[c].g_val = startBoundaryCosts[i];
			corners_[c].parent = startBoundaryNext[i];
			corners_[startBoundaryNext[i]].parent = start;

			if (corners_[c].g_val + corners_[c].h_val + EPSILON < blocks_[startBlock].g_val + blocks_[startBlock].h_val)
			{
				blocks_[startBlock].g_val = corners_[c].g_val;
				blocks_[startBlock].h_val = corners_[c].h_val;
			}
		}
		AddToOpen(startBlock);
	}

	GenerateCorner(goal, goal);
	for (unsigned int i = 0; i < corners_[goal].belongs_to.size(); i++)
	{
		GenerateBlock(corners_[goal].belongs_to[i], goal);
	}
	if (!OnBlockEdge(goal))	// If goal is on a block edge, its g-value will be updated automatically, no special attention needed
		blocks_[corners_[goal].belongs_to[0]].is_goal_block = true;

	// Begin search

	std::vector<unsigned int> ingressIndices;	// Will be used for storing the ingress vertices of the block being expanded
	ingressIndices.reserve(block_cell_width_*4);

	while(!open_list_.empty() && corners_[goal].g_val > GetMin().f_val + EPSILON)
	{
		blockId curr = GetMin().id;
		Block* b = &blocks_[curr];
		PopMin();

		#ifdef ANY_ANGLE_STATISTICS
			num_expansions_++;
		#endif

		// Figure out the ingress vertices
		ingressIndices.clear();
		for (unsigned int i = 0; i < b->num_boundary; i++)
			if (corners_[b->visibility_corners[i]].g_val + EPSILON < b->last_g_val[i])	// If the vertex' g-value is decreased since last expansion
				ingressIndices.push_back(i);

		// Update g-values
		for (unsigned int i = 0; i < ingressIndices.size(); i++)	// Go over the ingress vertices
		{
			unsigned int ind = ingressIndices[i];			// Index of the ingress vertex in vCorners array
			cost gVal = corners_[b->visibility_corners[ind]].g_val;		// g-value of the ingress vertex
			for (unsigned int j = 0; j < b->num_boundary; j++)	// Go over egress vertices
			{
				Corner* egress = &corners_[b->visibility_corners[j]];	// Pointer to the egress vertex
				if (gVal + b->dist[ind][j] + EPSILON < egress->g_val)	// If a shorter path is found to the egress vertex through the ingress vertex
				{
					egress->g_val = gVal + b->dist[ind][j];	// Update g-value
					egress->parent = b->visibility_corners[ind];		// Update parent
				}
			}
		}

		// Look at the updated g-values and generate successors accordingly
		for (unsigned int i = 0; i < b->num_boundary; i++)		// Go over all egress vertices
		{
			Corner* egress = &corners_[b->visibility_corners[i]];	// Pointer to the current egress vertex
			if(egress->g_val + EPSILON < b->last_g_val[i])	// If the egress vertex' g-value has decreased since the last expansion
			{
				b->last_g_val[i] = egress->g_val;	// Set its last expansion g-value to the current one

				for (unsigned int j = 0; j < egress->belongs_to.size(); j++)	// Go over all the blocks that the egress vertex belongs to
				{
					blockId succ = egress->belongs_to[j];	// Successor block
					if (succ != curr)	// Ignore the current block
					{
						GenerateBlock(succ, goal);	// Generate the successor block

						int ind = egress->index_in_block[j];	// index of the egress vertex in the successor block

						if (egress->g_val + EPSILON < blocks_[succ].last_g_val[ind])	// If the g-value of the egress vertex has decreased since the successor block's last expansion
							if(egress->g_val + egress->h_val + EPSILON < blocks_[succ].g_val + blocks_[succ].h_val)	// If the egress vertex' f-value is smaller than the sucessor block's f-value
							{
								blocks_[succ].g_val = egress->g_val;	// Update its g and h values (and therefore, its f-value)
								blocks_[succ].h_val = egress->h_val;
								AddToOpen(succ);	// Add it to the OPEN list, or update its position in the OPEN list
							}
					}
				}
			}
		}

		// If the block contains the goal, update goal's g-value
		if (b->is_goal_block)
		{
			//std::cout<<"Goal check at expansion "<<nExpansions<<std::endl;
			for (unsigned int i = 0; i < goalBoundaryCorners.size(); i++)	// Go over the boundary corners with paths to the goal vertex
			{

				cost newGVal = corners_[goalBoundaryCorners[i]].g_val + goalBoundaryCosts[i];
				//std::cout<<corners[goalBoundaryCorners[i]].gVal<<"\t"<<goalBoundaryCosts[i]<<std::endl;

				if (newGVal + EPSILON < corners_[goal].g_val)
				{
					corners_[goal].g_val = newGVal;
					corners_[goal].parent = goalBoundaryCorners[i];
				}
			}
		}
	}

#ifdef BLOCK_A_STAR_DEBUG
	std::cout<<"Path cost: "<<corners_[goal].g_val<<std::endl;
#endif

	// Search complete, extract path
	if (corners_[goal].g_val < INFINITE_COST)	// Then, we have a path
	{
		cornerId curr = goal;
		cornerId target = corners_[goal].parent;

		// Special case if the goal is in the interior
		if(!OnBlockEdge(goal))
		{
			thePath.push_back(corner_locations_[curr]);

			// Get the visibility corner that goal used to connect to its parent on the boundary corner
			unsigned int i = 0;
			while (i < goalBoundaryCorners.size() && goalBoundaryCorners[i] != target)
				i++;

			curr = goalBoundaryNext[i];	// visibility corner to move to in the goal block
		}

		thePath.push_back(corner_locations_[curr]);

		// Move from parent to parent and lookup the path in-between parents
		while (OnBlockEdge(target) || IsConvexCorner(target))
		{

			#ifdef BLOCK_A_STAR_DEBUG
			//printf("Moving from (%d, %d) to (%d, %d)\n", cornerLocs[curr].x, cornerLocs[curr].y, cornerLocs[target].x, cornerLocs[target].y);
			#endif

			// Choose the best block to traverse
			blockId b;
			cost bestCost = INFINITE_COST;

			for(unsigned int i = 0; i < corners_[curr].belongs_to.size(); i++)
				for(unsigned int j = 0; j < corners_[target].belongs_to.size(); j++)
					if(corners_[curr].belongs_to[i] == corners_[target].belongs_to[j])
					{
						Block* bp = &blocks_[corners_[curr].belongs_to[i]];
						cost newCost = bp->dist[bp->getCornerIndex(curr)][bp->getCornerIndex(target)];
						if (newCost < bestCost)
						{
							bestCost = newCost;
							b = corners_[curr].belongs_to[i];
						}
					}


			// Get the indices for the current and target vertex
			int currInd = blocks_[b].getCornerIndex(curr);
			int targetInd = blocks_[b].getCornerIndex(target);

			while (currInd != targetInd)
			{
				currInd = blocks_[b].next[currInd][targetInd];
				thePath.push_back(corner_locations_[blocks_[b].visibility_corners[currInd]]);
			}

			if (target == start)	// If our target in this iteration is the start vertex, break when we are as close to start as possible
				break;

			curr = target;
			target = corners_[target].parent;
		}

		thePath.push_back(corner_locations_[start]);
		std::reverse(thePath.begin(), thePath.end());
	}

#ifdef BLOCK_A_STAR_DEBUG
	std::cout<<"Path extracted!"<<std::endl;
#endif

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

	return corners_[goal].g_val;
}

void BlockAStar::ResetSearch()
{
	// Set last search and generated values to 0, so that when the search is incremented, all states are un-generated
	search_ = 0;
	for (cornerId i = 0; i < corners_.size(); i++)
		corners_[i].generated = 0;
	for (blockId i = 0; i < blocks_.size(); i++)
		blocks_[i].generated = 0;
}

// Heap operations.
void BlockAStar::AddToOpen(const blockId id)
{
    // If it is already in the open list, update its location.
	if (blocks_[id].list == OPEN_LIST) {
		int i = blocks_[id].heap_index;
		open_list_[i].g_val = blocks_[id].g_val;
		open_list_[i].f_val = blocks_[id].g_val + blocks_[id].h_val;
		PercolateUp(i);
	}
	// Otherwise, add it to the open list.
	else {
		blocks_[id].list = OPEN_LIST;
		blocks_[id].heap_index = open_list_.size();
		open_list_.push_back(HeapElement(id, blocks_[id].g_val, blocks_[id].g_val + blocks_[id].h_val));
		PercolateUp(open_list_.size()-1);
	}
}
void BlockAStar::PopMin()
{
	blocks_[open_list_[0].id].list = CLOSED_LIST;
	blocks_[open_list_[0].id].g_val = INFINITE_COST;
	blocks_[open_list_[0].id].h_val = INFINITE_COST;
	open_list_[0] = open_list_.back();
	blocks_[open_list_[0].id].heap_index = 0;
	open_list_.pop_back();
	PercolateDown(0);
}
void BlockAStar::PercolateUp(int index)
{
	HeapElement elem = open_list_[index];
	int parent = (index-1) >> 1;

	while(index > 0 && elem < open_list_[parent])
	{
		open_list_[index] = open_list_[parent];
		blocks_[open_list_[index].id].heap_index = index;

		index = parent;
		parent = (index-1) >> 1;

#ifdef ANY_ANGLE_STATISTICS
        num_percolations_++;
#endif
	}

	open_list_[index] = elem;
	blocks_[elem.id].heap_index = index;
}
void BlockAStar::PercolateDown(int index)
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
				blocks_[open_list_[index].id].heap_index = index;
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
			blocks_[open_list_[index].id].heap_index = index;
			index = child2;

#ifdef ANY_ANGLE_STATISTICS
			num_percolations_++;
#endif
		}
		else
			break;
	}

	open_list_[index] = elem;
	blocks_[elem.id].heap_index = index;
}

/*
#ifdef ANY_ANGLE_RUNNING_IN_HOG
void BlockAStar::OpenGLDraw(const MapEnvironment *env)
{
	bool displayBlockCorners = 0;
	bool displayBlockBoundaries = 0;
	bool displayConvexCorners = 0;
	bool displayEdges = 0;
	bool displayBlockAStarPath = 0;
	bool displayBlockAStarPSPath = 1;
	bool displayBlockTopLeft = 0;

	env->SetColor(0,0,1);

	if (displayBlockTopLeft)
		for (blockId b = 0; b < blocks.size(); b++)
			env->OpenGLDraw(blocks[b].topLeft);

	if (displayBlockCorners)
		for (cornerId i = 0; i < corner_locations_.size(); i++)
			if (OnBlockCorner(i))
				env->OpenGLDraw(corner_locations_[i]);

	if (displayBlockBoundaries)
		for (cornerId i = 0; i < corner_locations_.size(); i++)
			if (OnBlockEdge(i))
				env->OpenGLDraw(corner_locations_[i]);

	if (displayConvexCorners)
		for (cornerId i = 0; i < corner_locations_.size(); i++)
			if (IsConvexCorner(i))
				env->OpenGLDraw(corner_locations_[i]);

	if (displayBlockAStarPath)
	{
		env->SetColor(0, 140.0/255, 71.0/255);
		DisplayPath(env, blockAStarPath, false);
	}
	if (displayBlockAStarPSPath)
	{
		env->SetColor(0, 140.0/255, 71.0/255);
		DisplayPath(env, blockAStarPSPath, false);
	}
}
#endif
*/
