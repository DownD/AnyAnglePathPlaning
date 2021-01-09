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

#include "FieldAStar.h"

#ifdef ANY_ANGLE_RUNNING_IN_HOG
FieldAStar::FieldAStar(MapEnvironment *env) :
    AnyAngleAlgorithm::AnyAngleAlgorithm(env)
{
    FieldAStarInitialize();
}
#endif

FieldAStar::FieldAStar(std::vector<bool> &bits, int _width, int _height) :
    AnyAngleAlgorithm::AnyAngleAlgorithm(bits, _width, _height)
{
    FieldAStarInitialize();
}
void FieldAStar::FieldAStarInitialize()
{
    // Generate the edges between cells, for easy look-up
    corners_.resize(corner_locations_.size());

    /*
     for (cornerId i = 0; i < cornerLocs.size(); i++)
     {
     corners[i].neighbors = GetNeighbors(i);
     for (unsigned int j = 0; j < corners[i].neighbors.size(); j++)
     corners[i].neighborDist.push_back(EuclideanDistance(i, corners[i].neighbors[j]));
     }
     /*/

    for (cornerId c = 0; c < corner_locations_.size(); c++) {
        xyLoc l = corner_locations_[c];
        corners_[c].neighbors.clear();

        // North
        if (IsTraversable(NorthWestCell(l)) || IsTraversable(NorthEastCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x][l.y - 1]);
        else
            corners_[c].neighbors.push_back(-1);

        // NorthEast
        if (IsTraversable(NorthEastCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x + 1][l.y - 1]);
        else
            corners_[c].neighbors.push_back(-1);

        // East
        if (IsTraversable(NorthEastCell(l)) || IsTraversable(SouthEastCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x + 1][l.y]);
        else
            corners_[c].neighbors.push_back(-1);

        // SouthEast
        if (IsTraversable(SouthEastCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x + 1][l.y + 1]);
        else
            corners_[c].neighbors.push_back(-1);

        // South
        if (IsTraversable(SouthWestCell(l)) || IsTraversable(SouthEastCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x][l.y + 1]);
        else
            corners_[c].neighbors.push_back(-1);

        // SouthWest
        if (IsTraversable(SouthWestCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x - 1][l.y + 1]);
        else
            corners_[c].neighbors.push_back(-1);

        // West
        if (IsTraversable(NorthWestCell(l)) || IsTraversable(SouthWestCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x - 1][l.y]);
        else
            corners_[c].neighbors.push_back(-1);

        // NorthWest
        if (IsTraversable(NorthWestCell(l)))
            corners_[c].neighbors.push_back(corner_ids_[l.x - 1][l.y - 1]);
        else
            corners_[c].neighbors.push_back(-1);

        for (unsigned int i = 0; i < corners_[c].neighbors.size(); i++)
            if (corners_[c].neighbors[i] == -1)
                corners_[c].neighbor_dist.push_back(INFINITE_COST);
            else
                corners_[c].neighbor_dist.push_back(EuclideanDistance(c, corners_[c].neighbors[i]));
    }
    //*/

    // Reset the generated values
    ResetSearch();
    open_list_.reserve(10000);

    // Set up statistics
#ifdef ANY_ANGLE_STATISTICS
    SetupStatistics();
#endif
}
FieldAStar::~FieldAStar()
{
}

void FieldAStar::ResetSearch()
{
    // Set last search and generated values to 0, so that when the search is incremented, all states are un-generated
    search_ = 0;
    for (cornerId i = 0; i < corners_.size(); i++)
        corners_[i].generated = 0;
}

void FieldAStar::GenerateState(cornerId s, cornerId goal)
{
    if (corners_[s].generated < search_) {
#ifdef ANY_ANGLE_STATISTICS
        num_generated_++;
#endif
        corners_[s].generated = search_;
        corners_[s].h_val = HeuristicDistance(s, goal);
        corners_[s].g_val = INFINITE_COST;
        corners_[s].list = NO_LIST;
    }
}
void FieldAStar::InterpolateContinuous(double a, double b, double cCost, double dCost, double & y, double & c)
{
    double f = (cCost - dCost) / b;

    if (f <= EPSILON) // Cardinal neighbor has a smaller g-value, no need to interpolate
    {
        y = 0;
        c = cCost + a;
    }

    else if (f > 1 - EPSILON) // Diagonal neighbor is much more preferable, use it
    {
        y = b;
        c = dCost + sqrt(a * a + b * b);
    }
    else {
        // Now, we have 0 < f < 1 and the best point lies between l1 and l2, so we need to interpolate
        y = a * f / sqrt(1 - f * f);
        y = y < b ? y : b;
        c = sqrt(a * a + y * y) - f * y + cCost;
    }

    //printf("a: %g, b: %g, cCost: %g, dCost: %g f: %g, y: %g, c: %g \n", a, b, cCost, dCost, f, y, c);
}
cost FieldAStar::ComputeCost(cornerId s, cornerId sc, cornerId sd)
{

    // We know that sd is a valid neighbor of s (from the search function), so the cell in question must be unblocked

    double f = corners_[sc].g_val - corners_[sd].g_val;

    //std::cout<<f<<std::endl;

    if (f <= EPSILON) // Cardinal neighbor has a smaller g-value, no need to interpolate
    {
        return 1 + corners_[sc].g_val;
    }

    if (f > 1 - EPSILON) // Diagonal neighbor is much more preferable, use it
    {
        return DIAG_COST + corners_[sd].g_val;
    }
    // OPTIMIZATION: If the cell symmetric to the current cell wrt the lc-ld line is blocked, no need to interpolate?

    // Now, we have 0 < f < 1 and the best point lies between l1 and l2, so we need to interpolate
    double y = f / sqrt(1 - f * f);
    y = y < 1 ? y : 1;

    return sqrt(1 + y * y) + f * (1 - y) + corners_[sd].g_val;
}

bool IsInteger(float z)
{
    float d = z - (int) (z + EPSILON);
    return (d < EPSILON) && (d > -EPSILON);
}

void FieldAStar::GetNeighbors(xyLocCont l, std::vector<xyLocCont> & neighbors, std::vector<cost> & gValues, std::vector<bool> & exists, int lookahead)
{
    neighbors.clear();
    gValues.clear();
    exists.clear();

    float xLocs[3] = { ceil(l.x - 1 - EPSILON), l.x, floor(l.x + 1 + EPSILON) };
    float yLocs[3] = { ceil(l.y - 1 - EPSILON), l.y, floor(l.y + 1 + EPSILON) };

    // Fill the neighbors array with the expected locations of neighbors
    neighbors.push_back(xyLocCont(xLocs[1], yLocs[0])); // North
    neighbors.push_back(xyLocCont(xLocs[2], yLocs[0])); // North East
    neighbors.push_back(xyLocCont(xLocs[2], yLocs[1])); // Clockwise
    neighbors.push_back(xyLocCont(xLocs[2], yLocs[2]));
    neighbors.push_back(xyLocCont(xLocs[1], yLocs[2]));
    neighbors.push_back(xyLocCont(xLocs[0], yLocs[2]));
    neighbors.push_back(xyLocCont(xLocs[0], yLocs[1]));
    neighbors.push_back(xyLocCont(xLocs[0], yLocs[0]));

    // Check if the neighbors are valid
    bool nw = IsTraversable(SouthEastCell(xyLocCont(xLocs[0], yLocs[0]).getCorner())); // Is the southeast cell of the northwest neighbor is traversable?
    bool ne = IsTraversable(SouthWestCell(xyLocCont(xLocs[2], yLocs[0]).getCorner()));
    bool sw = IsTraversable(NorthEastCell(xyLocCont(xLocs[0], yLocs[2]).getCorner()));
    bool se = IsTraversable(NorthWestCell(xyLocCont(xLocs[2], yLocs[2]).getCorner()));

    // 0 means it is a neighbor, INFINITE_COST means it is not
    exists.push_back(nw || ne); // North
    exists.push_back(ne); // North East
    exists.push_back(ne || se); // Clockwise
    exists.push_back(se);
    exists.push_back(se || sw);
    exists.push_back(sw);
    exists.push_back(sw || nw);
    exists.push_back(nw);

    // Go over the neighbors and lookup the g-values of those that are actual corners and not points on an axis
    for (int d = 0; d < 8; d++)
        if (exists[d] && IsInteger(neighbors[d].x) && IsInteger(neighbors[d].y)) {
            xyLoc l1 = neighbors[d].getCorner();
            cornerId c = corner_ids_[l1.x][l1.y];
            GenerateState(c, corner_ids_[to_.x][to_.y]);

            if (lookahead > 0) {
                cost gVal;
                xyLocCont dummy;
                GetBestCostAndParent(neighbors[d], dummy, gVal, lookahead - 1);

                if (gVal < corners_[c].g_val)
                    corners_[c].g_val = gVal;
            }

            gValues.push_back(corners_[c].g_val);
        }
        else
            gValues.push_back(INFINITE_COST);

    // If there are non-corner places, we need to interpolate
    for (int d = 0; d < 8; d++)
        if (exists[d] && !(IsInteger(neighbors[d].x) && IsInteger(neighbors[d].y))) {
            // d's cw and ccw neighbors
            int cw = (d + 1) & 7;
            int ccw = (d + 7) & 7;

            double y;

            if (d == 0)
                y = neighbors[1].x - neighbors[0].x;
            else if (d == 4)
                y = neighbors[4].x - neighbors[5].x;
            else if (d == 2)
                y = neighbors[3].y - neighbors[2].y;
            else
                // d == 6
                y = neighbors[6].y - neighbors[7].y;

            gValues[d] = y * gValues[ccw] + (1 - y) * gValues[cw];
        }

    /*
     printf("\n");
     printf("l: (%g, %g)\n", l.x, l.y);

     for (int d = 0; d < 8; d++)
     {
     if (exists[d])
     {
     printf("(%g, %g)\t%g\n", neighbors[d].x, neighbors[d].y, gValues[d]);
     }
     else
     printf("(%g, %g)\t-\n", neighbors[d].x, neighbors[d].y);
     }
     //*/
}

void FieldAStar::GetBestCostAndParent(xyLocCont l, xyLocCont & neighbor, cost & gVal, int lookahead)
{
    std::vector<xyLocCont> neighbors;
    std::vector<cost> gValues;
    std::vector<bool> exists;

    GetNeighbors(l, neighbors, gValues, exists, lookahead);

    gVal = INFINITE_COST;

    for (unsigned int d = 1; d < neighbors.size(); d += 2) // Iterate over diagonal neighbors
    {
        if (exists[d]) {
            int cw = (d + 1) & 7; // clockwise neighbor of curr, among succ's neighbors
            int ccw = (d + 7) & 7; // counter-clockwise neighbor of curr, among succ's neighbors

            xyLocCont dLoc = neighbors[d];
            xyLocCont cwLoc = neighbors[cw];
            xyLocCont ccwLoc = neighbors[ccw];

            double dx = dLoc.x - l.x;
            double dy = dLoc.y - l.y;

            double ycw, yccw, a, b;

            // cw neighbor
            if (d == 1 || d == 5) // cw direction in the y axis
            {
                b = dy > 0 ? dy : -dy; // interpolate along the y axis
                a = dx > 0 ? dx : -dx;
            }
            else {
                b = dx > 0 ? dx : -dx; // interpolate along the x axis
                a = dy > 0 ? dy : -dy;
            }

            double cwCost;
            double ccwCost;

            InterpolateContinuous(a, b, gValues[cw], gValues[d], ycw, cwCost);

            // for ccw neighbor, the directions are different
            InterpolateContinuous(b, a, gValues[ccw], gValues[d], yccw, ccwCost);

            if (d == 1) // Diagonal is at NorthEast
            {
                cwLoc = xyLocCont(cwLoc.x, cwLoc.y - ycw); // Clockwise is to the East, therefore y moves towards North
                ccwLoc = xyLocCont(ccwLoc.x + yccw, ccwLoc.y); // CounterClockwise is to the North, therefore, y moves towards East
            }
            else if (d == 3) {
                cwLoc = xyLocCont(cwLoc.x + ycw, cwLoc.y);
                ccwLoc = xyLocCont(ccwLoc.x, ccwLoc.y + yccw);
            }
            else if (d == 5) {
                cwLoc = xyLocCont(cwLoc.x, cwLoc.y + ycw);
                ccwLoc = xyLocCont(ccwLoc.x - yccw, ccwLoc.y);
            }
            else //if (d == 7)
            {
                cwLoc = xyLocCont(cwLoc.x - ycw, cwLoc.y);
                ccwLoc = xyLocCont(ccwLoc.x, ccwLoc.y - yccw);
            }

            if (cwCost < gVal) {
                gVal = cwCost;
                neighbor = cwLoc;
            }

            if (ccwCost < gVal) {
                gVal = ccwCost;
                neighbor = ccwLoc;
            }
        }
    }

}

cost FieldAStar::ExtractPath(xyLoc startLoc, xyLoc goalLoc, std::vector<xyLocCont> & thePath)
{
    thePath.clear();
    cost pathCost = 0;

    xyLocCont curr = xyLocCont(goalLoc);
    xyLocCont start = xyLocCont(startLoc);
    thePath.push_back(curr);

    while (!(curr == start)) {
        //std::cout<<curr.x<<" "<<curr.y<<std::endl;
        cost gVal;
        xyLocCont next;
        this->GetBestCostAndParent(curr, next, gVal, FIELD_A_STAR_LOOKAHEAD);
        thePath.push_back(next);
        pathCost += EuclideanDistance(curr, next);
        curr = next;
    }
    std::reverse(thePath.begin(), thePath.end());

    //std::cout<<"Path cost: "<<pathCost<<std::endl;
    return pathCost;
}

cost FieldAStar::FindXYLocContPath(xyLoc from, xyLoc to, std::vector<xyLocCont> &thePath)
{
    if (search_ >= MAX_SEARCH)
        ResetSearch();

    search_++;
    thePath.clear();

    cornerId start = corner_ids_[from.x][from.y];
    cornerId goal = corner_ids_[to.x][to.y];

    open_list_.clear();
    GenerateState(start, goal);
    GenerateState(goal, goal);

    corners_[start].g_val = 0;
    AddToOpen(start);

    while (!open_list_.empty() && corners_[goal].g_val > GetMin().f_val + EPSILON) {
        cornerId curr = GetMin().id;
        PopMin();

#ifdef ANY_ANGLE_STATISTICS
        num_expansions_++;
#endif

        for (unsigned int d = 0; d < corners_[curr].neighbors.size(); d++) {
            cornerId succ = corners_[curr].neighbors[d];
            if (succ != -1) {
                GenerateState(succ, goal);

                if (corners_[succ].list != CLOSED_LIST) {
                    unsigned int revD = (d + 4) & 7; // reversed direction of d, & 7 = % 8 = mod 8
                    cornerId cw = corners_[succ].neighbors[(revD + 1) & 7]; // clockwise neighbor of curr, among succ's neighbors
                    cornerId ccw = corners_[succ].neighbors[(revD + 7) & 7]; // counter-clockwise neighbor of curr, among succ's neighbors

                    cost cwCost = INFINITE_COST;
                    cost ccwCost = INFINITE_COST;

                    if (cw != -1) // If clockwise neighbor exists
                    {
                        GenerateState(cw, goal);
                        if (d & 1) // If d (and therefore revD) is a diagonal direction (direction 0 = North, then clockwise with 7 = NorthWest, so (d mod 2) == 1 is a diagonal direction)
                            cwCost = ComputeCost(succ, cw, curr); // curr is a diagonal neighbor of succ and cw is a cardinal neighbor

                        else
                            cwCost = ComputeCost(succ, curr, cw); // otherwise, curr is a cardinal neighbor of succ and cw is a diagona neighbor
                    }

                    if (ccw != -1) // If clockwise neighbor exists
                    {
                        GenerateState(ccw, goal);
                        if (d & 1)
                            ccwCost = ComputeCost(succ, ccw, curr);
                        else
                            ccwCost = ComputeCost(succ, curr, ccw);
                    }

                    bool updated = false;
                    if (cwCost < corners_[succ].g_val - EPSILON) {
                        corners_[succ].g_val = cwCost;
                        corners_[succ].parent = curr; // curr's clockwise neighbor in succ's neighbor list is cw
                        updated = true;
                    }

                    if (ccwCost < corners_[succ].g_val - EPSILON) {
                        corners_[succ].g_val = ccwCost;
                        corners_[succ].parent = ccw; // ccw's clockwise neighbor in succ's neighbor list is ccw
                        updated = true;
                    }

                    if (updated)
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

    if (!(corners_[goal].g_val < INFINITE_COST))
        return INFINITE_COST;

    return ExtractPath(from, to, thePath);
}

// Heap operations.
void FieldAStar::AddToOpen(const cornerId id)
{
    // If it is already in the open list, update its location.
    if (corners_[id].list == OPEN_LIST)
    {
        int i = corners_[id].heap_index;
        open_list_[i].g_val = corners_[id].g_val;
        open_list_[i].f_val = corners_[id].g_val + corners_[id].h_val;
        PercolateUp(i);
    }
    // Otherwise, add it to the open list.
    else
    {
        corners_[id].list = OPEN_LIST;
        corners_[id].heap_index = open_list_.size();
        open_list_.push_back(HeapElement(id, corners_[id].g_val, corners_[id].g_val + corners_[id].h_val));
        PercolateUp(open_list_.size() - 1);
    }
}
void FieldAStar::PopMin()
{
    corners_[open_list_[0].id].list = CLOSED_LIST;
    open_list_[0] = open_list_.back();
    corners_[open_list_[0].id].heap_index = 0;
    open_list_.pop_back();
    PercolateDown(0);
}
void FieldAStar::PercolateUp(int index)
{
    HeapElement elem = open_list_[index];
    int parent = (index - 1) >> 1;

    while (index > 0 && elem < open_list_[parent]) {
        open_list_[index] = open_list_[parent];
        corners_[open_list_[index].id].heap_index = index;

        index = parent;
        parent = (index - 1) >> 1;

#ifdef ANY_ANGLE_STATISTICS
        num_percolations_++;
#endif
    }

    open_list_[index] = elem;
    corners_[elem.id].heap_index = index;
}
void FieldAStar::PercolateDown(int index)
{
    HeapElement elem = open_list_[index];
    int maxSize = open_list_.size();

    while (true) {
        int child1 = (index << 1) + 1;
        if (child1 >= maxSize)
            break;

        int child2 = child1 + 1;

        // If the first child has smaller key
        if (child2 == maxSize || open_list_[child1] < open_list_[child2]) {
            if (open_list_[child1] < elem) {
                open_list_[index] = open_list_[child1];
                corners_[open_list_[index].id].heap_index = index;
                index = child1;

#ifdef ANY_ANGLE_STATISTICS
                num_percolations_++;
#endif
            }
            else
                break;
        }

        else if (open_list_[child2] < elem) {
            open_list_[index] = open_list_[child2];
            corners_[open_list_[index].id].heap_index = index;
            index = child2;

#ifdef ANY_ANGLE_STATISTICS
            num_percolations_++;
#endif
        }
        else
            break;
    }

    open_list_[index] = elem;
    corners_[elem.id].heap_index = index;
}

