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

#include "ThetaStar.h"

#ifdef ANY_ANGLE_RUNNING_IN_HOG
ThetaStar::ThetaStar(MapEnvironment *env, int search_method) :
    AnyAngleAlgorithm::AnyAngleAlgorithm(env)
{
    ThetaStarInitialize(search_method);
}
#endif

ThetaStar::ThetaStar(std::vector<bool> &bits, int _width, int _height, int search_method) :
    AnyAngleAlgorithm::AnyAngleAlgorithm(bits, _width, _height)
{
    ThetaStarInitialize(search_method);
}
void ThetaStar::ThetaStarInitialize(const int search_method)
{
    search_method_ = search_method;
    // Generate the edges between cells, for easy look-up.
    corners_.resize(corner_locations_.size());
    for (cornerId i = 0; i < corner_locations_.size(); i++) {
        corners_[i].neighbors = GetNeighbors(i);
        for (unsigned int j = 0; j < corners_[i].neighbors.size(); j++)
            corners_[i].neighbor_dist.push_back(EuclideanDistance(i, corners_[i].neighbors[j]));
    }

    // Reset the generated values.
    ResetSearch();
    open_list_.reserve(10000);

#ifdef ANY_ANGLE_STATISTICS
    // Set up statistics.
    SetupStatistics();
#endif
}
ThetaStar::~ThetaStar()
{
}

const std::string ThetaStar::GetName() const
{
    switch (search_method_) {

    case A_STAR_EUC:
        return "A_EUC";

    case A_STAR_OCT:
        return "A_OCT";

    case THETA_STAR:
        return "T";

    case LAZY_THETA_STAR:
        return "L";

    default:
        return "";
    }
}

cost ThetaStar::FindXYLocPath(xyLoc from, xyLoc to, std::vector<xyLoc> &path)
{
    switch (search_method_) {

    case A_STAR_EUC:
        UseEuclideanDistance();
        return AStarSearch(from, to, path);

    case A_STAR_OCT:
        UseOctileDistance();
        return AStarSearch(from, to, path);

    case THETA_STAR:
        UseEuclideanDistance();
        return ThetaStarSearch(from, to, path);

    case LAZY_THETA_STAR:
        UseEuclideanDistance();
        return LazyThetaStarSearch(from, to, path);

    default:
        return INFINITE_COST;
    }
}

void ThetaStar::ValidateParent(const cornerId s, const cornerId goal)
{
    // Lazy Theta* assumes that there is always line-of-sight from the parent of an expanded state to a successor state.
    // When expanding a state, check if this is true.
    if (!LineOfSight(corners_[s].parent, s)) {

        // Since the previous parent is invalid, set g-value to infinity.
        corners_[s].g_val = INFINITE_COST;

        // Go over potential parents and update its parent to the parent that yields the lowest g-value for s.
        for (unsigned int i = 0; i < corners_[s].neighbors.size(); i++) {
            cornerId new_parent = corners_[s].neighbors[i];
            GenerateState(new_parent, goal);
            if (corners_[new_parent].list == CLOSED_LIST) {
                cost new_g_val = corners_[new_parent].g_val + corners_[s].neighbor_dist[i];
                if (new_g_val < corners_[s].g_val) {
                    corners_[s].g_val = new_g_val;
                    corners_[s].parent = new_parent;
                }
            }
        }
    }
}

cost ThetaStar::AStarSearch(const xyLoc from, const xyLoc to, std::vector<xyLoc> & path)
{
    if (search_ >= MAX_SEARCH)
        ResetSearch();

    // Initialize the search.
    search_++;
    path.clear();

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

        for (unsigned int i = 0; i < corners_[curr].neighbors.size(); i++) {
            cornerId succ = corners_[curr].neighbors[i];
            GenerateState(succ, goal);

            if (corners_[succ].list != CLOSED_LIST) {
                cost new_g_val = corners_[curr].g_val + corners_[curr].neighbor_dist[i];

                if (new_g_val + EPSILON < corners_[succ].g_val) {
                    corners_[succ].g_val = new_g_val;
                    corners_[succ].parent = curr;
                    AddToOpen(succ);
                }
            }
        }
    }

    // Extract the path.
    if (corners_[goal].g_val < INFINITE_COST) {
        cornerId curr = goal;
        while (curr != start) {
            path.push_back(corner_locations_[curr]);
            curr = corners_[curr].parent;
        }

        path.push_back(corner_locations_[curr]);
        std::reverse(path.begin(), path.end());
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

    return corners_[goal].g_val;
}
cost ThetaStar::ThetaStarSearch(const xyLoc from, const xyLoc to, std::vector<xyLoc> & path)
{
    if (search_ >= MAX_SEARCH)
        ResetSearch();

    // Initialize the search.
    search_++;
    path.clear();

    cornerId start = corner_ids_[from.x][from.y];
    cornerId goal = corner_ids_[to.x][to.y];

    open_list_.clear();
    GenerateState(start, goal);
    GenerateState(goal, goal);
    corners_[start].g_val = 0;

    // Set 'start's parent as itself. When 'start' is expanded for the first time, the grandparent of 'start's successor will be 'start' as well.
    corners_[start].parent = start;

    AddToOpen(start);

    while (!open_list_.empty() && corners_[goal].g_val > GetMin().f_val + EPSILON) {
        cornerId curr = GetMin().id;
        PopMin();

#ifdef ANY_ANGLE_STATISTICS
        num_expansions_++;
#endif

        for (unsigned int i = 0; i < corners_[curr].neighbors.size(); i++) {
            cornerId succ = corners_[curr].neighbors[i];
            GenerateState(succ, goal);

            if (corners_[succ].list != CLOSED_LIST) {
                cost new_g_val;
                cornerId new_parent;

                // First, check if there is LOS from the grandparent.
                if (LineOfSight(corners_[curr].parent, succ)) {
                    new_parent = corners_[curr].parent;
                    new_g_val = corners_[new_parent].g_val + EuclideanDistance(new_parent, succ);
                }
                else {
                    new_parent = curr;
                    new_g_val = corners_[new_parent].g_val + corners_[curr].neighbor_dist[i];
                }

                if (new_g_val + EPSILON < corners_[succ].g_val) {
                    corners_[succ].g_val = new_g_val;
                    corners_[succ].parent = new_parent;
                    AddToOpen(succ);
                }
            }
        }
    }

    // Extract the path.
    if (corners_[goal].g_val < INFINITE_COST)
    {
        cornerId curr = goal;
        while (curr != start) {
            path.push_back(corner_locations_[curr]);
            curr = corners_[curr].parent;
        }

        path.push_back(corner_locations_[curr]);
        std::reverse(path.begin(), path.end());
    }

#ifdef ANY_ANGLE_ASSERTIONS
    // Make sure that the heap property is satisfied for the open list
    for (int i = 0; i < open_list_.size(); i++) {
        int child1 = (i << 1) + 1;
        int child2 = (i << 1) + 2;

        if (child2 < open_list_.size()) {
            // if (!(open_list_[i].f_val < open_list_[child1].f_val + EPSILON) || !(open_list_[i].f_val < open_list_[child2].f_val + EPSILON)) {
            //     std::cout<<std::endl<<i<<":\t"<<open_list_[i].f_val<<"\t"<<open_list_[child1].f_val<<"\t"<<open_list_[child2].f_val<<std::endl;
            // }
            assert(open_list_[i].f_val < open_list_[child1].f_val + 100*EPSILON);
            assert(open_list_[i].f_val < open_list_[child2].f_val + 100*EPSILON);
        }
    }
#endif

    return corners_[goal].g_val;
}
cost ThetaStar::LazyThetaStarSearch(const xyLoc from, const xyLoc to, std::vector<xyLoc> & path)
{
    if (search_ >= MAX_SEARCH)
        ResetSearch();

    // Initialize the search.
    search_++;
    path.clear();

    cornerId start = corner_ids_[from.x][from.y];
    cornerId goal = corner_ids_[to.x][to.y];

    open_list_.clear();
    GenerateState(start, goal);
    GenerateState(goal, goal);
    corners_[start].g_val = 0;

    // Set 'start's parent as itself. When 'start' is expanded for the first time, the grandparent of 'start's successor will be 'start' as well.
    corners_[start].parent = start;

    AddToOpen(start);

    while (!open_list_.empty() && corners_[goal].g_val > GetMin().f_val + EPSILON) {
        cornerId curr = GetMin().id;
        PopMin();

#ifdef ANY_ANGLE_STATISTICS
        num_expansions_++;
#endif

        // First, check if we have LOS from the parent and update the parent and g-value of curr if necessary.
        ValidateParent(curr, goal);

        // Assume there is LOS from the parent of the state being expanded to the successor
        cornerId new_parent = corners_[curr].parent;

        for (unsigned int i = 0; i < corners_[curr].neighbors.size(); i++) {
            cornerId succ = corners_[curr].neighbors[i];
            GenerateState(succ, goal);

            if (corners_[succ].list != CLOSED_LIST) {
                cost new_g_val = corners_[new_parent].g_val + EuclideanDistance(new_parent, succ);
                if (new_g_val + EPSILON < corners_[succ].g_val) {
                    corners_[succ].g_val = new_g_val;
                    corners_[succ].parent = new_parent;
                    AddToOpen(succ);
                }
            }
        }
    }

    // Extract the path.
    if (corners_[goal].g_val < INFINITE_COST)
    {
        ValidateParent(goal, goal);
        cornerId curr = goal;
        while (curr != start) {
            path.push_back(corner_locations_[curr]);
            curr = corners_[curr].parent;
        }

        path.push_back(corner_locations_[curr]);
        std::reverse(path.begin(), path.end());
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

    return corners_[goal].g_val;
}

void ThetaStar::GenerateState(cornerId s, cornerId goal)
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
void ThetaStar::ResetSearch()
{
    // Set last search and generated values to 0, so that when the search is incremented, all states are un-generated
    search_ = 0;
    for (cornerId i = 0; i < corners_.size(); i++)
        corners_[i].generated = 0;
}

// Heap operations.
void ThetaStar::AddToOpen(cornerId id)
{
    // If it is already in the open list, update its location.
    if (corners_[id].list == OPEN_LIST) {
        int i = corners_[id].heap_index;
        open_list_[i].g_val = corners_[id].g_val;
        open_list_[i].f_val = corners_[id].g_val + corners_[id].h_val;
        PercolateUp(i);
    }
    // Otherwise, add it to the open list
    else {
        corners_[id].list = OPEN_LIST;
        corners_[id].heap_index = open_list_.size();
        open_list_.push_back(HeapElement(id, corners_[id].g_val, corners_[id].g_val + corners_[id].h_val));
        PercolateUp(open_list_.size() - 1);
    }
}
void ThetaStar::PopMin()
{
    corners_[open_list_[0].id].list = CLOSED_LIST;
    open_list_[0] = open_list_.back();
    corners_[open_list_[0].id].heap_index = 0;
    open_list_.pop_back();
    PercolateDown(0);
}
void ThetaStar::PercolateUp(int index)
{
    HeapElement elem = open_list_[index];
    int parent = (index - 1) >> 1;

    while (index > 0 && elem < open_list_[parent]) {
        open_list_[index] = open_list_[parent];
        corners_[open_list_[index].id].heap_index = index;

#ifdef ANY_ANGLE_STATISTICS
        num_percolations_++;
#endif

        index = parent;
        parent = (index - 1) >> 1;
    }

    open_list_[index] = elem;
    corners_[elem.id].heap_index = index;
}
void ThetaStar::PercolateDown(int index)
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

/*
 #ifdef ANY_ANGLE_RUNNING_IN_HOG
 void ThetaStar::OpenGLDraw(const MapEnvironment *env)
 {
 bool displayCorners = 0;
 bool displayEdges = 0;
 bool displaySearchTree = 0;
 bool displayAStarPath = 0;
 bool displayAStarOctPath = 1;
 bool displayAStarPSPath = 0;
 bool displayAStarOctPSPath = 0;
 bool displayThetaStarPath = 0;
 bool displayThetaStarPSPath = 0;
 bool displayLazyThetaStarPath = 0;

 env->SetColor(0,0,1);

 for (cornerId i = 0; i < corners.size(); i++) //  && i < 1000
 {
 if (displayCorners)
 env->OpenGLDraw(corner_locations_[i]);

 for (unsigned int j = 0; j < corners[i].neighbors.size(); j++)
 {
 cornerId i2 = corners[i].neighbors[j];
 if (i > i2)	// Do not draw the same edge twice, from different corners
 if (displayEdges)
 env->GLDrawLine(corner_locations_[i], corner_locations_[i2]);
 }

 if (displaySearchTree)
 {
 if (search_ != 0 && corners[i].generated == search_)
 {
 cornerId parent = corners[i].parent;
 DrawLine(env, corner_locations_[i], corner_locations_[parent]);
 }
 }
 }

 if (displayAStarPath)
 {
 env->SetColor(0, 0, 0);
 DisplayPath(env, aStarPath, 0);
 }
 if (displayAStarOctPath)
 {
 env->SetColor(0, 0, 0);
 DisplayPath(env, aStarOctPath, 0);
 }
 if (displayAStarPSPath)
 {
 env->SetColor(0, 1, 0);
 DisplayPath(env, aStarPSPath,0);
 }

 if (displayAStarOctPSPath)
 {
 env->SetColor(179.0/255, 56.0/255, 147.0/255);
 DisplayPath(env, aStarOctPSPath,0);
 }

 if (displayThetaStarPath)
 {
 env->SetColor(THETA_STAR_COLOR);
 DisplayPath(env, thetaStarPath, 0);
 }

 if (displayThetaStarPSPath)
 {
 env->SetColor(THETA_STAR_COLOR);
 DisplayPath(env, thetaStarPSPath, 0);
 }

 if (displayLazyThetaStarPath)
 {
 env->SetColor(1, 0, 0);
 DisplayPath(env, lazyThetaStarPath, 0);
 }
 }
 #endif
 */

