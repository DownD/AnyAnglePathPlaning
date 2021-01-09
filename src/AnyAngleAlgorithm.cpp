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

#include "AnyAngleAlgorithm.h"

#ifdef ANY_ANGLE_RUNNING_IN_HOG
AnyAngleAlgorithm::AnyAngleAlgorithm(const MapEnvironment *env)
{
    // Read the height and width of the map.
    Map* map = env->GetMap();
    width_ = map->GetMapWidth() + 2; // Add +2 for the frame of obstacles.
    height_ = map->GetMapHeight() + 2;

    // Generate the cells with the padding.
    cell_traversable_ = new bool*[width_];
    for (unsigned int x = 0; x < width_; x++) {
        cell_traversable_[x] = new bool[height_];
        for (unsigned int y = 0; y < height_; y++)
            cell_traversable_[x][y] = false;    // Block all cells (including the padded cells.
    }

    // Make unblocked cells traversable.
    int num_traversable = 0;
    for (unsigned int x = 0; x < width_ - 2; x++) {
        for (unsigned int y = 0; y < height_ - 2; y++) {
            if (map->GetTerrainType(x, y) == kGround) {
                cell_traversable_[x + 1][y + 1] = true;
                num_traversable++;
            }
        }
    }

    debug_loc_ = xyLoc(1,1);

#ifdef ANY_ANGLE_VERBOSE
    std::cout << "Traversable cells: " << num_traversable << std::endl;
#endif

    Initialize();
}
#endif
AnyAngleAlgorithm::AnyAngleAlgorithm(const std::vector<bool> &bits, int width, int height)
{
    // Read the height and width of the map.
    width_ = width + 2; // Add +2 for the frame of obstacles.
    height_ = height + 2;

    // Generate the cells with the padding.
    cell_traversable_ = new bool*[width_];
    for (unsigned int x = 0; x < width_; x++) {
        cell_traversable_[x] = new bool[height_];
        for (unsigned int y = 0; y < height_; y++)
            cell_traversable_[x][y] = false;    // Block all cells (including the padded cells.
    }

    // Make unblocked cells traversable.
    int num_traversable = 0;
    for (unsigned int x = 0; x < width_ - 2; x++) {
        for (unsigned int y = 0; y < height_ - 2; y++) {
            if (bits[y * width + x]) {
                cell_traversable_[x + 1][y + 1] = true;
                num_traversable++;
            }
        }
    }

#ifdef ANY_ANGLE_VERBOSE
    std::cout << "Traversable cells: " << num_traversable << std::endl;
#endif

    Initialize();
}
AnyAngleAlgorithm::~AnyAngleAlgorithm()
{
    for (unsigned int x = 0; x < width_ - 1; x++) {
        delete[] cell_traversable_[x];
        delete[] corner_ids_[x];
    }
    delete[] cell_traversable_[width_ - 1];

    delete[] cell_traversable_;
    delete[] corner_ids_;

#ifdef ANY_ANGLE_STATISTICS
    if (statistics_ != NULL)
        delete statistics_;
    if (statistics_with_smoothing_ != NULL)
        delete statistics_with_smoothing_;
#endif
}
void AnyAngleAlgorithm::Initialize()
{
    // Generate 'corner_ids_' with the padding.
    corner_ids_ = new cornerId*[width_ - 1];
    for (unsigned int x = 0; x < width_ - 1; x++) {
        corner_ids_[x] = new cornerId[height_ - 1];
    }

    // Determine the corners and create the two-way mapping between cornerIds and corner locations.
    corner_locations_.clear();
    valid_corner_locations_.clear();

    for (unsigned int y = 0; y < height_ - 1; y++) {    // Ignore the last row and the last column
        for (unsigned int x = 0; x < width_ - 1; x++) {
            xyLoc loc = xyLoc(x, y); // Location of the corner, not the cell

            corner_ids_[x][y] = corner_locations_.size();
            corner_locations_.push_back(loc);

            // A corner is a 'valid' corner if one of its four surrounding cells is traversable
            if (IsTraversable(NorthEastCell(loc)) || IsTraversable(NorthWestCell(loc)) ||
                IsTraversable(SouthEastCell(loc)) || IsTraversable(SouthWestCell(loc))) {
                valid_corner_locations_.push_back(loc);
            }
        }
    }

    using_octile_distance_ = false;

#ifdef ANY_ANGLE_VERBOSE
    printf("Corners: %u\n", corner_locations_.size());
    printf("Valid corners: %u\n", valid_corner_locations_.size());
#endif

#ifdef ANY_ANGLE_STATISTICS
    statistics_ = NULL;
    statistics_with_smoothing_ = NULL;
#endif
}

cost AnyAngleAlgorithm::FindPath(xyLoc from, xyLoc to)
{

#ifdef ANY_ANGLE_STATISTICS
    StartStatistics(from, to);
#endif

    cost c;
    if (!UsingXYLocCont()) {
        xyloc_path_.clear();
        c = FindXYLocPath(from, to, xyloc_path_);
#ifdef ANY_ANGLE_STATISTICS
        ReportStatistics(xyloc_path_, statistics_);
#endif

        if (ShouldSmoothPaths()) {
            c = SmoothPath(xyloc_path_, smoothed_xyloc_path_);
#ifdef ANY_ANGLE_STATISTICS
            ReportStatistics(smoothed_xyloc_path_, statistics_with_smoothing_);
#endif
        }
    }
    else {
        xyloc_cont_path_.clear();
        c = FindXYLocContPath(from, to, xyloc_cont_path_);
#ifdef ANY_ANGLE_STATISTICS
        ReportStatistics(xyloc_cont_path_, statistics_);
#endif
        if (ShouldSmoothPaths()) {
            c = SmoothPath(xyloc_cont_path_, smoothed_xyloc_cont_path_);
#ifdef ANY_ANGLE_STATISTICS
            ReportStatistics(smoothed_xyloc_cont_path_, statistics_with_smoothing_);
#endif
        }
    }

    return c;
}

std::vector<cornerId> AnyAngleAlgorithm::GetNeighbors(const cornerId c) const
{
    xyLoc l = corner_locations_[c];
    std::vector<cornerId> neighbors;

    // NorthWest
    if (IsTraversable(NorthWestCell(l)))
        neighbors.push_back(corner_ids_[l.x - 1][l.y - 1]);

    // NorthEast
    if (IsTraversable(NorthEastCell(l)))
        neighbors.push_back(corner_ids_[l.x + 1][l.y - 1]);

    // SouthWest
    if (IsTraversable(SouthWestCell(l)))
        neighbors.push_back(corner_ids_[l.x - 1][l.y + 1]);

    // SouthEast
    if (IsTraversable(SouthEastCell(l)))
        neighbors.push_back(corner_ids_[l.x + 1][l.y + 1]);

    // North
    if (IsTraversable(NorthWestCell(l)) || IsTraversable(NorthEastCell(l)))
        neighbors.push_back(corner_ids_[l.x][l.y - 1]);

    // South
    if (IsTraversable(SouthWestCell(l)) || IsTraversable(SouthEastCell(l)))
        neighbors.push_back(corner_ids_[l.x][l.y + 1]);

    // West
    if (IsTraversable(NorthWestCell(l)) || IsTraversable(SouthWestCell(l)))
        neighbors.push_back(corner_ids_[l.x - 1][l.y]);

    // East
    if (IsTraversable(NorthEastCell(l)) || IsTraversable(SouthEastCell(l)))
        neighbors.push_back(corner_ids_[l.x + 1][l.y]);

    return neighbors;
}
std::vector<xyLoc> AnyAngleAlgorithm::GetNeighbors(const xyLoc l) const
{
    std::vector<xyLoc> neighbors;

    // NorthWest
    if (IsTraversable(NorthWestCell(l)))
        neighbors.push_back(xyLoc(l.x - 1, l.y - 1));

    // NorthEast
    if (IsTraversable(NorthEastCell(l)))
        neighbors.push_back(xyLoc(l.x + 1, l.y - 1));

    // SouthWest
    if (IsTraversable(SouthWestCell(l)))
        neighbors.push_back(xyLoc(l.x - 1, l.y + 1));

    // SouthEast
    if (IsTraversable(SouthEastCell(l)))
        neighbors.push_back(xyLoc(l.x + 1, l.y + 1));

    // North
    if (IsTraversable(NorthWestCell(l)) || IsTraversable(NorthEastCell(l)))
        neighbors.push_back(xyLoc(l.x, l.y - 1));

    // South
    if (IsTraversable(SouthWestCell(l)) || IsTraversable(SouthEastCell(l)))
        neighbors.push_back(xyLoc(l.x, l.y + 1));

    // West
    if (IsTraversable(NorthWestCell(l)) || IsTraversable(SouthWestCell(l)))
        neighbors.push_back(xyLoc(l.x - 1, l.y));

    // East
    if (IsTraversable(NorthEastCell(l)) || IsTraversable(SouthEastCell(l)))
        neighbors.push_back(xyLoc(l.x + 1, l.y));

    return neighbors;
}

void AnyAngleAlgorithm::GetRandomProblem(xyLoc &l1, xyLoc &l2) const
{
    l1 = valid_corner_locations_[rand() % valid_corner_locations_.size()];
    l2 = valid_corner_locations_[rand() % valid_corner_locations_.size()];
}

bool AnyAngleAlgorithm::LineOfSight(xyLoc l1, xyLoc l2)
{
#ifdef ANY_ANGLE_STATISTICS
    num_los_checks_++;
#endif

    // This line of sight check uses only integer values. First it checks whether the movement along the x or the y axis is longer and moves along the longer
    // one cell by cell. dx and dy specify how many cells to move in each direction. Suppose dx is longer and we are moving along the x axis. For each
    // cell we pass in the x direction, we increase variable f by dy, which is initially 0. When f >= dx, we move along the y axis and set f -= dx. This way,
    // after dx movements along the x axis, we also move dy moves along the y axis.

    // x and y values correspond to corners, not cells.
    int x1 = l1.x; // Originate from this cell.
    int y1 = l1.y;

    int x2 = l2.x; // Move to this cell.
    int y2 = l2.y;

    int dy = l2.y - l1.y;
    int dx = l2.x - l1.x;

    int f = 0;
    int sy, sx; // Direction of movement. Value can be either 1 or -1.

    // The x and y locations correspond to corners, not cells. We might need to check different surrounding cells depending on the direction we do the
    // line of sight check. The following values are usedto determine which cell to check to see if it is unblocked.
    int x_offset, y_offset;

    if (dy < 0) {
        dy = -dy;
        sy = -1;
        y_offset = 0; // Cell is to the North
    }
    else {
        sy = 1;
        y_offset = 1; // Cell is to the South
    }

    if (dx < 0) {
        dx = -dx;
        sx = -1;
        x_offset = 0; // Cell is to the West
    }
    else {
        sx = 1;
        x_offset = 1; // Cell is to the East
    }

    if (dx >= dy) { // Move along the x axis and increment/decrement y when f >= dx.
        while (x1 != x2) {
            f = f + dy;
            if (f >= dx) {  // We are changing rows, we might need to check two cells this iteration.
                if (!IsTraversable(xyLoc(x1 + x_offset, y1 + y_offset)))
                    return false;

                y1 = y1 + sy;
                f = f - dx;
            }

            if (f != 0) {   // If f == 0, then we are crossing the row at a corner point and we don't need to check both cells.
                if (!IsTraversable(xyLoc(x1 + x_offset, y1 + y_offset)))
                    return false;
            }

            if (dy == 0) {  // If we are moving along a horizontal line, either the north or the south cell should be unblocked.
                if (!IsTraversable(xyLoc(x1 + x_offset, y1)) && !IsTraversable(xyLoc(x1 + x_offset, y1 + 1)))
                    return false;
            }

            x1 += sx;
        }
    }

    else {  //if (dx < dy). Move along the y axis and increment/decrement x when f >= dy.
        while (y1 != y2) {
            f = f + dx;
            if (f >= dy) {
                if (!IsTraversable(xyLoc(x1 + x_offset, y1 + y_offset)))
                    return false;

                x1 = x1 + sx;
                f = f - dy;
            }

            if (f != 0) {
                if (!IsTraversable(xyLoc(x1 + x_offset, y1 + y_offset)))
                    return false;
            }

            if (dx == 0) {
                if (!IsTraversable(xyLoc(x1, y1 + y_offset)) && !IsTraversable(xyLoc(x1 + 1, y1 + y_offset)))
                    return false;
            }

            y1 += sy;
        }
    }
    return true;
}
bool AnyAngleAlgorithm::LineOfSight(xyLocCont l1, xyLocCont l2)
{
#ifdef ANY_ANGLE_STATISTICS
    num_los_checks_++;
#endif

    // Similar implementation to LineOfSight(xyLoc, xyLoc). But, this time, we might start and end with a non-zero 'f' and 'f' needs to be a float.
    // This function assumes that l1 and l2 do not lie in the interiors of obstacles.

    double dy = l2.y - l1.y;
    double dx = l2.x - l1.x;

    // Decide along which axis we iterate.
    bool iterate_over_x = fabs(dx) >= fabs(dy);

    // Swap l1 and l2 if necessary to make sure that l1.x <= l2.x if we are iterating over the x axis, or l1.y <= l2.y if we are iterating over the y axis.
    if ((iterate_over_x && dx < 0) || (!iterate_over_x && dy < 0)) {
        xyLocCont temp = l2;
        l2 = l1;
        l1 = temp;
        dy = -dy;
        dx = -dx;
    }

    if (iterate_over_x) {
        double df = fabs(dy / dx); // Each time x is increased by 1, f will be adjusted by this amount.

        // We want to iterate over integer values of x. However, x1 and/or x2 might not be integers. Therefore, we try to stretch the endpoints of the line
        // (l1,l2) to make x1 and x2 integers. However, if stretching an endpoint causes the line to intersect with the interior of a blocked cell it does not
        // originally intersect, the line-of-sight check would fail, even if l1 and l2 have line-of-sight. If this is the case, we shorten the line at the
        // problematic endpoint(s).

        // Stretch the beginning of the line.
        int x1 = floor(l1.x + EPSILON); // If l1.x is already an integer, then x1 = l1.x.
        double y1 = GetIntersectingY(l1.x, l1.y, l2.x, l2.y, x1);

        // If stretching the beginning of the line caused the line to intersect with another cell, shorten the beginning of the line instead.
        if (IntervalContainsInteger(y1, l1.y)) {
            x1 = ceil(l1.x);
            y1 = GetIntersectingY(l1.x, l1.y, l2.x, l2.y, x1);
        }

        // Stretch the end of the line.
        int x2 = ceil(l2.x - EPSILON); // If l2.x is already an integer, then x2 = l2.x.
        double y2 = GetIntersectingY(l1.x, l1.y, l2.x, l2.y, x2);

        // If stretching the end of the line caused the line to intersect with another cell, shorten the end of the line instead.
        if (IntervalContainsInteger(y2, l2.y)) {
            x2 = floor(l2.x);
            y2 = GetIntersectingY(l1.x, l1.y, l2.x, l2.y, x2);
        }

        // Now, x1 and x2 are integers. The rest is similar to line-of-sight checks with xyLocs, except for the following:
        // Determine the starting value of f and y, and direction of movement for y
        double f;
        int y;
        int sy; // The direction we update the y coordinates (-1 or +1).
        int y_offset; // Used to determine which cell to check to see if it is unblocked.

        if (dy < 0) {
            sy = -1;
            y_offset = 0; // Cell is to the North
            y = ceil(y1 - EPSILON);
        }
        else {
            sy = 1;
            y_offset = 1; // Cell is to the South
            y = floor(y1 + EPSILON);
        }
        f = fabs(y1 - (double) y);

        while (x1 != x2) {
            f = f + df;
            if (f >= 1) {   // We are changing rows, we might need to check two cells this iteration.
                if (!IsTraversable(xyLoc(x1 + 1, y + y_offset)))
                    return false;

                y = y + sy;
                f = f - 1;
            }

            if (f > EPSILON) {  // If f == 0, then we are crossing the row at a corner point and we don't need to check both cells.
                if (!IsTraversable(xyLoc(x1 + 1, y + y_offset)))
                    return false;
            }

            if (fabs(dy) < EPSILON) {   // If we are moving along a horizontal line, either the north or the south cell should be unblocked.
                if (!IsTraversable(xyLoc(x1 + 1, y)) && !IsTraversable(xyLoc(x1 + 1, y + 1)))
                    return false;
            }

            x1 += 1;
        }
    }

    else {  //if we are iterating over y
        // Same as the case with x. Except we swap x and y.

        double df = fabs(dx / dy); // Each time y is increased by 1, f will be adjusted by this amount

        // Stretch the beginning of the line.
        int y1 = floor(l1.y + EPSILON); // If l1.y is already an integer, then y1 = l1.y.
        double x1 = GetIntersectingX(l1.x, l1.y, l2.x, l2.y, y1);

        // If stretching the beginning of the line caused the line to intersect with another cell, shorten the beginning of the line instead.
        if (IntervalContainsInteger(x1, l1.x)) {
            y1 = ceil(l1.y);
            x1 = GetIntersectingX(l1.x, l1.y, l2.x, l2.y, y1);
        }

        // Stretch the end of the line.
        int y2 = ceil(l2.y - EPSILON); // If l2.y is already an integer, then y2 = l2.y.
        double x2 = GetIntersectingX(l1.x, l1.y, l2.x, l2.y, y2);

        // If stretching the end of the line caused the line to intersect with another cell, shorten the end of the line instead.
        if (IntervalContainsInteger(x2, l2.x)) {
            y2 = floor(l2.y);
            x2 = GetIntersectingX(l1.x, l1.y, l2.x, l2.y, y2);
        }

        // Determine the starting value of f and x, and direction of movement for x.
        double f;
        int x;
        int sx; // The direction we update the x coordinates (-1 or +1).
        int x_offset; // Used to determine which cell to check to see if it is unblocked.

        if (dx < 0) {
            sx = -1;
            x_offset = 0; // Cell is to the North
            x = ceil(x1 - EPSILON);
        }
        else {
            sx = 1;
            x_offset = 1; // Cell is to the South
            x = floor(x1 + EPSILON);
        }
        f = fabs(x1 - (double) x);

        while (y1 != y2) {
            f = f + df;
            if (f >= 1) {   // We are changing columns, we might need to check two cells this iteration.
                if (!IsTraversable(xyLoc(x + x_offset, y1 + 1)))
                    return false;

                x = x + sx;
                f = f - 1;
            }

            if (f > EPSILON) {  // If f == 0, then we are crossing the column at a corner point and we don't need to check both cells.
                if (!IsTraversable(xyLoc(x + x_offset, y1 + 1)))
                    return false;
            }

            if (fabs(dx) < EPSILON) {   // If we are moving along a vertical line, either the north or the south cell should be unblocked.
                if (!IsTraversable(xyLoc(x, y1 + 1)) && !IsTraversable(xyLoc(x + 1, y1 + 1)))
                    return false;
            }

            y1 += 1;
        }
    }
    return true;
}

cost AnyAngleAlgorithm::SmoothPath(const std::vector<xyLoc> &path, std::vector<xyLoc> &smoothed_path)
{
    smoothed_path.clear();

    if (path.empty())
        return INFINITE_COST;

    // Add the start.
    smoothed_path.push_back(path[0]);

    // Go over all locations on the original path, in order.
    for (unsigned int i = 1; i < path.size(); i++) {
        if (!LineOfSight(smoothed_path.back(), path[i])) // If there is no line of sight to the last location that is added to the smoothed path..
        {
            smoothed_path.push_back(path[i - 1]); // ..add the i-1st location on the original path to the smoothed path.
        }
    }

    // Add the goal.
    if (!(smoothed_path.back().x == path.back().x && smoothed_path.back().y == path.back().y))
        smoothed_path.push_back(path.back());

    // Compute path cost.
    cost c = 0;
    for (unsigned int i = 1; i < smoothed_path.size(); i++)
        c += EuclideanDistance(smoothed_path[i-1], smoothed_path[i]);

    return c;
}
cost AnyAngleAlgorithm::SmoothPath(const std::vector<xyLocCont> & path, std::vector<xyLocCont> &smoothed_path)
{
    smoothed_path.clear();

    if (path.empty())
        return INFINITE_COST;

    // Add the start.
    smoothed_path.push_back(path[0]);

    // Go over all locations on the original path, in order.
    for (unsigned int i = 1; i < path.size(); i++) {
        if (!LineOfSight(smoothed_path.back(), path[i])) // If there is no line of sight to the last location that is added to the smoothed path..
        {
            smoothed_path.push_back(path[i - 1]); // ..add the i-1st location on the original path to the smoothed path.
        }
    }

    // Add the goal.
    if (fabs(smoothed_path.back().x - path.back().x) > EPSILON || fabs(smoothed_path.back().y - path.back().y) > EPSILON)
        smoothed_path.push_back(path.back());

    // Compute path cost.
    cost c = 0;
    for (unsigned int i = 1; i < smoothed_path.size(); i++)
        c += EuclideanDistance(smoothed_path[i-1], smoothed_path[i]);

    return c;
}


#ifdef ANY_ANGLE_STATISTICS
void AnyAngleAlgorithm::SetupStatistics()
{
    statistics_ = new AnyAngleStatistics(GetName());

    // Add the suffix '_PS' for reporting statistics with smoothing.
    std::string id_with_smoothing = GetName() + "_PS";
    statistics_with_smoothing_ = new AnyAngleStatistics(id_with_smoothing);
}
void AnyAngleAlgorithm::SetStatisticsFiles(std::string mapname)
{
    if (statistics_ != NULL)
        statistics_->OpenOutputFiles(mapname);
    if (statistics_with_smoothing_ != NULL)
        statistics_with_smoothing_->OpenOutputFiles(mapname);
}
void AnyAngleAlgorithm::PrintStatistics()
{
    if (statistics_ != NULL) {
        PrintAdditionalStatistics(statistics_);
        statistics_ -> ReportAllStatistics();
    }
    if (statistics_with_smoothing_ != NULL) {
        PrintAdditionalStatistics(statistics_with_smoothing_);
        statistics_with_smoothing_ -> ReportAllStatistics();
    }
}

cost AnyAngleAlgorithm::EvaluatePath(const std::vector<xyLoc> & path, const bool validate_path)
{
    int num_los_checks_so_far = num_los_checks_;    // Do not include the line-of-sight check for path validation in statistics.
    bool valid_path = true;

    if (validate_path) {
        if (path.size() == 0) {
            valid_path = false;
        }
        else {
            // First and last locations on the path should be start and goal
            if (path[0].x != from_.x || path[0].y != from_.y)
                valid_path = false;

            if (path.back().x != to_.x || path.back().y != to_.y)
                valid_path = false;
        }
    }

    // Compute the path cost and validate the path.
    cost c = 0;
    for (unsigned int i = 0; i + 1 < path.size(); i++) {
        if (validate_path && !LineOfSight(path[i], path[i + 1]))
            valid_path = false;
        c += EuclideanDistance(path[i], path[i + 1]);
    }

    num_los_checks_ = num_los_checks_so_far; // Do not include the line-of-sight check for path validation in statistics.

    // Count the number of heading changes
    num_heading_changes_ = 0;
    num_freespace_heading_changes_ = 0;
    num_taut_corner_heading_changes_ = 0;
    num_non_taut_corner_heading_changes_ = 0;

    for (unsigned int i = 1; i + 1 < path.size(); i++) {
        if (!CoLinear(path[i - 1].x, path[i - 1].y, path[i].x, path[i].y, path[i + 1].x, path[i + 1].y)) {
            num_heading_changes_++;

            // If the turning point is not at a convex corner of an obstacle
            if (!IsConvexCorner(path[i]))
                num_freespace_heading_changes_++;

            // If the turn at the convex corner of an obstacle produces a taut path
            else if (IsTautCornerTurn(path[i - 1].x, path[i - 1].y, path[i].x, path[i].y, path[i + 1].x, path[i + 1].y))
                num_taut_corner_heading_changes_++;

            // If the turn at the convex corner of an obstacle does not produce a taut path
            else
                num_non_taut_corner_heading_changes_++;
        }
    }

    if (!valid_path)
        return INFINITE_COST;
    else
        return c;
}
cost AnyAngleAlgorithm::EvaluatePath(const std::vector<xyLocCont> &path, const bool validate_path)
{
    int num_los_checks_so_far = num_los_checks_;    // Do not include the line-of-sight check for path validation in statistics.
    bool valid_path = true;

    if (validate_path) {
        if (path.size() == 0) {
            valid_path = false;
        }
        else {
            // First and last locations on the path should be start and goal
            if (!(path[0] == xyLocCont(from_)))
                valid_path = false;

            if (!(path.back() == xyLocCont(to_)))
                valid_path = false;
        }
    }

    // Compute the path cost and validate the path.
    cost c = 0;
    for (unsigned int i = 0; i + 1 < path.size(); i++) {
        if (validate_path && !LineOfSight(path[i], path[i + 1]))
            valid_path = false;
        c += EuclideanDistance(path[i], path[i + 1]);
    }

    num_los_checks_ = num_los_checks_so_far;    // Do not include the line-of-sight check for path validation in statistics.

    // Count the number of direction changes
    num_heading_changes_ = 0;
    num_freespace_heading_changes_ = 0;
    num_taut_corner_heading_changes_ = 0;
    num_non_taut_corner_heading_changes_ = 0;

    for (unsigned int i = 1; i + 1 < path.size(); i++) {
        if (!CoLinear(path[i - 1].x, path[i - 1].y, path[i].x, path[i].y, path[i + 1].x, path[i + 1].y)) {
            num_heading_changes_++;

            // If the turning point is not at a corner
            if (fabs((double) round(path[i].x) - path[i].x) > EPSILON || fabs((double) round(path[i].y) - path[i].y) > EPSILON)
                num_freespace_heading_changes_++;

            // If the turning point is not at a convex corner of an obstacle
            else if (!IsConvexCorner(round(path[i].x), round(path[i].y)))
                num_freespace_heading_changes_++;

            // If the turn at the convex corner of an obstacle produces a taut path
            else if (IsTautCornerTurn(path[i - 1].x, path[i - 1].y, round(path[i].x), round(path[i].y), path[i + 1].x, path[i + 1].y))
                num_taut_corner_heading_changes_++;

            // If the turn at the convex corner of an obstacle does not produce a taut path
            else
                num_non_taut_corner_heading_changes_++;
        }
    }

    if (!valid_path)
        return INFINITE_COST;
    else
        return c;
}

double AnyAngleAlgorithm::NormalizeAngle(const double theta) const
{
    double normalized_theta = fmod(theta, 360);
    if (normalized_theta < 0)
        normalized_theta += 360;

    return normalized_theta;
}
double AnyAngleAlgorithm::GetAngle(const double x1, const double y1, const double x2, const double y2) const
{
    // Normalize so that the vector is from (0,0) to (x,y)
    double x = x2 - x1;
    double y = y2 - y1;

    // -y because the y coordinate for the map grows as we move South
    double theta = atan2(-y,x) * 180 / 3.14159265358979323846;

    return NormalizeAngle(theta);
}
bool AnyAngleAlgorithm::IsTautCornerTurn(const double x1, const double y1, const int x2, const int y2, const double x3, const double y3) const
{
//   std::cout<<x2<<"\t"<<y2<<std::endl;

    // Compute the angle of the vector (x2, y2) -> (x1, y1)
    double theta_1 = GetAngle(x2, y2, x1, y1);
//    std::cout<<"Theta1 = "<<theta_1<<std::endl;

    // Compute the angle of the vector (x2, y2) -> (x3, y3)
    double theta_2 = GetAngle(x2, y2, x3, y3);
//    std::cout<<"Theta2 = "<<theta_2<<std::endl;

    double theta_diff = NormalizeAngle(theta_2 - theta_1);
//    std::cout<<"Theta diff = "<<theta_diff<<std::endl;

    double theta_bisector;

    if (theta_diff < 180)
        theta_bisector = NormalizeAngle(theta_1 + theta_diff/2);
    else
        theta_bisector = NormalizeAngle(theta_1 + theta_diff/2 + 180);

//    std::cout<<"Bisector = "<<theta_bisector<<std::endl;

    switch ((int) theta_bisector/90) {

    case 0:
        return !IsTraversable(NorthEastCell(x2,y2));

    case 1:
        return !IsTraversable(NorthWestCell(x2,y2));

    case 2:
        return !IsTraversable(SouthWestCell(x2,y2));

    case 3:
        return !IsTraversable(SouthEastCell(x2,y2));

    default:
        return false;
    }

    return false;
}

void AnyAngleAlgorithm::StartStatistics(const xyLoc from, const xyLoc to)
{
    // Store the start and goal corners for path validation.
    from_ = from;
    to_ = to;

    // Reset all statistics
    timer_.StartTimer();
    num_expansions_ = 0;
    num_generated_ = 0;
    num_percolations_ = 0;
    num_los_checks_ = 0;
    num_heading_changes_ = 0;
    num_freespace_heading_changes_ = 0;
    num_taut_corner_heading_changes_ = 0;
    num_non_taut_corner_heading_changes_ = 0;
    elapsed_time_ = 0;
}
void AnyAngleAlgorithm::ReportStatistics(const std::vector<xyLoc> &path, AnyAngleStatistics* stats, const bool validate_path)
{
    elapsed_time_ += timer_.EndTimer(); // Time passed since StartStatistics is called

    // Validate the path; does not increment line-of-sight checks.
    cost c = EvaluatePath(path, validate_path);

    // Report a failed search.
    if (!(c < INFINITE_COST)) {
        printf("Failed search from (%d, %d) to (%d, %d)\n", from_.x, from_.y, to_.x, to_.y);
    }

    // Report search statistics.
    stats->AddSearchData(elapsed_time_, c, num_expansions_, num_generated_, num_percolations_, num_los_checks_,
            num_heading_changes_, num_freespace_heading_changes_, num_taut_corner_heading_changes_, num_non_taut_corner_heading_changes_);

    timer_.StartTimer();    // Start the timer again in case we also want to report smoothing statistics.
}
void AnyAngleAlgorithm::ReportStatistics(const std::vector<xyLocCont> &path, AnyAngleStatistics* stats, const bool validate_path)
{
    elapsed_time_ += timer_.EndTimer(); // Time passed since StartStatistics is called

    // Validate the path; does not increment line-of-sight checks.
    cost c = EvaluatePath(path, validate_path);

    // Report a failed search.
    if (!(c < INFINITE_COST)) {
        printf("Failed search from (%d, %d) to (%d, %d)\n", from_.x, from_.y, to_.x, to_.y);
    }

    // Report search statistics.
    stats->AddSearchData(elapsed_time_, c, num_expansions_, num_generated_, num_percolations_, num_los_checks_,
            num_heading_changes_, num_freespace_heading_changes_, num_taut_corner_heading_changes_, num_non_taut_corner_heading_changes_);

    timer_.StartTimer();    // Start the timer again in case we also want to report smoothing statistics.
}
#endif

#ifdef ANY_ANGLE_RUNNING_IN_HOG
void AnyAngleAlgorithm::DrawLine(const MapEnvironment *env, const float x1, const float y1, const float x2, const float y2) const
{
    env->GLDrawColoredLine(x1 - 0.5, y1 - 0.5, x2 - 0.5, y2 - 0.5);
}
void AnyAngleAlgorithm::DrawPoint(const MapEnvironment *env, const float x, const float y, const int priority) const
{
    env->OpenGLPriorityDraw(x - 0.5, y - 0.5, priority);
}
void AnyAngleAlgorithm::DrawPath(const MapEnvironment *env, const std::vector<xyLoc> &path, const bool draw_mid_points, const int priority) const
{
    if (draw_mid_points)
        for (unsigned int i = 0; i < path.size(); i++)
            DrawPoint(env, path[i], priority);
    for (unsigned int i = 1; i < path.size(); i++)
        DrawLine(env, path[i - 1], path[i]);
}
void AnyAngleAlgorithm::DrawPath(const MapEnvironment *env, const std::vector<xyLocCont> &path, const bool draw_mid_points, const int priority) const
{
    if (draw_mid_points)
        for (unsigned int i = 0; i < path.size(); i++)
            DrawPoint(env, path[i], priority);
    for (unsigned int i = 1; i < path.size(); i++)
        DrawLine(env, path[i - 1], path[i]);
}

void AnyAngleAlgorithm::ShowPath(const MapEnvironment *env, float r, float g, float b)
{
    env->SetColor(r, g, b);

    if (!UsingXYLocCont()) {
        DrawPath(env, this->xyloc_path_, true);
    }

    else {
        DrawPath(env, this->xyloc_cont_path_, true);
    }
}
void AnyAngleAlgorithm::ShowSmoothedPath(const MapEnvironment *env, float r, float g, float b)
{
    if (!ShouldSmoothPaths())
        return;

    env->SetColor(r, g, b);

    if (!UsingXYLocCont()) {
        DrawPath(env, this->smoothed_xyloc_path_, false);
    }

    else {
        DrawPath(env, this->smoothed_xyloc_cont_path_, false);
    }
}
#endif
