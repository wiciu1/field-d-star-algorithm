import tkinter as tk
import heapq
import math
from dataclasses import dataclass

# Field D* Algorithm Implementation
# Based on D* Lite proposed by Sebastian (https://github.com/SebastianPPP/d_star_algorithm)
# Also based on documentation: https://www.ri.cmu.edu/pub_files/pub4/ferguson_david_2005_3/ferguson_david_2005_3.pdf

# constants
TOLERANCE = 0.01
EPS = 1e-9


@dataclass(frozen=True)
class InterpChoice:
    """Result of ComputeCost for one consecutive neighbor pair."""
    cost: float
    edge_param: float
    mode: str


@dataclass(frozen=True)
class SuccessorChoice:
    """Best successor geometry for a corner node."""
    cost: float
    pair_index: int
    edge_param: float
    mode: str
    s1: tuple[int, int]
    s2: tuple[int, int]
    c_cell: tuple[int, int]
    b_cell: tuple[int, int]


@dataclass(frozen=True)
class BoundaryCandidate:
    """Candidate waypoint on a cell boundary used during path extraction."""
    point: tuple[float, float]
    point_cost: float
    total_cost: float
    cell: tuple[int, int]
    edge: tuple[tuple[int, int], tuple[int, int]]

class Cell:
    """Grid cell with a uniform traversal cost or an obstacle flag."""
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.is_obstacle = False
        self.rect_id = None


class Node:
    """Corner node used by the Field D* planner."""
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.g = math.inf
        self.rhs = math.inf
        self.key = (math.inf, math.inf)
        self.index = (x, y)
        self.open_version = 0

    # Overloading
    def __lt__(self, other):
        return self.key < other.key

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        return f"Node({self.x},{self.y})"


# Field D* Algorithm class
class FieldDStar:
    """Field D* planner on a uniform 2D grid with corner nodes."""
    def __init__(self, cells, nodes, start_pos, goal_node):
        self.cells = cells
        self.nodes = nodes
        self.cols = len(cells)
        self.rows = len(cells[0])

        self.start_pos = start_pos
        self.goal_node = goal_node

        self.U = []         # PQ (min-heap)
        self.open_entries = {}
        self.path = []      # Path (Float coordinates)
        self.is_optimal = False
        self.k_m = 0.0

        # format: ((dx1, dy1), (dx2, dy2), (dx3, dy3), (dx4, dy4))
        # 1st (dx1, dy1): perpendicular neighbor node (s1)
        # 2nd (dx2, dy2): diagonal neighbor node (s2)
        # 3rd (dx3, dy3): primary cell (c) that the path cuts through (bounded by nodes s, s1, s2)
        # 4th (dx4, dy4): boundary cell (b) that shares the perpendicular edge (s-s1) with cell c
        self.pairs = [
            ((1, 0), (1, 1), (0, 0), (0, -1)),
            ((0, 1), (1, 1), (0, 0), (-1, 0)),
            ((-1, 0), (-1, 1), (-1, 0), (-1, -1)),
            ((0, 1), (-1, 1), (-1, 0), (0, 0)),
            ((-1, 0), (-1, -1), (-1, -1), (-1, 0)),
            ((0, -1), (-1, -1), (-1, -1), (0, -1)),
            ((1, 0), (1, -1), (0, -1), (0, 0)),
            ((0, -1), (1, -1), (0, -1), (-1, -1))
        ]

    def get_cell_cost(self, x, y):
        """Returns the traversal cost per unit length for a cell."""
        if 0 <= x < self.cols and 0 <= y < self.rows:
            return math.inf if self.cells[x][y].is_obstacle else 1.0
        return math.inf

    def is_finite(self, value):
        """Treats infinities and NaNs as invalid costs."""
        return math.isfinite(value) and not math.isnan(value)

    # Implementation of ComputeCost pseudocode from document
    def compute_cost_choice(self, c, b, g1, g2):
        """
        Calculates the traversal cost from node s through the edge formed by nodes s1 and s2.

        This mirrors the paper's ComputeCost routine and returns both the scalar cost
        and the edge parameter of the optimal interpolation point.

        :param c: Cost of the primary cell that the path cuts through
        :param b: Cost of the boundary cell
        :param g1: Computed cost to reach the perpendicular node (s1)
        :param g2: Computed cost to reach the diagonal node (s2)
        """

        if any(math.isnan(v) for v in (c, b, g1, g2)):
            return InterpChoice(math.inf, 0.0, "nan")
        if min(c, b) == math.inf:
            return InterpChoice(math.inf, 0.0, "blocked")
        if g1 <= g2 + EPS:
            return InterpChoice(min(c, b) + g1, 0.0, "straight")

        f = g1 - g2
        if f <= b + EPS:
            # Going diagonally
            if c <= f + EPS:
                return InterpChoice(c * math.sqrt(2) + g2, 1.0, "diag")
            else:
                # Calculate absolute minimum of function (best case)
                val = c ** 2 - f ** 2
                y = f / math.sqrt(val) if val > EPS else 1.0
                y = min(y, 1.0)
                return InterpChoice(c * math.sqrt(1 + y ** 2) + f * (1 - y) + g2, y, "cut")
        else:
            if c <= b + EPS:
                return InterpChoice(c * math.sqrt(2) + g2, 1.0, "diag")
            else:
                val = c ** 2 - b ** 2
                x_sub = b / math.sqrt(val) if val > EPS else 1.0
                x = 1.0 - min(x_sub, 1.0)
                return InterpChoice(c * math.sqrt(1 + (1 - x) ** 2) + b * x + g2, x, "boundary")

    def calculate_key(self, s):
        """Calculates the D* Lite priority key for a node."""
        h = math.hypot(s.x - self.start_pos[0], s.y - self.start_pos[1])
        return (min(s.g, s.rhs) + h + self.k_m, min(s.g, s.rhs))

    def best_successor_choice(self, s):
        """Returns the best interpolated successor configuration for node s."""
        if s == self.goal_node:
            return SuccessorChoice(
                cost=0.0,
                pair_index=-1,
                edge_param=0.0,
                mode="goal",
                s1=(s.x, s.y),
                s2=(s.x, s.y),
                c_cell=(-1, -1),
                b_cell=(-1, -1),
            )
        best = None

        for pair_index, ((dx1, dy1), (dx2, dy2), (dcx, dcy), (dbx, dby)) in enumerate(self.pairs):
            s1x, s1y = s.x + dx1, s.y + dy1
            s2x, s2y = s.x + dx2, s.y + dy2

            if (0 <= s1x <= self.cols and
                0 <= s1y <= self.rows and
                0 <= s2x <= self.cols and
                0 <= s2y <= self.rows):
                g1 = self.nodes[s1x][s1y].g
                g2 = self.nodes[s2x][s2y].g
                c_cost = self.get_cell_cost(s.x + dcx, s.y + dcy)
                b_cost = self.get_cell_cost(s.x + dbx, s.y + dby)

                choice = self.compute_cost_choice(c_cost, b_cost, g1, g2)
                if best is None or choice.cost + EPS < best.cost:
                    best = SuccessorChoice(
                        cost=choice.cost,
                        pair_index=pair_index,
                        edge_param=choice.edge_param,
                        mode=choice.mode,
                        s1=(s1x, s1y),
                        s2=(s2x, s2y),
                        c_cell=(s.x + dcx, s.y + dcy),
                        b_cell=(s.x + dbx, s.y + dby),
                    )
        if best is None:
            return SuccessorChoice(
                cost=math.inf,
                pair_index=-1,
                edge_param=0.0,
                mode="unreachable",
                s1=(s.x, s.y),
                s2=(s.x, s.y),
                c_cell=(-1, -1),
                b_cell=(-1, -1),
            )
        return best

    def compute_min_rhs(self, s):
        """Computes the minimum rhs value for a given node."""
        return self.best_successor_choice(s).cost

    def get_neighbors(self, s):
        """Returns the 8-neighborhood of a corner node."""
        neighbbrs = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0: continue
                nx, ny = s.x + dx, s.y + dy
                if 0 <= nx <= self.cols and 0 <= ny <= self.rows:
                    neighbbrs.append(self.nodes[nx][ny])
        return neighbbrs

    def predecessors(self, s):
        """On this grid the local predecessor and successor sets are identical."""
        return self.get_neighbors(s)

    def push_open(self, s):
        """Inserts or refreshes a node in OPEN using logical deletion."""
        s.key = self.calculate_key(s)
        s.open_version += 1
        self.open_entries[s.index] = (s.key, s.open_version)
        heapq.heappush(self.U, (s.key, s.open_version, s))

    def remove_from_open(self, s):
        """Marks a node as absent from OPEN."""
        self.open_entries.pop(s.index, None)

    def top_key(self):
        """Returns the current minimum valid key in OPEN."""
        while self.U:
            key, version, node = self.U[0]
            active = self.open_entries.get(node.index)
            if active == (key, version):
                return key
            heapq.heappop(self.U)
        return (math.inf, math.inf)

    def pop_valid(self):
        """Pops the minimum valid OPEN entry and discards stale heap entries."""
        while self.U:
            key, version, node = heapq.heappop(self.U)
            active = self.open_entries.get(node.index)
            if active == (key, version):
                self.open_entries.pop(node.index, None)
                return key, node
        return (math.inf, math.inf), None

    def update_node(self, s):
        """Recomputes rhs and synchronizes the node with OPEN."""
        if s != self.goal_node:
            s.rhs = self.compute_min_rhs(s)
        if abs(s.g - s.rhs) <= EPS or (math.isinf(s.g) and math.isinf(s.rhs)):
            self.remove_from_open(s)
        else:
            self.push_open(s)

    def Initialize(self):
        """Initializes the planner state for a fresh search."""
        self.U = []
        self.open_entries = {}
        self.is_optimal = False
        self.k_m = 0.0

        for x in range(self.cols + 1):
            for y in range(self.rows + 1):
                self.nodes[x][y].g = math.inf
                self.nodes[x][y].rhs = math.inf
                self.nodes[x][y].key = (math.inf, math.inf)
                self.nodes[x][y].open_version = 0

        self.goal_node.rhs = 0.0
        self.push_open(self.goal_node)

    def compute_shortest_path(self):
        """
        Repairs inconsistent nodes until the start reference is locally consistent.

        The loop follows the D* Lite stopping condition instead of running until OPEN
        becomes empty.
        """
        start_node = self.closest_node_to_start()
        while self.top_key() < self.calculate_key(start_node) or abs(start_node.g - start_node.rhs) > EPS:
            old_key, u = self.pop_valid()
            if u is None:
                break
            curr_key = self.calculate_key(u)

            if old_key < curr_key:
                self.push_open(u)
            elif u.g > u.rhs:
                u.g = u.rhs
                for n in self.predecessors(u):
                    self.update_node(n)
            else:
                u.g = math.inf
                self.update_node(u)
                for n in self.predecessors(u):
                    self.update_node(n)

        self.is_optimal = True
        self.extract_path()

    def closest_node_to_start(self):
        """Returns the corner node currently used as the D* Lite start reference."""
        sx, sy = self.start_pos
        nx = min(max(int(round(sx)), 0), self.cols)
        ny = min(max(int(round(sy)), 0), self.rows)
        return self.nodes[nx][ny]

    def point_on_edge(self, edge, t):
        """Interpolates a point on an axis-aligned edge."""
        (x1, y1), (x2, y2) = edge
        return (x1 + (x2 - x1) * t, y1 + (y2 - y1) * t)

    def interpolate_edge_cost(self, edge, t):
        """Interpolates cost-to-go on an edge from its endpoint g-values."""
        (x1, y1), (x2, y2) = edge
        g1 = self.nodes[x1][y1].g
        g2 = self.nodes[x2][y2].g
        if not self.is_finite(g1) and not self.is_finite(g2):
            return math.inf
        if not self.is_finite(g1):
            return g2
        if not self.is_finite(g2):
            return g1
        return g1 + (g2 - g1) * t

    def edge_minimizer(self, point, cell_cost, edge):
        """
        Analytically minimizes travel-to-edge plus linear edge interpolation.
        Kept as a local utility for point evaluation experiments and debugging.
        """
        if not self.is_finite(cell_cost):
            return None

        (x1, y1), (x2, y2) = edge
        g1 = self.nodes[x1][y1].g
        g2 = self.nodes[x2][y2].g
        if not self.is_finite(g1) and not self.is_finite(g2):
            return None

        px, py = point
        if x1 == x2:
            base_coord = y1
            target_coord = py
            orth_dist = abs(px - x1)
        else:
            base_coord = x1
            target_coord = px
            orth_dist = abs(py - y1)

        slope = (g2 - g1) if self.is_finite(g1) and self.is_finite(g2) else 0.0
        anchor = target_coord - base_coord
        candidates = [0.0, 1.0]

        if abs(slope) < cell_cost - EPS and orth_dist > EPS:
            denom = cell_cost ** 2 - slope ** 2
            if denom > EPS:
                delta = abs(slope) * orth_dist / math.sqrt(denom)
                coord = anchor - math.copysign(delta, slope)
                candidates.append(coord)
        elif orth_dist <= EPS:
            if slope > 0:
                candidates.append(0.0)
            elif slope < 0:
                candidates.append(1.0)
            else:
                candidates.append(min(max(anchor, 0.0), 1.0))

        best = None
        for raw_t in candidates:
            t = min(max(raw_t, 0.0), 1.0)
            edge_point = self.point_on_edge(edge, t)
            edge_cost = self.interpolate_edge_cost(edge, t)
            if not self.is_finite(edge_cost):
                continue
            total = cell_cost * math.hypot(edge_point[0] - px, edge_point[1] - py) + edge_cost
            if best is None or total + EPS < best[0]:
                best = (total, t, edge_point, edge_cost)
        return best

    def cells_containing_point(self, point):
        """Returns all cells incident to a fractional point."""
        px, py = point
        cells = []
        min_cx = int(max(0, math.floor(px - TOLERANCE)))
        max_cx = int(min(self.cols - 1, math.floor(px + TOLERANCE)))
        min_cy = int(max(0, math.floor(py - TOLERANCE)))
        max_cy = int(min(self.rows - 1, math.floor(py + TOLERANCE)))
        for cx in range(min_cx, max_cx + 1):
            for cy in range(min_cy, max_cy + 1):
                cells.append((cx, cy))
        return cells

    def cell_edges(self, cx, cy):
        """Returns the four boundary edges of a cell."""
        return [
            ((cx, cy), (cx + 1, cy)),
            ((cx + 1, cy), (cx + 1, cy + 1)),
            ((cx, cy + 1), (cx + 1, cy + 1)),
            ((cx, cy), (cx, cy + 1)),
        ]

    def cell_corners(self, cx, cy):
        """Returns the four corner nodes of a cell."""
        return [
            self.nodes[cx][cy],
            self.nodes[cx + 1][cy],
            self.nodes[cx + 1][cy + 1],
            self.nodes[cx][cy + 1],
        ]

    def point_in_cell(self, point, cell):
        """Checks whether a point lies in the closed cell footprint."""
        px, py = point
        cx, cy = cell
        return (cx - TOLERANCE <= px <= cx + 1 + TOLERANCE and
                cy - TOLERANCE <= py <= cy + 1 + TOLERANCE)

    def choice_edge(self, choice):
        """Returns the edge associated with a successor choice."""
        return (choice.s1, choice.s2)

    def choice_point(self, choice):
        """Converts a successor choice into a concrete boundary waypoint."""
        return self.point_on_edge(self.choice_edge(choice), choice.edge_param)

    def choice_point_cost(self, choice):
        """Returns the interpolated cost-to-go of the waypoint selected by a choice."""
        return self.interpolate_edge_cost(self.choice_edge(choice), choice.edge_param)

    def direct_goal_cost(self, point):
        """Returns the direct in-cell travel cost to the goal when visible in one cell."""
        for cell in self.cells_containing_point(point):
            if self.point_in_cell((self.goal_node.x, self.goal_node.y), cell):
                cell_cost = self.get_cell_cost(*cell)
                if self.is_finite(cell_cost):
                    return cell_cost * math.hypot(point[0] - self.goal_node.x, point[1] - self.goal_node.y)
        return math.inf

    def local_successor_candidates(self, point):
        """Builds boundary waypoint candidates from corner-node successor geometry.

        Instead of sampling edges, this reuses the interpolation decisions already
        computed for corner nodes by the planner.
        """
        candidates = []
        seen = set()

        for cell in self.cells_containing_point(point):
            cell_cost = self.get_cell_cost(*cell)
            if not self.is_finite(cell_cost):
                continue

            for corner in self.cell_corners(*cell):
                choice = self.best_successor_choice(corner)
                if not self.is_finite(choice.cost) or choice.pair_index < 0:
                    continue

                choice_point = self.choice_point(choice)
                if not self.point_in_cell(choice_point, cell):
                    continue

                if choice.c_cell != cell and choice.b_cell != cell:
                    continue

                # Travel from the current fractional point to the waypoint uses the
                # traversal cost of the cell the segment lies in.
                segment_cost = cell_cost * math.hypot(choice_point[0] - point[0], choice_point[1] - point[1])
                point_cost = self.choice_point_cost(choice)
                if not self.is_finite(point_cost):
                    continue

                key = (round(choice_point[0], 6), round(choice_point[1], 6), cell)
                if key in seen:
                    continue
                seen.add(key)

                candidates.append(
                    BoundaryCandidate(
                        point=choice_point,
                        point_cost=point_cost,
                        total_cost=segment_cost + point_cost,
                        cell=cell,
                        edge=self.choice_edge(choice),
                    )
                )

        goal_cost = self.direct_goal_cost(point)
        if self.is_finite(goal_cost):
            candidates.append(
                BoundaryCandidate(
                    point=(self.goal_node.x, self.goal_node.y),
                    point_cost=0.0,
                    total_cost=goal_cost,
                    cell=(-1, -1),
                    edge=((self.goal_node.x, self.goal_node.y), (self.goal_node.x, self.goal_node.y)),
                )
            )

        return candidates

    def point_cost(self, point, forbidden_edges=None):
        """Returns a one-step lookahead cost-to-goal estimate for a fractional point."""
        forbidden_edges = forbidden_edges or set()
        best = math.inf
        for candidate in self.local_successor_candidates(point):
            normalized_edge = tuple(sorted(candidate.edge))
            if normalized_edge in forbidden_edges:
                continue
            best = min(best, candidate.total_cost)
        return best

    def next_boundary_candidates(self, point):
        """Returns waypoint candidates scored with an extra one-step lookahead."""
        candidates = []
        for candidate in self.local_successor_candidates(point):
            if math.hypot(candidate.point[0] - point[0], candidate.point[1] - point[1]) <= TOLERANCE:
                continue
            forbidden = {tuple(sorted(candidate.edge))}
            lookahead = self.point_cost(candidate.point, forbidden_edges=forbidden)
            total_cost = candidate.total_cost if not self.is_finite(lookahead) else (
                math.hypot(candidate.point[0] - point[0], candidate.point[1] - point[1]) *
                self.get_cell_cost(*candidate.cell) + lookahead
            )
            candidates.append(
                BoundaryCandidate(
                    point=candidate.point,
                    point_cost=candidate.point_cost,
                    total_cost=total_cost,
                    cell=candidate.cell,
                    edge=candidate.edge,
                )
            )
        return candidates

    def extract_path(self):
        """Extracts a path by iteratively selecting boundary waypoints.

        The extractor works on fractional points on cell boundaries rather than
        snapping the path to corner nodes.
        """
        self.path = [(self.start_pos[0], self.start_pos[1])]
        curr_x, curr_y = self.start_pos
        visited = set()

        for _ in range(300):
            if math.hypot(curr_x - self.goal_node.x, curr_y - self.goal_node.y) < TOLERANCE:
                break

            candidates = self.next_boundary_candidates((curr_x, curr_y))
            if not candidates:
                break
            best = min(candidates, key=lambda candidate: (candidate.total_cost, candidate.point_cost))
            best_p = best.point

            rounded = (round(best_p[0], 2), round(best_p[1], 2))
            if rounded in visited:
                break
            visited.add(rounded)

            curr_x, curr_y = best_p
            self.path.append((curr_x, curr_y))

        if math.hypot(curr_x - self.goal_node.x, curr_y - self.goal_node.y) < TOLERANCE:
            if math.hypot(self.path[-1][0] - self.goal_node.x, self.path[-1][1] - self.goal_node.y) > TOLERANCE:
                self.path.append((self.goal_node.x, self.goal_node.y))
        elif self.is_finite(self.goal_node.g):
            goal_distance = math.hypot(curr_x - self.goal_node.x, curr_y - self.goal_node.y)
            if goal_distance <= math.sqrt(2) + TOLERANCE:
                self.path.append((self.goal_node.x, self.goal_node.y))
        return self.path


# GUI
class FieldDStarApp:
    def __init__(self, master, cols=20, rows=15, cell_size=50):
        self.master = master
        master.title("Field D* Simulation")

        self.COLS = cols
        self.ROWS = rows
        self.CELL_SIZE = cell_size

        self.canvas = tk.Canvas(master, width=self.COLS * self.CELL_SIZE, height=self.ROWS * self.CELL_SIZE, bg='white')
        self.canvas.pack()

        # Cells initialization
        self.cells = [[Cell(x, y) for y in range(self.ROWS)] for x in range(self.COLS)]
        # Nodes initialization
        self.nodes = [[Node(x, y) for y in range(self.ROWS + 1)] for x in range(self.COLS + 1)]

        self.start_pos = (1.5, self.ROWS - 1.5)
        self.goal_node = self.nodes[self.COLS - 2][2]


        self.dstar = FieldDStar(self.cells, self.nodes, self.start_pos, self.goal_node)
        self.path_line_id = None
        self.agent_marker = None

        self.draw_grid()
        self.dstar.Initialize()
        self.update_visualization()

        master.bind('<space>', lambda event: self.handle_space_press())
        self.canvas.bind('<Button-1>', self.toggle_obstacle)
        self.canvas.bind('<Button-3>', self.toggle_obstacle)

        tk.Label(master, text="Space: Plan / Agent's move | Click: Add / Remove obstacle").pack()

    def draw_grid(self):
        for x in range(self.COLS):
            for y in range(self.ROWS):
                c = self.cells[x][y]
                x1, y1 = x * self.CELL_SIZE, y * self.CELL_SIZE
                x2, y2 = x1 + self.CELL_SIZE, y1 + self.CELL_SIZE
                c.rect_id = self.canvas.create_rectangle(x1, y1, x2, y2, outline='#ddd', fill='white')

        # Draw Goal
        gx, gy = self.goal_node.x * self.CELL_SIZE, self.goal_node.y * self.CELL_SIZE
        self.canvas.create_oval(gx - 8, gy - 8, gx + 8, gy + 8, fill='purple')
        self.canvas.create_text(gx, gy - 15, text="GOAL", fill='purple', font=("Arial", 10, "bold"))

    def toggle_obstacle(self, event):
        cx, cy = event.x // self.CELL_SIZE, event.y // self.CELL_SIZE
        if 0 <= cx < self.COLS and 0 <= cy < self.ROWS:
            self.cells[cx][cy].is_obstacle = not self.cells[cx][cy].is_obstacle

            # Prevention of putting obstacle in the agent cell
            ax, ay = self.dstar.start_pos
            if cx <= ax <= cx + 1 and cy <= ay <= cy + 1 and self.cells[cx][cy].is_obstacle:
                self.cells[cx][cy].is_obstacle = False
                return

            # If cost of cell has changed, replan its 4 corners
            for nx in (cx, cx + 1):
                for ny in (cy, cy + 1):
                    self.dstar.update_node(self.nodes[nx][ny])

            self.dstar.is_optimal = False
            self.update_visualization()

    def move_agent(self):
        if math.hypot(self.dstar.start_pos[0] - self.dstar.goal_node.x,
                      self.dstar.start_pos[1] - self.dstar.goal_node.y) < 0.1:
            print("Goal achieved")
            return False

        if not hasattr(self.dstar, 'path') or len(self.dstar.path) < 2:
            print("Can't find the math to the goal")
            return False

        # Agent's vector move -> fractional point
        next_pos = self.dstar.path[1]

        # Update of traversed distance and start point in the core algorithm
        dist = math.hypot(next_pos[0] - self.dstar.start_pos[0], next_pos[1] - self.dstar.start_pos[1])
        self.dstar.k_m += dist
        self.dstar.start_pos = next_pos

        print(f"Agents move to fractional vector: ({next_pos[0]:.2f}, {next_pos[1]:.2f})")

        # Generating new path
        self.dstar.extract_path()
        return True

    def handle_space_press(self):
        if not self.dstar.is_optimal:
            self.dstar.compute_shortest_path()
        elif self.dstar.is_optimal:
            self.move_agent()

        self.update_visualization()

    def update_visualization(self):
        for x in range(self.COLS):
            for y in range(self.ROWS):
                c = self.cells[x][y]
                if c.is_obstacle:
                    self.canvas.itemconfig(c.rect_id, fill='#654321')  # Obstacle
                else:
                    self.canvas.itemconfig(c.rect_id, fill='white')

        if self.path_line_id:
            self.canvas.delete(self.path_line_id)

        if hasattr(self.dstar, 'path') and len(self.dstar.path) > 1:
            coords = []
            for px, py in self.dstar.path:
                coords.extend([px * self.CELL_SIZE, py * self.CELL_SIZE])
            self.path_line_id = self.canvas.create_line(*coords, fill='blue', width=4, smooth=False)

        # Agent
        if self.agent_marker:
            self.canvas.delete(self.agent_marker)

        sx, sy = self.dstar.start_pos[0] * self.CELL_SIZE, self.dstar.start_pos[1] * self.CELL_SIZE
        self.agent_marker = self.canvas.create_oval(sx - 6, sy - 6, sx + 6, sy + 6, fill='green')


if __name__ == "__main__":
    ROOT = tk.Tk()
    APP = FieldDStarApp(ROOT, cols=20, rows=15, cell_size=50)
    ROOT.mainloop()
