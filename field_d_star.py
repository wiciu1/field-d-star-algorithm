import tkinter as tk
import heapq
import math

# Field D* Algorithm Implementation
# Based on D* Lite proposed by Sebastian (https://github.com/SebastianPPP/d_star_algorithm)
# Also based on documentation: https://www.ri.cmu.edu/pub_files/pub4/ferguson_david_2005_3/ferguson_david_2005_3.pdf

# constants
TOLERANCE = 0.01

class Cell:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.is_obstacle = False
        self.rect_id = None


class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.g = math.inf
        self.rhs = math.inf
        self.key = (math.inf, math.inf)

    # Overloading
    def __lt__(self, other):
        return self.key < other.key

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        return f"Node({self.x},{self.y})"


# Field D* Algorithm class
class FieldDStar:
    def __init__(self, cells, nodes, start_pos, goal_node):
        self.cells = cells
        self.nodes = nodes
        self.cols = len(cells)
        self.rows = len(cells[0])

        self.start_pos = start_pos
        self.goal_node = goal_node

        self.U = []         # PQ (min-heap)
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
        """Returns the traversal cost of the cell at the specified coordinates"""
        if 0 <= x < self.cols and 0 <= y < self.rows:
            return math.inf if self.cells[x][y].is_obstacle else 1.0
        return math.inf

    # Implementation of ComputeCost pseudocode from document
    def compute_cost(self, c, b, g1, g2):
        """
        Calculates the traversal cost from node s through the edge formed by nodes s1 and s2
        It uses linear interpolation and the Pythagorean theorem to optimize the
        crossing angle through cell c (or along the boundary shared with cell b)

        :param c: Cost of the primary cell that the path cuts through
        :param b: Cost of the boundary cell
        :param g1: Computed cost to reach the perpendicular node (s1)
        :param g2: Computed cost to reach the diagonal node (s2)
        """

        # If obstacle -> can't move
        if min(c, b) == math.inf: return math.inf
        # If going straight is >= going diagonal you shouldn't cut the path
        if g1 <= g2: return min(c, b) + g1

        f = g1 - g2
        if f <= b:
            # Going diagonally
            if c <= f:
                return c * math.sqrt(2) + g2
            else:
                # Calculate absolute minimum of function (best case)
                val = c ** 2 - f ** 2
                y = f / math.sqrt(val) if val > 0.00001 else 1.0 # preventing ZeroDivisionError
                y = min(y, 1.0)
                return c * math.sqrt(1 + y ** 2) + f * (1 - y) + g2
        else:
            if c <= b:
                return c * math.sqrt(2) + g2
            else:
                val = c ** 2 - b ** 2
                x_sub = b / math.sqrt(val) if val > 0.00001 else 1.0
                x = 1.0 - min(x_sub, 1.0)
                return c * math.sqrt(1 + (1 - x) ** 2) + b * x + g2

    def calculate_key(self, s):
        """
        Calculates the PQ key for a given node
        Returns a key (k1, k2)
        """
        h = math.hypot(s.x - self.start_pos[0], s.y - self.start_pos[1])
        return (min(s.g, s.rhs) + h + self.k_m, min(s.g, s.rhs))

    def compute_min_rhs(self, s):
        """ Computes the minimum rhs value for a given node """
        if s == self.goal_node: return 0.0
        min_rhs = math.inf

        for (dx1, dy1), (dx2, dy2), (dcx, dcy), (dbx, dby) in self.pairs:
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

                cost = self.compute_cost(c_cost, b_cost, g1, g2)
                if cost < min_rhs: min_rhs = cost
        return min_rhs

    def get_neighbors(self, s):
        """ Returns all neighbors of given node s """
        neighbbrs = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0: continue
                nx, ny = s.x + dx, s.y + dy
                if 0 <= nx <= self.cols and 0 <= ny <= self.rows:
                    neighbbrs.append(self.nodes[nx][ny])
        return neighbbrs

    def update_node(self, s):
        """ Updates the rhs value of a node and manages its status in the priority queue """
        s.rhs = self.compute_min_rhs(s)
        if s.g != s.rhs:
            s.key = self.calculate_key(s)
            heapq.heappush(self.U, (s.key, id(s), s))

    def Initialize(self):
        """ Starting initialization """
        self.U = []
        self.is_optimal = False
        self.k_m = 0.0

        for x in range(self.cols + 1):
            for y in range(self.rows + 1):
                self.nodes[x][y].g = math.inf
                self.nodes[x][y].rhs = math.inf
                self.nodes[x][y].key = (math.inf, math.inf)

        self.goal_node.rhs = 0.0
        self.goal_node.key = self.calculate_key(self.goal_node)
        heapq.heappush(self.U, (self.goal_node.key, id(self.goal_node), self.goal_node))

    def compute_shortest_path(self):
        """
        The main cost propagation loop of the algorithm
        Pops inconsistent nodes from the priority queue and updates their neighbors
        until the start position is consistent and the optimal path is found.
        """
        while self.U:
            k, _, u = self.U[0]
            curr_key = self.calculate_key(u)

            if k < curr_key:
                heapq.heappop(self.U)
                heapq.heappush(self.U, (curr_key, id(u), u))
                continue

            heapq.heappop(self.U)

            if u.g == u.rhs:
                continue

            if u.g > u.rhs:
                u.g = u.rhs
                for n in self.get_neighbors(u): self.update_node(n)
            else:
                u.g = math.inf
                self.update_node(u)
                for n in self.get_neighbors(u): self.update_node(n)

        self.is_optimal = True
        self.extract_path()

    def extract_path(self):
        """Edge sampling method to find the most optimal path """
        self.path = [(self.start_pos[0], self.start_pos[1])]
        curr_x, curr_y = self.start_pos
        visited = set()

        for _ in range(300):
            if math.hypot(curr_x - self.goal_node.x, curr_y - self.goal_node.y) < TOLERANCE:
                break

            best_p = None
            min_cost = math.inf

            # Find only the cells where the agent is physically located
            # TOLERANCE is used for cases where the agent is exactly on a cell boundary
            cells_to_check = []
            min_cx = int(max(0, math.floor(curr_x - TOLERANCE)))
            max_cx = int(min(self.cols - 1, math.floor(curr_x + TOLERANCE)))
            min_cy = int(max(0, math.floor(curr_y - TOLERANCE)))
            max_cy = int(min(self.rows - 1, math.floor(curr_y + TOLERANCE)))

            for cx in range(min_cx, max_cx + 1):
                for cy in range(min_cy, max_cy + 1):
                    cells_to_check.append((cx, cy))

            for cx, cy in cells_to_check:
                # Skip evaluating paths through obstacle cells
                c_cost = self.get_cell_cost(cx, cy)
                if c_cost == math.inf:
                    continue

                edges = [
                    ((cx, cy), (cx + 1, cy)), ((cx + 1, cy), (cx + 1, cy + 1)),
                    ((cx, cy + 1), (cx + 1, cy + 1)), ((cx, cy), (cx, cy + 1))
                ]

                for (x1, y1), (x2, y2) in edges:
                    n1, n2 = self.nodes[x1][y1], self.nodes[x2][y2]
                    if n1.g == math.inf and n2.g == math.inf:
                        continue

                    # Sample points along the edge to find the best path
                    STEPS = 20
                    for i in range(STEPS + 1):
                        px = x1 + (x2 - x1) * i / STEPS
                        py = y1 + (y2 - y1) * i / STEPS

                        dist = math.hypot(px - curr_x, py - curr_y)
                        if dist < 0.05:
                            continue

                        # Field D* interpolation equation
                        gp = n1.g + (n2.g - n1.g) * (i / STEPS)
                        cost = dist * c_cost + gp

                        if cost < min_cost:
                            min_cost = cost
                            best_p = (px, py)

            if not best_p:
                break

            # Prevention of infite loops (edge sampling)
            rounded = (round(best_p[0], 2), round(best_p[1], 2))
            if rounded in visited:
                break
            visited.add(rounded)

            curr_x, curr_y = best_p
            self.path.append((curr_x, curr_y))

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

        tk.Label(master, text="Spacja: Plan / Agent's move | Click: Add / Remove obstacle").pack()

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