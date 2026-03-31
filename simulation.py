import math
import tkinter as tk

from field_d_star import Cell, FieldDStar, Node, OBSTACLE_VAL


class FieldDStarApp:
    """Tkinter UI for editing a cost grid and stepping the Field D* planner."""
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

        self.controls = tk.Frame(master)
        self.controls.pack(fill='x', padx=8, pady=6)

        self.brush_mode = tk.StringVar(value='cost')
        tk.Radiobutton(self.controls, text='Paint cost', variable=self.brush_mode, value='cost').pack(side='left')
        tk.Radiobutton(self.controls, text='Obstacle (254)', variable=self.brush_mode, value='obstacle').pack(side='left')
        tk.Radiobutton(self.controls, text='Erase (0)', variable=self.brush_mode, value='erase').pack(side='left')

        self.cost_scale = tk.Scale(
            self.controls,
            from_=0,
            to=253,
            orient='horizontal',
            length=260,
            label='Cell cost (0-253)',
        )
        self.cost_scale.set(80)
        self.cost_scale.pack(side='left', padx=10)

        self.cost_label = tk.Label(self.controls, text='Brush value: 80')
        self.cost_label.pack(side='left', padx=6)
        self.cost_scale.configure(command=self.update_brush_label)

        master.bind('<space>', lambda event: self.handle_space_press())
        self.canvas.bind('<Button-1>', self.paint_cell)
        self.canvas.bind('<B1-Motion>', self.paint_cell)
        self.canvas.bind('<Button-3>', self.erase_cell)
        self.canvas.bind('<B3-Motion>', self.erase_cell)

        tk.Label(
            master,
            text='Space: plan / move agent | Left click-drag: paint | Right click-drag: erase to 0',
        ).pack()

    def draw_grid(self):
        """Creates the static grid widgets and goal marker."""
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

    def update_brush_label(self, _value=None):
        """Refreshes the label showing the currently selected brush value."""
        self.cost_label.config(text=f'Brush value: {int(self.cost_scale.get())}')

    def cell_display_color(self, cost_value):
        """Maps cell cost to a visible heat-style fill color."""
        if cost_value >= OBSTACLE_VAL:
            return '#654321'
        # white -> orange/red as cost rises
        ratio = max(0.0, min(1.0, cost_value / 253.0))
        red = 255
        green = int(255 - 135 * ratio)
        blue = int(255 - 255 * ratio)
        return f'#{red:02x}{green:02x}{blue:02x}'

    def apply_cell_value(self, cx, cy, value):
        """Writes a cost into one cell and invalidates affected planner nodes."""
        if not (0 <= cx < self.COLS and 0 <= cy < self.ROWS):
            return

        value = int(max(0, min(OBSTACLE_VAL, value)))

        # Prevent placing an obstacle on the agent footprint cell.
        ax, ay = self.dstar.start_pos
        if cx <= ax <= cx + 1 and cy <= ay <= cy + 1 and value >= OBSTACLE_VAL:
            return

        cell = self.cells[cx][cy]
        if int(cell.cost_value) == value:
            return

        cell.cost_value = value

        # If cost of cell has changed, replan its 4 corners.
        for nx in (cx, cx + 1):
            for ny in (cy, cy + 1):
                self.dstar.update_node(self.nodes[nx][ny])

        self.dstar.is_optimal = False
        self.update_visualization()

    def brush_value(self):
        """Returns the integer cost represented by the active brush mode."""
        mode = self.brush_mode.get()
        if mode == 'obstacle':
            return OBSTACLE_VAL
        if mode == 'erase':
            return 0
        return int(self.cost_scale.get())

    def paint_cell(self, event):
        """Paints the cell under the cursor using the active brush."""
        cx, cy = event.x // self.CELL_SIZE, event.y // self.CELL_SIZE
        self.apply_cell_value(cx, cy, self.brush_value())

    def erase_cell(self, event):
        """Resets the cell under the cursor to free-space cost."""
        cx, cy = event.x // self.CELL_SIZE, event.y // self.CELL_SIZE
        self.apply_cell_value(cx, cy, 0)

    def move_agent(self):
        """Moves the agent by one extracted waypoint and refreshes the path."""
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
        """Plans once, then advances the agent on subsequent space presses."""
        if not self.dstar.is_optimal:
            self.dstar.compute_shortest_path()
        elif self.dstar.is_optimal:
            self.move_agent()

        self.update_visualization()

    def update_visualization(self):
        """Redraws cell costs, extracted path and current agent position."""
        for x in range(self.COLS):
            for y in range(self.ROWS):
                c = self.cells[x][y]
                self.canvas.itemconfig(c.rect_id, fill=self.cell_display_color(int(c.cost_value)))

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
    """Launches the standalone Field D* grid editor and simulation window."""
    ROOT = tk.Tk()
    APP = FieldDStarApp(ROOT, cols=20, rows=15, cell_size=50)
    ROOT.mainloop()
