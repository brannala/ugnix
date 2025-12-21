#!/usr/bin/env python3
"""
Coalescent Tree Viewer for coalsim output

Displays phylograms from coalsim with:
- Trees proportional to generations (branch lengths)
- Chromosome diagram showing segment location in cM
- Scrollable interface for navigating between trees
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import re
import sys
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Optional

# Use Agg backend for matplotlib when embedded in tkinter
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle


@dataclass
class TreeNode:
    """Represents a node in the phylogenetic tree"""
    name: Optional[str] = None
    branch_length: float = 0.0
    children: List['TreeNode'] = None
    x: float = 0.0  # horizontal position for plotting
    y: float = 0.0  # vertical position (time/height)

    def __post_init__(self):
        if self.children is None:
            self.children = []

    def is_leaf(self) -> bool:
        return len(self.children) == 0


def parse_newick(newick: str) -> TreeNode:
    """Parse a Newick format string into a tree structure"""
    newick = newick.strip()
    if newick.endswith(';'):
        newick = newick[:-1]

    def parse_node(s: str, pos: int) -> Tuple[TreeNode, int]:
        node = TreeNode()

        if s[pos] == '(':
            # Internal node with children
            pos += 1  # skip '('
            while True:
                child, pos = parse_node(s, pos)
                node.children.append(child)
                if s[pos] == ',':
                    pos += 1
                elif s[pos] == ')':
                    pos += 1
                    break

        # Parse name and branch length
        name_length = ""
        while pos < len(s) and s[pos] not in '(),;':
            name_length += s[pos]
            pos += 1

        if ':' in name_length:
            parts = name_length.split(':')
            node.name = parts[0] if parts[0] else None
            try:
                node.branch_length = float(parts[1]) if parts[1] else 0.0
            except ValueError:
                node.branch_length = 0.0
        else:
            node.name = name_length if name_length else None

        return node, pos

    root, _ = parse_node(newick, 0)
    return root


def calculate_tree_height(node: TreeNode, current_height: float = 0.0) -> float:
    """Calculate height from root, then convert so tips are at y=0"""
    node.y = current_height
    if node.is_leaf():
        return current_height

    max_height = current_height
    for child in node.children:
        child_height = calculate_tree_height(child, current_height + child.branch_length)
        max_height = max(max_height, child_height)

    return max_height


def convert_to_time_from_present(node: TreeNode, max_height: float):
    """Convert y coordinates so tips are at 0 (present) and root is at TMRCA"""
    node.y = max_height - node.y
    for child in node.children:
        convert_to_time_from_present(child, max_height)


def assign_x_positions(node: TreeNode, leaf_count: List[int] = None) -> float:
    """Assign x positions to nodes, returning the x position of the node"""
    if leaf_count is None:
        leaf_count = [0]

    if node.is_leaf():
        node.x = leaf_count[0]
        leaf_count[0] += 1
        return node.x

    child_positions = []
    for child in node.children:
        child_positions.append(assign_x_positions(child, leaf_count))

    node.x = sum(child_positions) / len(child_positions)
    return node.x


def count_leaves(node: TreeNode) -> int:
    """Count the number of leaves in the tree"""
    if node.is_leaf():
        return 1
    return sum(count_leaves(child) for child in node.children)


def draw_tree(ax, node: TreeNode, axis_height: float):
    """Draw the tree as a proper phylogram.

    Tips are at y=0 (present), root at y=TMRCA (past).
    Standard phylogram: horizontal lines at parent height, vertical lines to children.
    axis_height is used for consistent label positioning.
    """
    for child in node.children:
        # Proper phylogram shape:
        # 1. Horizontal line at parent's height from parent.x to child.x
        ax.plot([node.x, child.x], [node.y, node.y], 'k-', linewidth=1.5)
        # 2. Vertical line at child's x from parent's height down to child's height
        ax.plot([child.x, child.x], [node.y, child.y], 'k-', linewidth=1.5)

        # Recursively draw children
        draw_tree(ax, child, axis_height)

    # Draw leaf markers and labels
    if node.is_leaf():
        label = node.name if node.name else ""
        # Circle at tip
        ax.plot(node.x, node.y, 'ko', markersize=6)
        # Label below tip
        ax.text(node.x, -axis_height * 0.03, label,
                ha='center', va='top', fontsize=10, fontweight='bold')


class CoalescentTreeViewer:
    """Main application class for the tree viewer"""

    def __init__(self, root: tk.Tk, trees: List[str], intervals: List[Tuple[float, float]],
                 chromosome_length_cm: float = 100.0):
        self.root = root
        self.root.title("Coalescent Tree Viewer")
        self.root.geometry("1200x900")

        self.trees = trees
        self.intervals = intervals
        self.chromosome_length_cm = chromosome_length_cm
        self.current_index = 0
        self.fixed_time_axis = tk.BooleanVar(value=False)

        # Pre-calculate max TMRCA across all trees for fixed axis mode
        self.global_max_height = self._calculate_global_max_height()

        self.setup_ui()
        self.display_tree(0)

    def _calculate_global_max_height(self) -> float:
        """Calculate the maximum TMRCA across all trees"""
        max_height = 0.0
        for tree_str in self.trees:
            try:
                tree = parse_newick(tree_str)
                height = calculate_tree_height(tree)
                max_height = max(max_height, height)
            except:
                pass
        return max_height

    def setup_ui(self):
        """Set up the user interface"""
        # Menu bar
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # View menu
        view_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="View", menu=view_menu)
        view_menu.add_checkbutton(
            label="Fixed Time Axis",
            variable=self.fixed_time_axis,
            command=self._on_axis_mode_change
        )

        # Main container
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Top frame for navigation
        nav_frame = ttk.Frame(main_frame)
        nav_frame.pack(fill=tk.X, pady=(0, 10))

        # Navigation buttons
        self.prev_btn = ttk.Button(nav_frame, text="< Previous", command=self.prev_tree)
        self.prev_btn.pack(side=tk.LEFT, padx=5)

        self.next_btn = ttk.Button(nav_frame, text="Next >", command=self.next_tree)
        self.next_btn.pack(side=tk.LEFT, padx=5)

        # Tree counter label
        self.counter_label = ttk.Label(nav_frame, text="", font=('Helvetica', 12))
        self.counter_label.pack(side=tk.LEFT, padx=20)

        # Slider for quick navigation
        self.slider = ttk.Scale(nav_frame, from_=0, to=max(0, len(self.trees)-1),
                                orient=tk.HORIZONTAL, command=self.on_slider)
        self.slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)

        # Jump to specific tree
        ttk.Label(nav_frame, text="Go to:").pack(side=tk.LEFT, padx=(10, 2))
        self.jump_entry = ttk.Entry(nav_frame, width=6)
        self.jump_entry.pack(side=tk.LEFT, padx=2)
        self.jump_entry.bind('<Return>', self.jump_to_tree)
        ttk.Button(nav_frame, text="Go", command=self.jump_to_tree).pack(side=tk.LEFT, padx=2)

        # Create figure with two subplots
        self.fig = Figure(figsize=(12, 9), dpi=100, facecolor='white')

        # Chromosome diagram at top (smaller)
        self.ax_chrom = self.fig.add_axes([0.1, 0.85, 0.8, 0.10])

        # Tree diagram (larger, below chromosome)
        self.ax_tree = self.fig.add_axes([0.1, 0.08, 0.8, 0.72])

        # Canvas for matplotlib
        self.canvas = FigureCanvasTkAgg(self.fig, master=main_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Toolbar for zoom/pan
        toolbar_frame = ttk.Frame(main_frame)
        toolbar_frame.pack(fill=tk.X)
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()

        # Keyboard bindings
        self.root.bind('<Left>', lambda e: self.prev_tree())
        self.root.bind('<Right>', lambda e: self.next_tree())
        self.root.bind('<Home>', lambda e: self.display_tree(0))
        self.root.bind('<End>', lambda e: self.display_tree(len(self.trees) - 1))

    def _on_axis_mode_change(self):
        """Called when time axis mode is toggled"""
        self.display_tree(self.current_index)

    def draw_chromosome(self, interval: Tuple[float, float]):
        """Draw simple chromosome diagram with highlighted segment"""
        self.ax_chrom.clear()

        # Simple rectangle, 3/4 width, centered
        chrom_width = 0.75  # 75% of axis width
        chrom_left = (1.0 - chrom_width) / 2  # centered
        chrom_height = 0.4
        chrom_y = 0.3

        # Convert interval from (0,1) to rectangle coordinates
        start_frac = interval[0]
        end_frac = interval[1]

        # Draw chromosome background (gray rectangle)
        chrom_rect = Rectangle((chrom_left, chrom_y),
                                chrom_width, chrom_height,
                                facecolor='#D0D0D0', edgecolor='#404040',
                                linewidth=2)
        self.ax_chrom.add_patch(chrom_rect)

        # Draw highlighted segment (red) within the chromosome
        seg_left = chrom_left + start_frac * chrom_width
        seg_width = (end_frac - start_frac) * chrom_width
        if seg_width > 0:
            highlight = Rectangle((seg_left, chrom_y),
                                   seg_width, chrom_height,
                                   facecolor='#E74C3C', edgecolor='#C0392B',
                                   linewidth=1.5)
            self.ax_chrom.add_patch(highlight)

        # Labels at ends showing cM
        self.ax_chrom.text(chrom_left, chrom_y - 0.08, '0',
                          ha='center', va='top', fontsize=9)
        self.ax_chrom.text(chrom_left + chrom_width, chrom_y - 0.08,
                          f'{self.chromosome_length_cm:.0f} cM',
                          ha='center', va='top', fontsize=9)

        # Segment position label above in kb (1 cM = 1 Mb = 1000 kb)
        start_kb = interval[0] * self.chromosome_length_cm * 1000
        end_kb = interval[1] * self.chromosome_length_cm * 1000
        seg_center = seg_left + seg_width / 2
        self.ax_chrom.text(seg_center, chrom_y + chrom_height + 0.08,
                          f'[{start_kb:.1f} kb, {end_kb:.1f} kb]',
                          ha='center', va='bottom', fontsize=10, fontweight='bold')

        # Set axis properties
        self.ax_chrom.set_xlim(0, 1)
        self.ax_chrom.set_ylim(0, 1)
        self.ax_chrom.set_aspect('auto')
        self.ax_chrom.axis('off')

    def display_tree(self, index: int):
        """Display tree at the given index"""
        if not self.trees:
            return

        index = max(0, min(index, len(self.trees) - 1))
        self.current_index = index

        # Update navigation state
        self.prev_btn.state(['!disabled'] if index > 0 else ['disabled'])
        self.next_btn.state(['!disabled'] if index < len(self.trees) - 1 else ['disabled'])
        self.counter_label.config(text=f"Tree {index + 1} of {len(self.trees)}")
        self.slider.set(index)

        # Draw chromosome diagram
        if self.intervals:
            self.draw_chromosome(self.intervals[index])

        # Parse and draw tree
        self.ax_tree.clear()

        try:
            tree = parse_newick(self.trees[index])
            tree_height = calculate_tree_height(tree)
            # Convert so tips are at y=0 (present) and root at y=TMRCA (past)
            convert_to_time_from_present(tree, tree_height)
            assign_x_positions(tree)
            num_leaves = count_leaves(tree)

            # Use fixed or variable axis scale
            if self.fixed_time_axis.get():
                axis_height = self.global_max_height
            else:
                axis_height = tree_height

            # Draw the tree
            draw_tree(self.ax_tree, tree, axis_height)

            # Set axis properties - tips at y=0, root at y=axis_height
            self.ax_tree.set_xlim(-0.5, num_leaves - 0.5)
            self.ax_tree.set_ylim(-axis_height * 0.08, axis_height * 1.02)

            # Y-axis: time from present (0 at bottom) to past (TMRCA at top)
            self.ax_tree.set_ylabel('Time (coalescent units)', fontsize=12, fontweight='bold')

            # Remove x-axis ticks (leaves are labeled)
            self.ax_tree.set_xticks([])
            self.ax_tree.set_xlabel('Samples', fontsize=12, fontweight='bold')

            # Add grid for time reference
            self.ax_tree.yaxis.grid(True, linestyle='--', alpha=0.3)
            self.ax_tree.set_axisbelow(True)

            # Title
            if self.intervals:
                start_cm = self.intervals[index][0] * self.chromosome_length_cm
                end_cm = self.intervals[index][1] * self.chromosome_length_cm
                self.ax_tree.set_title(
                    f'Gene Tree for Segment {start_cm:.2f} - {end_cm:.2f} cM (TMRCA = {tree_height:.3f})',
                    fontsize=12, fontweight='bold', pad=10)
            else:
                self.ax_tree.set_title(f'Gene Tree (TMRCA = {tree_height:.3f})',
                                       fontsize=12, fontweight='bold', pad=10)

        except Exception as e:
            self.ax_tree.text(0.5, 0.5, f'Error parsing tree:\n{str(e)}\n\nTree string:\n{self.trees[index][:100]}...',
                             ha='center', va='center', transform=self.ax_tree.transAxes,
                             fontsize=10, color='red')

        self.canvas.draw()

    def prev_tree(self):
        """Show previous tree"""
        if self.current_index > 0:
            self.display_tree(self.current_index - 1)

    def next_tree(self):
        """Show next tree"""
        if self.current_index < len(self.trees) - 1:
            self.display_tree(self.current_index + 1)

    def on_slider(self, value):
        """Handle slider movement"""
        index = int(float(value))
        if index != self.current_index:
            self.display_tree(index)

    def jump_to_tree(self, event=None):
        """Jump to a specific tree number"""
        try:
            index = int(self.jump_entry.get()) - 1  # Convert to 0-indexed
            if 0 <= index < len(self.trees):
                self.display_tree(index)
            else:
                messagebox.showwarning("Invalid Index",
                                      f"Please enter a number between 1 and {len(self.trees)}")
        except ValueError:
            messagebox.showwarning("Invalid Input", "Please enter a valid number")


def load_trees(tree_file: str) -> List[str]:
    """Load trees from file"""
    trees = []
    with open(tree_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                trees.append(line)
    return trees


def load_intervals(interval_file: str) -> List[Tuple[float, float]]:
    """Load intervals from file"""
    intervals = []
    pattern = re.compile(r'\[([0-9.]+),([0-9.]+)\]')
    with open(interval_file, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                intervals.append((float(match.group(1)), float(match.group(2))))
    return intervals


def main():
    parser = argparse.ArgumentParser(
        description='View coalescent trees from coalsim output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                           # Default: 100 cM chromosome
  %(prog)s -c 5                      # For coalsim -r0.05 (5 cM)
  %(prog)s -c 50 -t mytrees.txt      # Custom files, 50 cM

Chromosome length calculation:
  coalsim uses recombination rate in Morgans (-r option)
  Chromosome length (cM) = recRate * 100

  Examples:
    coalsim -r0.01  ->  use -c 1    (1 cM)
    coalsim -r0.05  ->  use -c 5    (5 cM)
    coalsim -r0.1   ->  use -c 10   (10 cM)
    coalsim -r1.0   ->  use -c 100  (100 cM, default)

The program expects:
  - Tree file: Newick format trees, one per line (from -g f)
  - Interval file: [start,end] format, one per line

Keyboard shortcuts:
  Left/Right arrows: Previous/Next tree
  Home/End: First/Last tree

Run coalsim with -g f option to generate these files.
        """)

    parser.add_argument('-t', '--trees', default='trees.txt',
                        help='File containing Newick trees (default: trees.txt)')
    parser.add_argument('-i', '--intervals', default='mrcaintv.txt',
                        help='File containing MRCA intervals (default: mrcaintv.txt)')
    parser.add_argument('-c', '--chromosome-length', type=float, default=None,
                        help='Chromosome length in cM (default: 100.0, or from -r)')
    parser.add_argument('-r', '--recrate', type=float, default=None,
                        help='Recombination rate from coalsim (computes cM = rate*100)')

    args = parser.parse_args()

    # Determine chromosome length
    if args.recrate is not None:
        chromosome_length = args.recrate * 100  # Convert Morgans to cM
    elif args.chromosome_length is not None:
        chromosome_length = args.chromosome_length
    else:
        chromosome_length = 100.0  # Default

    # Load data
    try:
        trees = load_trees(args.trees)
    except FileNotFoundError:
        print(f"Error: Tree file '{args.trees}' not found.")
        print("Run coalsim with -g f option to generate tree output.")
        sys.exit(1)

    try:
        intervals = load_intervals(args.intervals)
    except FileNotFoundError:
        print(f"Warning: Interval file '{args.intervals}' not found. Proceeding without intervals.")
        intervals = [(0, 1)] * len(trees)  # Default to full chromosome

    if not trees:
        print("Error: No trees found in tree file.")
        sys.exit(1)

    if len(intervals) != len(trees):
        print(f"Warning: Number of intervals ({len(intervals)}) doesn't match number of trees ({len(trees)}).")
        # Pad or truncate intervals
        if len(intervals) < len(trees):
            intervals.extend([(0, 1)] * (len(trees) - len(intervals)))
        else:
            intervals = intervals[:len(trees)]

    print(f"Loaded {len(trees)} trees from {args.trees}")
    print(f"Loaded {len(intervals)} intervals from {args.intervals}")
    print(f"Chromosome length: {chromosome_length} cM")

    # Create and run application
    root = tk.Tk()
    app = CoalescentTreeViewer(root, trees, intervals, chromosome_length)
    root.mainloop()


if __name__ == '__main__':
    main()
