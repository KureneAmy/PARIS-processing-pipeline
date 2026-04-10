#!/usr/bin/env python3
"""
intrxn_specificity.py - Multi-RNA Interaction Visualization Tool
Supports user-defined small RNA pair visualization

Usage:
    python intrxn_specificity.py <input_file> <output_dir> <config_json>

Example:
    python intrxn_specificity.py interactions.filtered ./plots rna_pairs.json
"""

import sys
import json
import re
import os
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from collections import defaultdict

class RNAInteractionVisualizer:
    """RNA Interaction Visualization Class"""
    
    def __init__(self, interactions_file, output_dir, rna_config):
        """
        Initialize the visualizer
        
        Args:
            interactions_file: Path to the interaction file
            output_dir: Output directory path
            rna_config: RNA configuration dictionary or JSON file path
        """
        self.interactions_file = interactions_file
        self.output_dir = output_dir
        self.rna_pairs = []
        self.rna_info = {}
        self.interactions = defaultdict(int)
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Load configuration
        if isinstance(rna_config, str):
            with open(rna_config, 'r') as f:
                config = json.load(f)
        else:
            config = rna_config
        
        self._parse_config(config)
    
    def _parse_config(self, config):
        """Parse RNA configuration"""
        if 'rna_pairs' in config:
            # Configuration contains list of RNA pairs
            for pair in config['rna_pairs']:
                rna1 = pair['rna1']
                rna2 = pair['rna2']
                size1 = pair.get('size1', 5070)  # Default to 28S size
                size2 = pair.get('size2', 13357)  # Default to 45S size
                
                self.rna_pairs.append({
                    'rna1': rna1,
                    'rna2': rna2,
                    'size1': size1,
                    'size2': size2,
                    'label1': pair.get('label1', rna1),
                    'label2': pair.get('label2', rna2)
                })
                
                self.rna_info[rna1] = {'size': size1}
                self.rna_info[rna2] = {'size': size2}
        else:
            raise ValueError("Config must contain 'rna_pairs' key")
        
        print(f"[INFO] Loaded {len(self.rna_pairs)} RNA pairs for visualization")
    
    def read_interactions(self):
        """Read interaction data from file"""
        print(f"[INFO] Reading interactions from: {self.interactions_file}")
        
        try:
            with open(self.interactions_file, 'r') as f:
                line_count = 0
                matched_count = 0
                
                for line in f:
                    line_count += 1
                    
                    if line.startswith("Group") or line.startswith("#"):
                        continue
                    
                    try:
                        parts = line.strip().split()
                        if len(parts) < 2:
                            continue
                        
                        intrxn_str = parts[1]
                        if "<=>Not found" in intrxn_str:
                            intrxn_str = intrxn_str.replace("<=>Not found", "<=>")
                        
                        interactions = intrxn_str.replace("<=>", " ").split()
                        
                        if len(interactions) < 2:
                            continue
                        
                        rna1_str, rna2_str = interactions[0], interactions[1]
                        
                        # Parse RNA names and coordinates
                        rna1_match = re.search(r'([A-Za-z0-9_.-]+)\|([+-]):(\d+)-(\d+)', rna1_str)
                        rna2_match = re.search(r'([A-Za-z0-9_.-]+)\|([+-]):(\d+)-(\d+)', rna2_str)
                        
                        if not (rna1_match and rna2_match):
                            continue
                        
                        rna1_name, strand1, pos1_s, pos1_e = rna1_match.groups()
                        rna2_name, strand2, pos2_s, pos2_e = rna2_match.groups()
                        
                        pos1_s, pos1_e = int(pos1_s), int(pos1_e)
                        pos2_s, pos2_e = int(pos2_s), int(pos2_e)
                        
                        # For each RNA pair, record the interactions
                        for pair in self.rna_pairs:
                            if ((pair['rna1'] in rna1_name and pair['rna2'] in rna2_name) or 
                                (pair['rna1'] in rna2_name and pair['rna2'] in rna1_name)):
                                
                                # Use midpoint of coordinates
                                if pair['rna1'] in rna1_name:
                                    key = (pos1_s + pos1_e) // 2
                                    self.interactions[f"{pair['rna1']}_{pair['rna2']}"] += 1
                                else:
                                    key = (pos2_s + pos2_e) // 2
                                    self.interactions[f"{pair['rna1']}_{pair['rna2']}"] += 1
                                
                                matched_count += 1
                        
                    except (IndexError, ValueError, AttributeError):
                        continue
            
            print(f"[INFO] Processed {line_count} lines, matched {matched_count} interactions")
            
        except FileNotFoundError:
            print(f"[ERROR] File not found: {self.interactions_file}")
            raise
    
    def plot_coverage_heatmap(self):
        """Generate interaction coverage heatmap"""
        print(f"[INFO] Generating interaction heatmaps...")
        
        for pair in self.rna_pairs:
            rna1 = pair['rna1']
            rna2 = pair['rna2']
            size1 = pair['size1']
            size2 = pair['size2']
            label1 = pair['label1']
            label2 = pair['label2']
            
            # Initialize 2D matrix to store interaction information
            interaction_matrix = [[0 for _ in range(size2)] for _ in range(size1)]
            
            try:
                with open(self.interactions_file, 'r') as f:
                    for line in f:
                        if line.startswith("Group") or line.startswith("#"):
                            continue
                        
                        try:
                            parts = line.strip().split()
                            if len(parts) < 2:
                                continue
                            
                            intrxn_str = parts[1]
                            if "<=>Not found" in intrxn_str:
                                intrxn_str = intrxn_str.replace("<=>Not found", "<=>")
                            
                            interactions = intrxn_str.replace("<=>", " ").split()
                            
                            if len(interactions) < 2:
                                continue
                            
                            rna1_str, rna2_str = interactions[0], interactions[1]
                            
                            rna1_match = re.search(r'([A-Za-z0-9_.-]+)\|([+-]):(\d+)-(\d+)', rna1_str)
                            rna2_match = re.search(r'([A-Za-z0-9_.-]+)\|([+-]):(\d+)-(\d+)', rna2_str)
                            
                            if not (rna1_match and rna2_match):
                                continue
                            
                            rna1_name, strand1, pos1_s, pos1_e = rna1_match.groups()
                            rna2_name, strand2, pos2_s, pos2_e = rna2_match.groups()
                            
                            pos1_s, pos1_e = int(pos1_s), int(pos1_e)
                            pos2_s, pos2_e = int(pos2_s), int(pos2_e)
                            
                            # Check if this matches the current RNA pair
                            if ((rna1 in rna1_name and rna2 in rna2_name) or 
                                (rna1 in rna2_name and rna2 in rna1_name)):
                                
                                if rna1 in rna1_name:
                                    pos1_mid = (pos1_s + pos1_e) // 2
                                    pos2_mid = (pos2_s + pos2_e) // 2
                                else:
                                    pos1_mid = (pos2_s + pos2_e) // 2
                                    pos2_mid = (pos1_s + pos1_e) // 2
                                
                                # Prevent out of bounds
                                if pos1_mid < size1 and pos2_mid < size2:
                                    interaction_matrix[pos1_mid][pos2_mid] += 1
                        
                        except (IndexError, ValueError, AttributeError):
                            continue
            
            except Exception as e:
                print(f"[WARNING] Error processing heatmap for {rna1} vs {rna2}: {e}")
                continue
            
            # Plot heatmap
            self._plot_single_heatmap(interaction_matrix, rna1, rna2, label1, label2, size1, size2)
    
    def _plot_single_heatmap(self, matrix, rna1, rna2, label1, label2, size1, size2):
        """Plot a single heatmap"""
        import numpy as np
        
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Convert to numpy array for plotting
        data = np.array(matrix)
        
        # Use log scale to highlight low-frequency interactions
        data_log = np.log2(data + 1)
        
        im = ax.imshow(data_log, cmap='YlOrRd', aspect='auto', interpolation='nearest')
        
        ax.set_xlabel(f'{label2} Position', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'{label1} Position', fontsize=12, fontweight='bold')
        ax.set_title(f'{label1} vs {label2} Interaction Heatmap', fontsize=14, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('log2(Interaction Count)', fontsize=10)
        
        # Set ticks
        if size1 > 1000:
            ax.set_xticks(np.linspace(0, size2-1, 5))
            ax.set_yticks(np.linspace(0, size1-1, 5))
            ax.set_xticklabels([f'{int(x*size2/5)}' for x in range(5)])
            ax.set_yticklabels([f'{int(y*size1/5)}' for y in range(5)])
        
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, f'{rna1}_vs_{rna2}_heatmap.pdf')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"[INFO] Saved heatmap: {output_file}")
    
    def plot_coverage_bars(self):
        """Generate coverage bar plots"""
        print(f"[INFO] Generating coverage bar plots...")
        
        for pair in self.rna_pairs:
            rna1 = pair['rna1']
            rna2 = pair['rna2']
            size1 = pair['size1']
            size2 = pair['size2']
            label1 = pair['label1']
            label2 = pair['label2']
            
            # Initialize coverage arrays
            coverage1 = [0] * size1
            coverage2 = [0] * size2
            
            try:
                with open(self.interactions_file, 'r') as f:
                    for line in f:
                        if line.startswith("Group") or line.startswith("#"):
                            continue
                        
                        try:
                            parts = line.strip().split()
                            if len(parts) < 2:
                                continue
                            
                            intrxn_str = parts[1]
                            if "<=>Not found" in intrxn_str:
                                intrxn_str = intrxn_str.replace("<=>Not found", "<=>")
                            
                            interactions = intrxn_str.replace("<=>", " ").split()
                            
                            if len(interactions) < 2:
                                continue
                            
                            rna1_str, rna2_str = interactions[0], interactions[1]
                            
                            rna1_match = re.search(r'([A-Za-z0-9_.-]+)\|([+-]):(\d+)-(\d+)', rna1_str)
                            rna2_match = re.search(r'([A-Za-z0-9_.-]+)\|([+-]):(\d+)-(\d+)', rna2_str)
                            
                            if not (rna1_match and rna2_match):
                                continue
                            
                            rna1_name, strand1, pos1_s, pos1_e = rna1_match.groups()
                            rna2_name, strand2, pos2_s, pos2_e = rna2_match.groups()
                            
                            pos1_s, pos1_e = int(pos1_s), int(pos1_e)
                            pos2_s, pos2_e = int(pos2_s), int(pos2_e)
                            
                            # Check if this matches the current RNA pair
                            if ((rna1 in rna1_name and rna2 in rna2_name) or 
                                (rna1 in rna2_name and rna2 in rna1_name)):
                                
                                if rna1 in rna1_name:
                                    # Coverage for RNA1
                                    for i in range(pos1_s, min(pos1_e, size1)):
                                        if i < size1:
                                            coverage1[i] += 1
                                    # Coverage for RNA2
                                    for i in range(pos2_s, min(pos2_e, size2)):
                                        if i < size2:
                                            coverage2[i] += 1
                                else:
                                    # Swapped case
                                    for i in range(pos2_s, min(pos2_e, size1)):
                                        if i < size1:
                                            coverage1[i] += 1
                                    for i in range(pos1_s, min(pos1_e, size2)):
                                        if i < size2:
                                            coverage2[i] += 1
                        
                        except (IndexError, ValueError, AttributeError):
                            continue
            
            except Exception as e:
                print(f"[WARNING] Error processing coverage for {rna1} vs {rna2}: {e}")
                continue
            
            # Plot bar charts
            self._plot_single_coverage(coverage1, label1, rna1, rna2, 'rna1')
            self._plot_single_coverage(coverage2, label2, rna1, rna2, 'rna2')
    
    def _plot_single_coverage(self, coverage, label, rna1, rna2, rna_type):
        """Plot coverage for a single RNA"""
        fig, ax = plt.subplots(figsize=(14, 4))
        
        ax.bar(range(len(coverage)), coverage, color='steelblue', edgecolor='none', width=1)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        
        ax.set_xlim(0, len(coverage))
        ax.set_xlabel(f'{label} Position', fontsize=12, fontweight='bold')
        ax.set_ylabel('Interaction Coverage', fontsize=12, fontweight='bold')
        
        if rna_type == 'rna1':
            title = f'{rna1} Coverage in Interactions with {rna2}'
        else:
            title = f'{rna2} Coverage in Interactions with {rna1}'
        
        ax.set_title(title, fontsize=13, fontweight='bold')
        
        plt.tight_layout()
        
        output_file = os.path.join(self.output_dir, f'{rna1}_vs_{rna2}_{rna_type}_coverage.pdf')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"[INFO] Saved coverage plot: {output_file}")
    
    def generate_summary_report(self):
        """Generate summary report"""
        report_file = os.path.join(self.output_dir, 'visualization_summary.txt')
        
        with open(report_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("RNA Interaction Visualization Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            f.write(f"Input file: {self.interactions_file}\n")
            f.write(f"Output directory: {self.output_dir}\n\n")
            
            f.write("RNA Pairs Analyzed:\n")
            f.write("-" * 60 + "\n")
            
            for i, pair in enumerate(self.rna_pairs, 1):
                f.write(f"{i}. {pair['label1']} (size: {pair['size1']}) vs {pair['label2']} (size: {pair['size2']})\n")
                f.write(f"   RNA identifiers: {pair['rna1']} <-> {pair['rna2']}\n")
                f.write(f"   Generated plots:\n")
                f.write(f"     - {pair['rna1']}_vs_{pair['rna2']}_heatmap.pdf\n")
                f.write(f"     - {pair['rna1']}_vs_{pair['rna2']}_rna1_coverage.pdf\n")
                f.write(f"     - {pair['rna1']}_vs_{pair['rna2']}_rna2_coverage.pdf\n")
                f.write("\n")
            
            f.write("=" * 60 + "\n")
            f.write("Visualization complete!\n")
        
        print(f"[INFO] Summary report saved: {report_file}")

def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    config_file = sys.argv[3]
    
    try:
        # Initialize visualizer
        visualizer = RNAInteractionVisualizer(input_file, output_dir, config_file)
        
        # Read interaction data
        visualizer.read_interactions()
        
        # Generate visualizations
        visualizer.plot_coverage_bars()
        visualizer.plot_coverage_heatmap()
        
        # Generate summary report
        visualizer.generate_summary_report()
        
        print("[SUCCESS] Visualization completed!")
        
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()