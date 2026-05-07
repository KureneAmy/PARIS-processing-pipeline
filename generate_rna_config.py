#!/usr/bin/env python3
"""
generate_rna_config.py - get config JSON from config.yaml
"""

import json
import sys
import yaml

def main():
    if len(sys.argv) < 3:
        print("Usage: python generate_rna_config.py <config.yaml> <output.json>")
        sys.exit(1)
    
    config_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    rna_viz = config.get('rna_intrxn_visualization', {})
    
    output_config = {
        'rna_pairs': rna_viz.get('rna_pairs', [])
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_config, f, indent=2)
    
    print(f"[INFO] Generated RNA config: {output_file}")

if __name__ == "__main__":
    main()