#!/usr/bin/env python3
"""
Create IGV batch script to generate alignment screenshots
Usage: python create_igv_batch.py [config.yaml]
"""

import sys
import os
import yaml

def create_igv_batch_script(regions_file, output_dir, genome_fasta, bam_file, 
                            problematic_file, window_padding=50, max_panel_height=2000):
    """Generate IGV batch script from BED file"""
    
    batch_commands = []
    
    # Initial setup
    batch_commands.append("new")
    batch_commands.append(f"genome {genome_fasta}")
    batch_commands.append(f"load {bam_file}")
    batch_commands.append(f"load {problematic_file}")
    batch_commands.append("snapshotDirectory " + output_dir)
    batch_commands.append(f"maxPanelHeight {max_panel_height}")
    
    # Read positions from BED file
    with open(regions_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            
            chrom, start, end, region_id = parts[0], int(parts[1]), int(parts[2]), parts[3]
            pos_1based = start + 1
            
            # Add context window
            window_start = max(1, pos_1based - window_padding)
            window_end = end + window_padding
            
            region = f"{chrom}:{window_start}-{window_end}"
            snapshot_name = f"{region_id}_{chrom}.png"
            
            batch_commands.append(f"goto {region}")
            batch_commands.append("sort base")
            batch_commands.append(f"snapshot {snapshot_name}")
    
    batch_commands.append("exit")
    
    return '\n'.join(batch_commands)

def load_config(config_file):
    """Load configuration from YAML file"""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    # Get config file from command line or use default
    config_file = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
    
    if not os.path.exists(config_file):
        print(f"Error: Config file '{config_file}' not found")
        print(f"Usage: {sys.argv[0]} [config.yaml]")
        sys.exit(1)
    
    # Load configuration
    config = load_config(config_file)
    paths = config['paths']
    settings = config['igv_settings']
    
    # Create output directory if it doesn't exist
    os.makedirs(paths['output_dir'], exist_ok=True)
    
    # Generate batch script
    batch_script = create_igv_batch_script(
        regions_file=paths['regions_file'],
        output_dir=paths['output_dir'],
        genome_fasta=paths['genome_fasta'],
        bam_file=paths['bam_file'],
        problematic_file=paths['problematic_file'],
        window_padding=settings['window_padding'],
        max_panel_height=settings['max_panel_height']
    )
    
    # Write batch script
    batch_file = paths['script_file']
    with open(batch_file, 'w') as f:
        f.write(batch_script)
    
    print(f"IGV batch script created: {batch_file}")
    print(f"Using config: {config_file}")
    print(f"\nTo run IGV in batch mode:")
    print(f"  igv.sh -b {batch_file}")
    print(f"\nScreenshots will be saved to: {paths['output_dir']}/")

if __name__ == "__main__":
    main()