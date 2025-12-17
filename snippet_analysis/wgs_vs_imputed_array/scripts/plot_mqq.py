#!/usr/bin/env python3
"""
Plot Manhattan and QQ plots from PLINK2 association results using gwaslab.
"""

import argparse
import logging
import gwaslab as gl
import matplotlib.pyplot as plt

# Suppress matplotlib warnings
logging.getLogger("matplotlib").setLevel(logging.ERROR)
plt.style.use('default')

# Set font to sans-serif to avoid Arial warnings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Bitstream Vera Sans', 'Computer Modern Sans Serif', 'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Arial', 'Helvetica', 'Avant Garde', 'sans-serif']


def main():
    parser = argparse.ArgumentParser(
        description='Generate Manhattan and QQ plots from PLINK2 GWAS results'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input PLINK2 .glm.logistic file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output plot filename (e.g., results.mqq.png)'
    )
    parser.add_argument(
        '--build',
        default='38',
        choices=['37', '38'],
        help='Genome build version (default: 38)'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=400,
        help='DPI for output image (default: 400)'
    )
    parser.add_argument(
        '--title',
        default='',
        help='Optional title for the plots'
    )
    parser.add_argument(
        '--sig-level',
        type=float,
        default=5e-8,
        help='Significance level threshold (default: 5e-8)'
    )
    
    args = parser.parse_args()
    
    print(f"Loading GWAS summary statistics from: {args.input}")
    
    # Load PLINK2 results
    sumstats = gl.Sumstats(
        args.input,
        fmt="plink2",
        build=args.build,
        ea='A1',
        nea='OMITTED',
        OR_95L='L95',
        OR_95U='U95'
    )
    
    print(f"Generating Manhattan and QQ plots...")
    print(f"Significance level: {args.sig_level}")
    
    # Plot MQQ (Manhattan + QQ)
    sumstats.plot_mqq(
        sig_level=args.sig_level,
        save=args.output,
        save_args={"dpi": args.dpi, "facecolor": "white"}
    )
    
    print(f"Plots saved to: {args.output}")


if __name__ == '__main__':
    main()
