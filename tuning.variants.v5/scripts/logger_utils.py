#!/usr/bin/env python3
"""
Logging utilities for pipeline scripts
Provides dual output to both stdout and log file
"""

import sys
import os

class DualLogger:
    """Logger that writes to both stdout and a file"""
    
    def __init__(self, log_file):
        self.log_file = log_file
        self.terminal = sys.stdout
        self.log = open(log_file, 'w')
        
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.flush()
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()
    
    def close(self):
        self.log.close()

def setup_logging(output_file):
    """
    Setup logging to write to both console and file
    Returns the log file path
    """
    # Generate log filename from output file
    if output_file.endswith('.tsv'):
        log_file = output_file.replace('.tsv', '.log')
    elif output_file.endswith('.txt'):
        log_file = output_file.replace('.txt', '.log')
    elif output_file.endswith('.png'):
        # For plot scripts, use prefix
        log_file = output_file.rsplit('.', 1)[0] + '.log'
    else:
        log_file = output_file + '.log'
    
    # Ensure directory exists
    log_dir = os.path.dirname(log_file)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)
    
    # Redirect stdout to dual logger
    sys.stdout = DualLogger(log_file)
    
    return log_file
