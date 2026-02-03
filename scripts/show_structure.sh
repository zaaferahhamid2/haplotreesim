#!/bin/bash
# Quick script to display the repository structure

echo "HaploTreeSim Repository Structure"
echo "=================================="
echo ""
echo "Total Files: 16 core files"
echo ""

tree -L 3 -I '__pycache__|*.pyc|*.egg-info' /mnt/user-data/outputs/haplotreesim/ || \
find /mnt/user-data/outputs/haplotreesim/ -type f -name "*.py" -o -name "*.md" -o -name "*.txt" -o -name "*.yml" | \
  grep -v __pycache__ | \
  grep -v .egg-info | \
  sort | \
  sed 's|/mnt/user-data/outputs/haplotreesim/||g' | \
  awk '{print "  " $0}'
