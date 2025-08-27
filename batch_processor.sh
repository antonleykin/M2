#!/bin/bash

# Script to find an available BATCH directory and launch gfan command
# Usage: ./batch_processor.sh [nBatch]

# Default number of batches if not provided
nBatch=${1:-10}

echo "Searching through BATCH-0 to BATCH-$((nBatch-1)) for available directory..."

# Search through BATCH directories
for i in $(seq 0 $((nBatch-1))); do
    batch_dir="BATCH-$i"
    
    # Check if directory exists
    if [ ! -d "$batch_dir" ]; then
        echo "Directory $batch_dir does not exist, skipping..."
        continue
    fi
    
    # Check if "done" or "processing" files exist
    if [ -f "$batch_dir/done" ] || [ -f "$batch_dir/processing" ]; then
        echo "Directory $batch_dir is busy (has 'done' or 'processing' file), skipping..."
        continue
    fi
    
    # Found an available directory
    echo "Found available directory: $batch_dir"
    echo "Creating 'processing' file to mark as busy..."
    touch "$batch_dir/processing"
    
    # Change to the batch directory
    cd "$batch_dir"
    
    # Check if input file exists
    if [ ! -f "_tropicalprevariety.input" ]; then
        echo "Error: _tropicalprevariety.input file not found in $batch_dir"
        rm -f "processing"
        cd ..
        continue
    fi
    
    echo "Launching gfan command in $batch_dir..."
    
    # Launch the gfan command
    gfan _tropicalprevariety -j32 --log1 --usevaluation --bits64 --halfopenrestrictions --matrixoutput < _tropicalprevariety.input > _tropicalprevariety.output 2> gfan.log
    
    # Check if command completed successfully
    if [ $? -eq 0 ]; then
        echo "gfan command completed successfully in $batch_dir"
        # Remove processing file and create done file
        rm -f "processing"
        touch "done"
        echo "Marked $batch_dir as completed"
    else
        echo "gfan command failed in $batch_dir"
        # Remove processing file on failure
        rm -f "processing"
    fi
    
    # Return to parent directory
    cd ..
    
    # Exit after processing one batch
    exit 0
done

echo "No available BATCH directories found"
exit 1