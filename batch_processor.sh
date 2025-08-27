#!/bin/bash

# Script to find available BATCH directories and launch gfan commands concurrently
# Usage: ./batch_processor.sh [nBatch] [max_concurrent]

# Default number of batches and max concurrent processes if not provided
nBatch=${1:-10}
max_concurrent=${2:-3}

echo "Searching through BATCH-0 to BATCH-$((nBatch-1)) for available directories..."
echo "Processing up to $max_concurrent batches concurrently..."

# Function to process a single batch
process_batch() {
    local batch_dir="$1"
    local batch_num="$2"
    
    echo "[BATCH-$batch_num] Starting processing..."
    
    # Change to the batch directory
    cd "$batch_dir"
    
    # Check if input file exists
    if [ ! -f "_tropicalprevariety.input" ]; then
        echo "[BATCH-$batch_num] Error: _tropicalprevariety.input file not found"
        rm -f "processing"
        cd ..
        return 1
    fi
    
    echo "[BATCH-$batch_num] Launching gfan command..."
    
    # Launch the gfan command
    gfan _tropicalprevariety -j32 --log1 --usevaluation --bits64 --halfopenrestrictions --matrixoutput < _tropicalprevariety.input > _tropicalprevariety.output 2> gfan.log
    
    # Check if command completed successfully
    if [ $? -eq 0 ]; then
        echo "[BATCH-$batch_num] gfan command completed successfully"
        # Remove processing file and create done file
        rm -f "processing"
        touch "done"
        echo "[BATCH-$batch_num] Marked as completed"
    else
        echo "[BATCH-$batch_num] gfan command failed"
        # Remove processing file on failure
        rm -f "processing"
    fi
    
    # Return to parent directory
    cd ..
    
    echo "[BATCH-$batch_num] Finished processing"
}

# Array to store background process PIDs
declare -a pids=()
declare -a batch_nums=()

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
    
    # Start processing this batch in background
    process_batch "$batch_dir" "$i" &
    pid=$!
    pids+=($pid)
    batch_nums+=($i)
    
    echo "Started BATCH-$i with PID $pid"
    
    # Check if we've reached max concurrent processes
    if [ ${#pids[@]} -ge $max_concurrent ]; then
        echo "Reached maximum concurrent processes ($max_concurrent), waiting for one to complete..."
        
        # Wait for any background process to complete (compatible approach)
        while true; do
            for j in "${!pids[@]}"; do
                if ! kill -0 "${pids[j]}" 2>/dev/null; then
                    echo "BATCH-${batch_nums[j]} (PID ${pids[j]}) has completed"
                    unset pids[j]
                    unset batch_nums[j]
                    # Reindex arrays
                    pids=("${pids[@]}")
                    batch_nums=("${batch_nums[@]}")
                    break 2  # Break out of both loops
                fi
            done
            sleep 1  # Wait a second before checking again
        done
    fi
done

# Wait for all remaining background processes to complete
if [ ${#pids[@]} -gt 0 ]; then
    echo "Waiting for remaining ${#pids[@]} processes to complete..."
    for i in "${!pids[@]}"; do
        wait "${pids[i]}"
        echo "BATCH-${batch_nums[i]} (PID ${pids[i]}) has completed"
    done
fi

echo "All batch processing completed"