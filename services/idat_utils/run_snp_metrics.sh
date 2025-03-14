#!/bin/bash
# Simple script to run test_snp_metrics.py with hardcoded paths

# Start timing
start_time=$(date +%s)

# Check if we should run in background mode
if [ "$1" == "--background" ]; then
  # Run with nohup to keep it running if terminal disconnects
  nohup bash -c "
    echo \"Starting process at \$(date)\" > output.log
    python test_snp_metrics.py >> output.log 2>&1
    exit_code=\$?
    end_time=\$(date +%s)
    runtime=\$((end_time - $start_time))
    
    echo \"\" >> output.log
    echo \"Process completed at \$(date)\" >> output.log
    echo \"Exit code: \$exit_code\" >> output.log
    echo \"Total runtime: \$runtime seconds (\$((runtime / 60)) minutes)\" >> output.log
  " &

  pid=$!
  echo "Process started with PID: $pid"
  echo "Monitor progress with: tail -f output.log"
else
  # Run in foreground
  echo "Starting process at $(date)" | tee output.log
  python test_snp_metrics.py 2>&1 | tee -a output.log
  exit_code=$?
  end_time=$(date +%s)
  runtime=$((end_time - start_time))
  
  echo "" | tee -a output.log
  echo "Process completed at $(date)" | tee -a output.log
  echo "Exit code: $exit_code" | tee -a output.log
  echo "Total runtime: $runtime seconds ($((runtime / 60)) minutes)" | tee -a output.log
  
  echo "==================================================="
  echo "Process complete! See output above for details."
  echo "==================================================="
fi

# Test processing a single sample with default settings
# python test_snp_metrics.py --mode single

# # Test parallel processing of multiple samples with verbose output
# python test_snp_metrics.py --mode parallel --verbose

# # Test with a custom configuration file and keep intermediate files
# python test_snp_metrics.py --mode single --config my_config.ini --no-cleanup



