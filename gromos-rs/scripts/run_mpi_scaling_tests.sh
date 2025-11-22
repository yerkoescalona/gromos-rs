#!/bin/bash
# Automated MPI scaling test script
#
# This script runs the MPI scaling benchmark across different numbers of
# processes and generates a summary report.
#
# Usage:
#   ./scripts/run_mpi_scaling_tests.sh [system_size] [steps]
#
# Examples:
#   ./scripts/run_mpi_scaling_tests.sh 1000 100    # Small system
#   ./scripts/run_mpi_scaling_tests.sh 5000 50     # Medium system
#   ./scripts/run_mpi_scaling_tests.sh 10000 20    # Large system

set -e

# Configuration
ATOMS=${1:-1000}
STEPS=${2:-100}
MAX_PROCS=${3:-8}
OUTPUT_DIR="scaling_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULT_FILE="${OUTPUT_DIR}/scaling_${ATOMS}atoms_${TIMESTAMP}.csv"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔═══════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║       GROMOS-RS MPI Scaling Test Suite                        ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check if MPI is available
if ! command -v mpirun &> /dev/null; then
    echo -e "${YELLOW}ERROR: mpirun not found. Please install OpenMPI or MPICH.${NC}"
    exit 1
fi

# Build the benchmark tool
echo -e "${GREEN}Building MPI scaling benchmark...${NC}"
cargo build --release --features use-mpi --bin mpi_scaling

if [ $? -ne 0 ]; then
    echo -e "${YELLOW}ERROR: Build failed. Make sure MPI libraries are installed.${NC}"
    exit 1
fi

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Write CSV header
echo "processes,atoms,pairs,steps,total_time_s,time_per_step_ms,pairs_per_sec" > "${RESULT_FILE}"

echo ""
echo -e "${GREEN}Configuration:${NC}"
echo "  Atoms:        ${ATOMS}"
echo "  Steps:        ${STEPS}"
echo "  Max procs:    ${MAX_PROCS}"
echo "  Output:       ${RESULT_FILE}"
echo ""

# Store baseline (1 process) time for speedup calculation
BASELINE_TIME=""

# Run tests with different numbers of processes
for NP in 1 2 4 8 16; do
    if [ ${NP} -gt ${MAX_PROCS} ]; then
        break
    fi

    echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
    echo -e "${BLUE}Testing with ${NP} process(es)${NC}"
    echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"

    # Run the benchmark
    RESULT=$(mpirun -np ${NP} ./target/release/mpi_scaling \
        --atoms ${ATOMS} \
        --steps ${STEPS} \
        --csv 2>/dev/null | grep -v "^$" | tail -1)

    # Store result
    echo "${RESULT}" >> "${RESULT_FILE}"

    # Parse and display
    TIME=$(echo ${RESULT} | cut -d',' -f5)
    TIME_PER_STEP=$(echo ${RESULT} | cut -d',' -f6)
    PAIRS_PER_SEC=$(echo ${RESULT} | cut -d',' -f7)

    echo "  Total time:       ${TIME} s"
    echo "  Time per step:    ${TIME_PER_STEP} ms"
    echo "  Pairs/sec:        ${PAIRS_PER_SEC}"

    # Calculate speedup vs baseline
    if [ ${NP} -eq 1 ]; then
        BASELINE_TIME=${TIME}
        echo "  Speedup:          1.00x (baseline)"
        echo "  Efficiency:       100.0%"
    else
        SPEEDUP=$(echo "scale=2; ${BASELINE_TIME} / ${TIME}" | bc)
        EFFICIENCY=$(echo "scale=1; ${SPEEDUP} * 100 / ${NP}" | bc)
        echo "  Speedup:          ${SPEEDUP}x"
        echo "  Efficiency:       ${EFFICIENCY}%"
    fi

    echo ""
done

# Generate summary report
echo -e "${GREEN}═══════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}Scaling Test Complete!${NC}"
echo -e "${GREEN}═══════════════════════════════════════════════════════════════${NC}"
echo ""
echo "Results saved to: ${RESULT_FILE}"
echo ""

# Generate a simple text report
REPORT_FILE="${OUTPUT_DIR}/scaling_${ATOMS}atoms_${TIMESTAMP}_report.txt"

cat > "${REPORT_FILE}" << EOF
GROMOS-RS MPI Scaling Report
Generated: ${TIMESTAMP}
System: ${ATOMS} atoms, ${STEPS} steps

Processes | Time (s) | Time/step (ms) | Speedup | Efficiency
----------|----------|----------------|---------|------------
EOF

# Parse results and calculate metrics
NP=1
while IFS=',' read -r procs atoms pairs steps time time_per_step pairs_per_sec; do
    # Skip header
    if [ "${procs}" == "processes" ]; then
        continue
    fi

    if [ ${procs} -eq 1 ]; then
        BASELINE_TIME=${time}
        SPEEDUP="1.00"
        EFFICIENCY="100.0"
    else
        SPEEDUP=$(echo "scale=2; ${BASELINE_TIME} / ${time}" | bc)
        EFFICIENCY=$(echo "scale=1; ${SPEEDUP} * 100 / ${procs}" | bc)
    fi

    printf "%9s | %8s | %14s | %7s | %10s%%\n" \
        "${procs}" "${time}" "${time_per_step}" "${SPEEDUP}x" "${EFFICIENCY}" >> "${REPORT_FILE}"

done < "${RESULT_FILE}"

cat "${REPORT_FILE}"
echo ""

# Check for good scaling
echo -e "${BLUE}Scaling Analysis:${NC}"
echo ""

# Read last line of results (highest process count)
LAST_LINE=$(tail -1 "${RESULT_FILE}")
LAST_PROCS=$(echo ${LAST_LINE} | cut -d',' -f1)
LAST_TIME=$(echo ${LAST_LINE} | cut -d',' -f5)

if [ -n "${BASELINE_TIME}" ] && [ -n "${LAST_TIME}" ]; then
    FINAL_SPEEDUP=$(echo "scale=2; ${BASELINE_TIME} / ${LAST_TIME}" | bc)
    FINAL_EFFICIENCY=$(echo "scale=1; ${FINAL_SPEEDUP} * 100 / ${LAST_PROCS}" | bc)

    echo "  Peak configuration: ${LAST_PROCS} processes"
    echo "  Peak speedup:       ${FINAL_SPEEDUP}x"
    echo "  Peak efficiency:    ${FINAL_EFFICIENCY}%"
    echo ""

    # Provide recommendations
    if (( $(echo "${FINAL_EFFICIENCY} > 70" | bc -l) )); then
        echo -e "  ${GREEN}✓ Excellent scaling! Consider using more processes.${NC}"
    elif (( $(echo "${FINAL_EFFICIENCY} > 50" | bc -l) )); then
        echo -e "  ${YELLOW}✓ Good scaling. ${LAST_PROCS} processes is near optimal.${NC}"
    else
        echo -e "  ${YELLOW}⚠ Poor scaling beyond ${LAST_PROCS} processes. System may be too small.${NC}"
        echo "    Consider:"
        echo "    - Increasing system size (--atoms)"
        echo "    - Reducing number of processes"
        echo "    - Checking network performance"
    fi
fi

echo ""
echo -e "${GREEN}All results saved to: ${OUTPUT_DIR}/${NC}"
echo ""
