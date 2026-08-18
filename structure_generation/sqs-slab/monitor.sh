#!/bin/bash

# ATAT MCSQS Run with Integrated Progress Monitoring
# Auto-generated script with embedded file creation
# Generated configuration: slab_Ni_5x6x12.vasp, 12-site primitive cell x 5x6x1 supercell, 360 atoms

# --- Configuration ---
LOG_FILE="mcsqs1.log"
PROGRESS_FILE="mcsqs_parallel_progress.csv"
DEFAULT_MCSQS_ARGS="-rc"
TIME_LIMIT_MINUTES=0
TIME_LIMIT_SECONDS=$((TIME_LIMIT_MINUTES * 60))

# --- Auto-generate ATAT Input Files ---
create_input_files() {
   echo "Creating ATAT input files..."

   cat > rndstr.in << 'EOF'
0.047541 0.047541 1.000000 90.00 90.00 60.00
1 0 0
0 1 0
0 0 1
0.333333 0.333333 0.286506 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.666667 0.666667 0.325323 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.000000 0.000000 0.364140 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.333333 0.333333 0.402957 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.666667 0.666667 0.441774 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.000000 0.000000 0.480591 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.333333 0.333333 0.519409 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.666667 0.666667 0.558226 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.000000 0.000000 0.597043 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.333333 0.333333 0.635860 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.666667 0.666667 0.674677 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
0.000000 0.000000 0.713494 Al=0.033333,B=0.033333,C=0.033333,Cr=0.066667,Fe=0.033333,Mo=0.100000,Ni=0.700000
EOF

   cat > sqscell.out << 'EOF'
1

5 0 0
0 6 0
0 0 1
EOF

   echo "✅ Input files created: rndstr.in, sqscell.out"
}

# --- Monitoring Functions ---

extract_latest_objective() {
   grep "Objective_function=" "$1" | tail -1 | sed 's/.*= *//' 2>/dev/null || echo ""
}

extract_latest_step() {
   grep -c "Objective_function=" "$1" 2>/dev/null || echo "0"
}

extract_latest_correlation() {
   grep "Correlations_mismatch=" "$1" | tail -1 | sed 's/.*= *//' | awk '{print $1}' 2>/dev/null || echo ""
}

count_correlations() {
   grep "Correlations_mismatch=" "$1" | tail -1 | awk -F'\t' '{print NF-1}' 2>/dev/null || echo "0"
}

is_mcsqs_running() {
   pgrep -f "mcsqs" > /dev/null
   return $?
}

convert_bestsqs_to_poscar() {
    local bestsqs_file="$1"
    local poscar_file="$2"
    
    if [ ! -f "$bestsqs_file" ]; then
        echo "⚠️  Warning: $bestsqs_file not found"
        return 1
    fi
    
    echo "🔄 Converting $bestsqs_file to $poscar_file..."
    
    python3 - "$bestsqs_file" "$poscar_file" << 'PYEOF'
import sys
import numpy as np

def parse_bestsqs(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    A = np.array([[float(x) for x in lines[i].split()] for i in range(3)])
    B = np.array([[float(x) for x in lines[i].split()] for i in range(3, 6)])
    
    A_scaled = A * 52.355002
    final_lattice = np.dot(B, A_scaled)
    
    atoms = []
    for i in range(6, len(lines)):
        line = lines[i].strip()
        if line:
            parts = line.split()
            if len(parts) >= 4:
                x, y, z, element = float(parts[0]), float(parts[1]), float(parts[2]), parts[3]
                if element.lower() in ['vac', "'vac", 'vacancy', 'x']:
                    continue
                cart_pos = np.dot([x, y, z], A_scaled)
                atoms.append((element, cart_pos))
    
    return final_lattice, atoms

def write_poscar(lattice, atoms, filename, comment):
    from collections import defaultdict
    
    element_groups = defaultdict(list)
    for element, pos in atoms:
        element_groups[element].append(pos)
    
    atomic_weights = {
        'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00,
        'F': 19.00, 'Ne': 20.18, 'Na': 22.99, 'Mg': 24.31, 'Al': 26.98, 'Si': 28.09, 'P': 30.97, 'S': 32.07,
        'Cl': 35.45, 'Ar': 39.95, 'K': 39.10, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.00,
        'Mn': 54.94, 'Fe': 55.85, 'Co': 58.93, 'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38, 'Ga': 69.72, 'Ge': 72.63,
        'As': 74.92, 'Se': 78.96, 'Br': 79.90, 'Kr': 83.80, 'Rb': 85.47, 'Sr': 87.62, 'Y': 88.91, 'Zr': 91.22,
        'Nb': 92.91, 'Mo': 95.96, 'Tc': 98.00, 'Ru': 101.1, 'Rh': 102.9, 'Pd': 106.4, 'Ag': 107.9, 'Cd': 112.4,
        'In': 114.8, 'Sn': 118.7, 'Sb': 121.8, 'Te': 127.6, 'I': 126.9, 'Xe': 131.3, 'Cs': 132.9, 'Ba': 137.3,
        'La': 138.9, 'Ce': 140.1, 'Pr': 140.9, 'Nd': 144.2, 'Pm': 145.0, 'Sm': 150.4, 'Eu': 152.0, 'Gd': 157.3,
        'Tb': 158.9, 'Dy': 162.5, 'Ho': 164.9, 'Er': 167.3, 'Tm': 168.9, 'Yb': 173.0, 'Lu': 175.0, 'Hf': 178.5,
        'Ta': 180.9, 'W': 183.8, 'Re': 186.2, 'Os': 190.2, 'Ir': 192.2, 'Pt': 195.1, 'Au': 197.0, 'Hg': 200.6,
        'Tl': 204.4, 'Pb': 207.2, 'Bi': 209.0, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0,
        'Ac': 227.0, 'Th': 232.0, 'Pa': 231.0, 'U': 238.0, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0, 'Cm': 247.0,
        'Bk': 247.0, 'Cf': 251.0, 'Es': 252.0, 'Fm': 257.0, 'Md': 258.0, 'No': 259.0, 'Lr': 262.0
    }
    
    elements = sorted(element_groups.keys(), key=lambda x: atomic_weights.get(x, 999.0))
    
    with open(filename, 'w') as f:
        f.write(f'{comment}\n')
        f.write('1.0\n')
        
        for vec in lattice:
            f.write(f'  {vec[0]:15.9f} {vec[1]:15.9f} {vec[2]:15.9f}\n')
        
        f.write(' '.join(elements) + '\n')
        f.write(' '.join(str(len(element_groups[el])) for el in elements) + '\n')
        
        f.write('Direct\n')
        inv_lattice = np.linalg.inv(lattice)
        for element in elements:
            for cart_pos in element_groups[element]:
                frac_pos = np.dot(cart_pos, inv_lattice)
                f.write(f'  {frac_pos[0]:15.9f} {frac_pos[1]:15.9f} {frac_pos[2]:15.9f}\n')

try:
    import sys
    bestsqs_file = sys.argv[1] if len(sys.argv) > 1 else "$bestsqs_file"
    poscar_file = sys.argv[2] if len(sys.argv) > 2 else "$poscar_file"
    
    comment = f"SQS from {bestsqs_file}"
    lattice, atoms = parse_bestsqs(bestsqs_file)
    write_poscar(lattice, atoms, poscar_file, comment)
    print(f"✅ Successfully converted {bestsqs_file} to {poscar_file}")
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYEOF
    
    return $?
}

start_parallel_monitoring_process() {
   local output_file="$1"
   local minute=0

   echo "Monitor started for 5 parallel runs. Waiting for 5 seconds to allow mcsqs to initialize..."
   sleep 5

   header="Minute,Timestamp"
   for i in $(seq 1 5); do
       header="$header,Run${i}_Steps,Run${i}_Objective,Run${i}_Status"
   done
   header="$header,Best_Overall_Objective,Best_Run"
   echo "$header" > "$output_file"

   echo "----------------------------------------"
   echo "Monitoring 5 parallel MCSQS runs every minute"
   echo "Log files: mcsqs1.log, mcsqs2.log, ..., mcsqs5.log"
   echo "----------------------------------------"

   while true; do
       minute=$((minute + 1))
       local current_time=$(date +"%m/%d/%Y %H:%M")

       row_data="$minute,$current_time"
       best_objective=""
       best_run=""
       any_running=false

       for i in $(seq 1 5); do
           local log_file="mcsqs${i}.log"
           local objective="N/A"
           local step_count="0"
           local status="STOPPED"

           if pgrep -f "mcsqs.*-ip=${i}" > /dev/null; then
               status="RUNNING"
               any_running=true
           fi

           if [ -f "$log_file" ]; then
               objective=$(extract_latest_objective "$log_file")
               step_count=$(extract_latest_step "$log_file")
               objective=${objective:-"N/A"}
               step_count=${step_count:-"0"}
           fi

           row_data="$row_data,$step_count,$objective,$status"

           if [ "$objective" != "N/A" ] && [ -n "$objective" ]; then
               if [ -z "$best_objective" ] || awk "BEGIN {exit !($objective < $best_objective)}" 2>/dev/null; then
                   best_objective="$objective"
                   best_run="Run$i"
               fi
           fi
       done

       best_objective=${best_objective:-"N/A"}
       best_run=${best_run:-"N/A"}
       row_data="$row_data,$best_objective,$best_run"

       echo "$row_data" >> "$output_file"

       printf "Minute %3d | Active runs: " "$minute"
        for i in $(seq 1 5); do
            if pgrep -f "mcsqs.*-ip=${i}" > /dev/null; then
                printf "R%d " "$i"
            else
                printf "%s " "--"  
            fi
        done
        printf "| Best: %s (%s)\n" "$best_objective" "$best_run"

       if [ "$any_running" = false ]; then
           echo "All parallel runs stopped. Collecting final data..."
           break
       fi

       sleep 60
   done

   echo "----------------------------------------"
   echo "Parallel monitoring process finished."
}

# --- Main Script Logic ---

check_prerequisites() {
   echo "Checking prerequisites..."

   create_input_files

   echo "Generating clusters with corrdump..."
   echo "Command: corrdump -l=rndstr.in -ro -noe -nop -clus -2=0.069"
   corrdump -l=rndstr.in -ro -noe -nop -clus -2=0.069
   if [ $? -ne 0 ]; then
       echo "ERROR: corrdump command failed!"
       exit 1
   fi
   echo "✅ Clusters generated successfully."
   echo "✅ All prerequisites satisfied."
}
cleanup() {
   echo ""
   echo "=========================================="
   echo "🛑 Interrupt signal received or process completed"
   echo "=========================================="
   
   echo "🧹 Stopping MCSQS processes..."
   if [ -n "$MCSQS_PID" ]; then kill "$MCSQS_PID" 2>/dev/null; fi
   if [ -n "$MONITOR_PID" ]; then kill "$MONITOR_PID" 2>/dev/null; fi
   if [ -n "$TIMER_PID" ]; then kill "$TIMER_PID" 2>/dev/null; fi
   pkill -9 -f "mcsqs" 2>/dev/null || true
   sleep 2
   
   echo ""
   echo "=========================================="
   echo "📄 Converting bestsqs*.out files to POSCAR format..."
   echo "=========================================="
   
   found_files=0
   best_run=""
   best_objective=""
   
   if [ -f "$PROGRESS_FILE" ]; then
       last_line=$(tail -1 "$PROGRESS_FILE")
       best_objective=$(echo "$last_line" | cut -d',' -f$((3 + 3 * 5)))
       best_run=$(echo "$last_line" | cut -d',' -f$((4 + 3 * 5)) | sed 's/Run//')
   fi
   
   for outfile in bestsqs*.out; do
       if [ -f "$outfile" ]; then
           found_files=1
           basename="${outfile%.out}"
           poscar_file="${basename}_POSCAR"
           
           if convert_bestsqs_to_poscar "$outfile" "$poscar_file"; then
               echo "  ✅ $outfile → $poscar_file"
           else
               echo "  ❌ Failed to convert $outfile"
           fi
       fi
   done
   
   if [ $found_files -eq 0 ]; then
       echo "  ⚠️  No bestsqs*.out files found"
   else
       if [ -n "$best_run" ] && [ -n "$best_objective" ]; then
           if [ 5 -gt 1 ]; then
               best_poscar="bestsqs${best_run}_POSCAR"
           else
               best_poscar="bestsqs_POSCAR"
           fi
           
           if [ -f "$best_poscar" ]; then
               cp "$best_poscar" "POSCAR_best_overall"
               echo ""
               echo "🏆 Best structure (objective: $best_objective) saved as POSCAR_best_overall"
               if [ 5 -gt 1 ]; then
                   echo "    Source: Run $best_run (bestsqs${best_run}.out)"
               else
                   echo "    Source: bestsqs.out"
               fi
           else
               echo ""
               echo "⚠️  Could not find best POSCAR file: $best_poscar"
           fi
       else
           echo ""
           echo "⚠️  Could not determine best structure (no progress data found)"
       fi
       
       echo ""
       echo "=========================================="
       echo "✅ Conversion complete!"
       echo "=========================================="
   fi
   
   exit 0
}

trap cleanup SIGINT SIGTERM

# --- Execution ---
echo "================================================"
echo "    ATAT MCSQS with Integrated Monitoring"
echo "================================================"
echo "Configuration:"
echo "  - Structure: slab_Ni_5x6x12.vasp"
echo "  - Supercell: 5×6×1 of 12-site primitive cell (360 atoms)"
echo "  - Parallel runs: 5"
echo "  - Command: mcsqs -rc -ip=1 & mcsqs -rc -ip=2 & ... (parallel execution)"
echo "  - Time limit: None (manual stop)"
echo "  - Log file: $LOG_FILE"
echo "  - Progress file: $PROGRESS_FILE"
echo "================================================"

check_prerequisites

rm -f "$LOG_FILE" "$PROGRESS_FILE" mcsqs*.log
echo ""
echo "Starting ATAT MCSQS optimization and progress monitor..."

mcsqs -rc -ip=1 > mcsqs1.log 2>&1 &
mcsqs -rc -ip=2 > mcsqs2.log 2>&1 &
mcsqs -rc -ip=3 > mcsqs3.log 2>&1 &
mcsqs -rc -ip=4 > mcsqs4.log 2>&1 &
mcsqs -rc -ip=5 > mcsqs5.log 2>&1 &
MCSQS_PID=$!



start_parallel_monitoring_process "$PROGRESS_FILE" &
MONITOR_PID=$!

echo "✅ MCSQS started"
echo "✅ Monitor started (PID: $MONITOR_PID)"
echo ""
echo "Real-time progress logged to: $PROGRESS_FILE"
if [ $TIME_LIMIT_MINUTES -gt 0 ]; then
    echo "⏱️  Will auto-stop after $TIME_LIMIT_MINUTES minutes"
    echo "Press Ctrl+C to stop earlier and auto-convert to POSCAR."
else
    echo "Press Ctrl+C to stop optimization and auto-convert to POSCAR."
fi
echo "================================================"

wait
MCSQS_EXIT_CODE=$?

echo ""
if [ $MCSQS_EXIT_CODE -eq 124 ]; then
    echo "⏱️  Time limit reached ($TIME_LIMIT_MINUTES minutes). MCSQS stopped automatically."
else
    echo "MCSQS process finished with exit code: $MCSQS_EXIT_CODE."
fi

echo "Allowing monitor to capture final data..."
sleep 5

kill $MONITOR_PID 2>/dev/null
wait $MONITOR_PID 2>/dev/null

echo ""
echo "=========================================="
echo "🔄 Converting bestsqs files to POSCAR..."
echo "=========================================="

found_files=0
best_run=""
best_objective=""

if [ -f "$PROGRESS_FILE" ]; then
    last_line=$(tail -1 "$PROGRESS_FILE")
    best_objective=$(echo "$last_line" | cut -d',' -f$((3 + 3 * 5)))
    best_run=$(echo "$last_line" | cut -d',' -f$((4 + 3 * 5)) | sed 's/Run//')
fi

for outfile in bestsqs*.out; do
    if [ -f "$outfile" ]; then
        found_files=1
        basename="${outfile%.out}"
        poscar_file="${basename}_POSCAR"
        
        if convert_bestsqs_to_poscar "$outfile" "$poscar_file"; then
            echo "  ✅ $outfile → $poscar_file"
        else
            echo "  ❌ Failed to convert $outfile"
        fi
    fi
done

if [ $found_files -eq 0 ]; then
    echo "  ⚠️  No bestsqs*.out files found"
else
    if [ -n "$best_run" ] && [ -n "$best_objective" ]; then
        if [ 5 -gt 1 ]; then
            best_poscar="bestsqs${best_run}_POSCAR"
        else
            best_poscar="bestsqs_POSCAR"
        fi
        
        if [ -f "$best_poscar" ]; then
            cp "$best_poscar" "POSCAR_best_overall"
            echo ""
            echo "🏆 Best structure (objective: $best_objective) saved as POSCAR_best_overall"
            if [ 5 -gt 1 ]; then
                echo "    Source: Run $best_run (bestsqs${best_run}.out)"
            else
                echo "    Source: bestsqs.out"
            fi
        else
            echo ""
            echo "⚠️  Could not find best POSCAR file: $best_poscar"
        fi
    else
        echo ""
        echo "⚠️  Could not determine best structure (no progress data found)"
    fi
fi

echo ""
echo "================================================"
echo "              Optimization Complete"
echo "================================================"


if [ -f "$PROGRESS_FILE" ]; then
   echo "Progress Summary:"
   echo "  - Total monitoring time:   ~$(tail -1 "$PROGRESS_FILE" | cut -d',' -f1) minutes"
   echo "  - Best overall objective:  $(tail -1 "$PROGRESS_FILE" | cut -d',' -f$((3 + 3 * 5)))"
fi

echo "================================================"
