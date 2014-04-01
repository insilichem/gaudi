#!/usr/bin/sh
echo "This will run a script in Chimera in debug mode (no GUI)"
chimera="/home/jrodriguez/.local/UCSF-Chimera-1.8.1"

echo "Running script and args: $1"
echo "in session $2"

cd "$chimera"
"${chimera}/bin/chimera" --nogui --script "$1" "$2"