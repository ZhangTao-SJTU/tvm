# Pick a path that won't move
LEGACY_DIR="$HOME/.venvs/legacy_tool"
python3.10.14 -m venv "$LEGACY_DIR"         # X.Y = the legacy Python that works
"$LEGACY_DIR/bin/pip" install "pyvoro"

