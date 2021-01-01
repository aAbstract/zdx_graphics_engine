rm ./bin/engine
echo "[BASH-INFO]: Removed old binary"
echo "[BASH-INFO]: Compiling..."
./compile.sh ./bin/engine main.cpp
echo "[BASH-INFO]: Compiled new binary"
echo "[BASH-INFO]: Testing..."
./run.sh ./bin/engine