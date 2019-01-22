git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug
make 
cd ../..
git clone https://github.com/lh3/wgsim.git
cd wgsim
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
