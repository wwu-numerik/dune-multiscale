echo "execute: make $* 2>&1 | repl.sh"
echo ""
make $* 2>&1 | ~/MyDune/DUNE/dune-multiscale/src/repl.sh
