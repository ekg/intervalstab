# intervalstab

FastStabbing based on the STL

## building

`cmake -H. -Bbuild && cmake --build build -- -j 4`

## testing

`bin/intervalstab -T x -s 20000 -M 200 -m 10 -D 0 -S 233282`

## acknowledgements

This is a fork of code produced for this paper on optimal structures to solve the interval stabbing problem.

    Schmidt, Jens M. "Interval stabbing problems in small integer ranges."
    International Symposium on Algorithms and Computation. Springer, Berlin, Heidelberg, 2009.
