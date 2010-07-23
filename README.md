Description
===========
`vecplot` is a simple vector field plot generator.
It generates nice-looking PostScript output for quick field visualization.

Usage
=====
    vecplot [-a] [-k kmeans-attempts] [-n num-arrows] data-file > output.ps

Flags
=====
* `-a`: Make one arrow for each data point specified in the file (Use (a)ll the data).
* `-k`: Specify how many k-means clustering attempts to perform. Defaults to 1.
* `-n`: Specify the number of arrows to use. Defaults to 100. If there are fewer than this many data points in the input file, then one arrow is created for each point.

Input format
============
The input format should have one field sample per line.
Each field sample has the format

    x-coord  y-coord  field-x-component  field-y-component

In other words, this is just like Gnuplot's input format for `vec`, or Matlab's format for `quiver`.
The principal difference is that the points need not be on a regular grid.

Compiling
=========
    g++ -O2 *.cpp kmpp/*.cpp -o vecplot

License
=======
The code uses David Arthur's k-means clustering code, which has no license specified.
For everything else, say ... BSD.

Todo
====
* Detect when data is on a regular grid and perform a simple bi-cubic (or similar) interpolation of the data instead of using k-means clustering for arrow placement.
* Draw tick marks, nicely.
* Supply a scale bar for the vector field arrow sizes.
