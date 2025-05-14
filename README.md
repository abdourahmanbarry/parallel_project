# parallel_project guide

To run the parallel version of this project, write: <br>
`g++ -fopenmp parallel_Fem.cpp -o pfem ` <br>
`./pfem`

To run the serial version of this project, write: <br>
`g++ -fopenmp serial_Fem.cpp -o sfem` <br>
`./sfem`

The pde problem that is being solved is of the form `-au'' = f`. <br>
you can go to main method in either file to vary `a` and `f`, which are simple lambda functions.
For this problem, `a` must always be positive.<br>

You can then chose which point you want to evaluate by assigning that value to the `s` variable.

