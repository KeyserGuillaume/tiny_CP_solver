# tiny_CP_solver

This is a school exercise consisting in implementing a CP solver for small problems with binary variables.

Right now a possible improvement would be to add to the class ConstraintArc an access to the following information: "if I set u to its ith value, what are the values of v which become impossible ?"

In this way I would avoid iterating like an idiot over every value that v has, even though in both the n-queens problem and the graph coloring problem I already know the values which become impossible. The pb is that a value x which becomes impossible for v might already not be in v's domain so I will have to assume v's domain values are sorted and make the efficient search in a sorted list.
