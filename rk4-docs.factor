USING: help.markup help.syntax rk4.examples ;
IN: rk4

HELP: <rk4>
{ $description "Simple runge-kutta implementation for generating approximated solutions for a set of first order differential equations" }
{ $examples
    "A lorenz attractor is a popular system to model with this: "
    { $code "USING: rk4 rk4.examples ;" "lorenz." }
    "note that the produced chart is a 2 dimensional representation of a 3 dimensional solution. "
    "Similarly, the rabinovich-fabrikant system (stable alpha-gamma limit cycle): "
    { $code "USING: rk4 rk4.examples ;" "rabinovich-fabrikant." }
} ;

