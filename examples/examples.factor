USING: kernel io accessors arrays rk4 sequences math math.matrices 
       ui.gadgets ui.gadgets.charts ui.gadgets.panes ui.gadgets.charts.lines ui.theme ;
IN: rk4.examples

: lorenz-dx/dt ( tx..n -- dx )
    1 tail first2 ! discard t and z
    swap - 10 * ;

: lorenz-dy/dt ( tx..n -- dy )
    1 tail first3 ! discard t
    28 swap - [ swap ] dip * swap - ;

: lorenz-dz/dt ( tx..n -- dz )
    1 tail first3 ! discard t
    [ * ] dip 8/3 * - ;

: <lorenz> ( -- dx..n/dt delta tx..n t-limit )
    { [ lorenz-dx/dt ] [ lorenz-dy/dt ] [ lorenz-dz/dt ] } 0.01 { 0 0 1 21/20 } 150 ;


: lorenz. ( -- )
    chart new { { -20 20 } { -20 20 } } >>axes
    line new link-color >>color <lorenz> <rk4> { 0 3 } cols-except >>data add-gadget
    gadget. ;


:: rf-dx/dt ( tx..n gamma -- dx )
    tx..n 1 tail first3 :> z :> y :> x
    y z 1 - x sq + * gamma x * + ;

:: rf-dy/dt ( tx..n gamma -- dy )
    tx..n 1 tail first3 :> z :> y :> x
    x 3 z * 1 + x sq - * gamma y * + ;

:: rf-dz/dt ( tx..n alpha -- dz )
    tx..n 1 tail first3 :> z :> y :> x
    2 neg z * alpha x y * + * ;

: <rabinovich-fabrikant> ( gamma alpha -- dx..n/dt delta tx..n t-limit )
    [ [ '[ _ rf-dx/dt ] ] keep '[ _ rf-dy/dt ] ] dip
    '[ _ rf-dz/dt ]
    3array
    0.01 { 0 -1 0 0.5 } 150 ;


: rabinovich-fabrikant. ( -- )
    chart new { { -2 2 } { -2 2 } } >>axes
    line new link-color >>color 0.1 0.14 <rabinovich-fabrikant> <rk4> { 0 3 } cols-except >>data add-gadget
    gadget. ;
