USING: kernel accessors fry combinators sequences sequences.generalizations
       arrays math math.vectors math.functions ;
IN: rk4

: runge-kutta-2-transform ( rk1 tx..n delta -- tx..n' delta )
    [ swap [ [ v*n ] keep prefix 1/2 v*n ] dip v+ 2 v*n ] keep ;

: runge-kutta-3-transform ( rk2 tx..n delta -- tx..n' delta )
    runge-kutta-2-transform ;

: runge-kutta-4-transform ( rk3 tx..n delta -- tx..n' delta )
    [ [ swap ] dip [ v*n ] keep prefix v+ ] keep ;

: (runge-kutta) ( delta tx..n dx..n/dt -- rk )
    [ swap ] dip dup length>> [ cleave ] dip narray swap v*n
    ; inline

: runge-kutta-differentials ( dx..n/dt -- seq )
    '[ _ (runge-kutta) ] ;

: runge-kutta-transforms ( tx..n delta dx..n/dt -- seq )
    -rot swap
    [ { [ ]
      [ runge-kutta-2-transform ]
      [ runge-kutta-3-transform ]
      [ runge-kutta-4-transform ] } ] dip
    '[ _ runge-kutta-differentials compose ] map
    [ 2curry ] 2with map ;

: runge-kutta-4 ( dx..n/dt delta tx..n -- seq )
    [
        ! set up the set of 4 equations
        ! NOTE differential functions and timestep are curried
        runge-kutta-transforms [ [ dup ] compose ] map

        ! reduce the set of equations into one call
        [ ] [ compose ] reduce

        ! this will always produce 4 outputs with a duplication of the last result
        ! NOTE the dup is kept in the transform so that function
        ! can be used in higher-order estimation
        call( -- x x x x x ) drop 4array 

        ! reduce the results to the estimated change for the timestep
        { 0 0 0 } [ v+ ] reduce
        1/6 v*n
    ] keep

    ! add the resulting vector with the previous state vector
    [ 1 tail v+ ] keep

    ! prefix the new time+delta-time
    first prefix
    ;

: time+delta ( dx..n/dt delta tx..n -- tx..n dx..n/dt delta )
    [ swap [ 0 ] 2dip '[ _ + ] change-nth ] 2keep ; inline

: time-limit-predicate ( t-limit -- quot: ( tx..n -- ? ) )
    '[ dup first _ <= ] ; inline

: <rk4> ( dxn..n/dt delta initial-state t-limit -- seq )
    time-limit-predicate [ [ time+delta runge-kutta-4 ] [ 3drop f ] if ] compose 2with follow ; inline
