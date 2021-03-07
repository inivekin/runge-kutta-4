! Copyright (C) 2020 Your name.
! See http://factorcode.org/license.txt for BSD license.
USING: kernel math sequences arrays fry vectors math.vectors 
        classes.struct raylib.ffi grouping assocs accessors 
        namespaces see formatting combinators alien.enums
        locals math.functions 
        ui.gadgets ui.gadgets.panes ui.gadgets.charts ui.gadgets.charts.lines 
        ui.theme
        ;
IN: sigvis

: countdown ( int -- array_int )
    [ dup 1 > 
        [ 1 - ] [ drop f ] if 
    ] follow ;

! apply a set of three elements on a set of 2
: second_order ( pairs i idp -- final acc )
    '[ swap dup last swap reverse 
         { 1 } append _ v* 0 [ + ] reduce 2array swap drop 
     ] accumulate ;

: setup-camera ( -- camera ) 
    S{ Camera2D { offset S{ Vector2 f 300.0 0.0 } }
                { target S{ Vector2 f 0 0 } }
                { rotation 0.0 }
                { zoom 10.0 }
     } ;

: setup-3d-camera ( -- cam )
    S{ Camera3D 
        { position S{ Vector3 f 2 2 0 } }
        { target S{ Vector3 f 0 0 0 } }
        { up S{ Vector3 f 0 1 0 } }
        { fovy 60.0 } { type 0 } 
     } 
     CAMERA_PERSPECTIVE enum>number >>type
     ;

: 2nd-order>coords ( pairs -- enumerated-pairs ) 
    [ first ] map zip-index
    [ first2 Vector2 <struct> rot >>x swap >>y ] map
    2 clump ;
: prep-window ( pairs -- quot )
    '[ [ window-should-close ] 
       [ begin-drawing BLACK clear-background dup
            begin-mode-2d _ 
                [ dup last swap first LIME draw-line-v ] each
            end-mode-2d end-drawing
       ] 
     ] ; inline

: draw-window ( seq -- quot: ( camera -- ) ) 
    '[ 
         [ window-should-close ] 
         [ 
             begin-drawing 
                 BLACK clear-background dup 
                 begin-mode-2d _ 
                     [ dup last swap first LIME draw-line-v ]
                     each 
                 end-mode-2d
             end-drawing 
         ]
         until drop close-window 
     ] ; inline

: default-input ( -- final seq )
    50 { -1 1000 } <repetition>
    50 { 1000.2 -100.3 } <repetition> append 
    { 10.01 -20.2 } { 0.01 0.9 -1.9 } 
    second_order 2nd-order>coords ;

: default-init ( -- final )
    default-input 640 480 "test" init-window setup-camera 
    60 set-target-fps
    swap draw-window call ;

: dynamic-input ( pid -- quot: ( -- final seq ) )
    '[ 50 { 0 0 } <repetition> 
       50 { 0 0 } <repetition> append
       { 10.01 -20.2 } _ second_order 2nd-order>coords
     ] ; inline
: draw-dynamic-line ( final seq -- )
    [ dup last swap first LIME draw-line-v ] each drop ;
: draw-dynamic-window ( cam -- )
    { 0 0 0 } swap [ window-should-close ] 
    [ 
        begin-drawing 
            BLACK clear-background dup
            begin-mode-2d 
                swap
                75 is-key-down [ [ 0.01 + ] map ] when
                74 is-key-down [ [ 0.01 - ] map ] when
                264 is-key-down [ 0 over nth 0.01 - swap 0 swap remove-nth 0 swap insert-nth ] when
                265 is-key-down [ 0 over nth 0.01 + swap 0 swap remove-nth 0 swap insert-nth ] when
                262 is-key-down [ 1 over nth 0.01 + swap 1 swap remove-nth 1 swap insert-nth ] when
                263 is-key-down [ 1 over nth 0.01 - swap 1 swap remove-nth 1 swap insert-nth ] when
                72 is-key-down [ 2 over nth 0.01 + swap 2 swap remove-nth 2 swap insert-nth ] when
                76 is-key-down [ 2 over nth 0.01 - swap 2 swap remove-nth 2 swap insert-nth ] when
                dup "%f %f %f" vsprintf 0 0 4 RED draw-text
                dup
                dynamic-input call
                draw-dynamic-line
                swap
            end-mode-2d
        end-drawing 
    ]
    until drop close-window drop ;
: dynamic-init ( -- )
    640 480 "test" init-window setup-camera
    60 set-target-fps
    draw-dynamic-window ;

: constant:s ( -- s )
    10 ;

: constant:b ( -- b )
    8/3 ;

: constant:r ( -- r )
    28 ;

: constant:ts ( -- ts )
    ! 2/15 ;
    1/100 ;

: dx/dt ( txyz  -- dx )
    1 tail first2 ! discard t and z
    swap - constant:s * ;

: dy/dt ( txyz -- dy )
    1 tail first3 ! discard t
    constant:r swap - [ swap ] dip * swap - ;

: discard-time ( txyz -- x y z )
    1 tail first3 ! discard t
    ;

: dz/dt ( txyz -- dz )
    1 tail first3 ! discard t
    [ * ] dip constant:b * - ;

: (runge-kutta) ( diff-funcs txyz timestep -- rk1 )
    [ swap cleave 3array ] dip v*n ; inline

: runge-kutta-4-transform ( txyz rk3 timestep -- rk4 )
    ! prefix v*
    prefix v+ ;
    
: runge-kutta-2-transform ( txyz rk1 timestep -- timestep txyz' )
    ! dup [ runge-kutta-4-transform 1/2 v*n ] dip ;
    dup [ prefix 1/2 v*n v+ ] dip ;

: runge-kutta-3-transform ( txyz rk2 timestep -- timestep txyz' )
    runge-kutta-2-transform ;

: runge-kutta ( diff-functions txyz timestep -- solution )
    ! store original time - timestep at back
    [ dup first ] dip [ - -rot ] keep
    ! first pass
    3dup (runge-kutta)
    ! second pass
    [ 3dup ] dip
    [ swap runge-kutta-2-transform (runge-kutta) ] keep
    ! third pass
    [ 3dup ] 2dip
    [  [ drop ] dip swap runge-kutta-3-transform (runge-kutta) ] 2keep
    ! fourth pass
    [ [ 2drop ] dip swap dup [ runge-kutta-4-transform ] dip (runge-kutta) ] 3keep
    ! 3 [ v+ ] times 1/6 swap v*n ! same as below but doesn't follow stack effects
    4array { 0 0 0 } [ v+ ] reduce 1/2 v*n ! FIXME should be 1/6 not 1/2
    ! prefix with original time - timestep
    swap prefix
    ; inline

! : runge-kutta ( diff-functions txyz timestep -- solution )
!     ! store original time - timestep at back
!     [ dup first ] dip [ - -rot ] keep
!     [ (runge-kutta) ] 3keep
!     [ [ pick ] dip runge-kutta-2-transform (runge-kutta) ] 3keep
!     [ [ pick ] dip runge-kutta-3-transform (runge-kutta) ] 3keep
!     [ pick ] dip runge-kutta-4-transform (runge-kutta)
!     4array { 0 0 0 } [ v+ ] reduce 1/2 v*n ! FIXME should be 1/6 not 1/2
!     ! prefix with original time - timestep
!     swap prefix
!     ; inline

! :: runge-kutta ( diff-functions txyz timestep -- solution )
    ! B
    ! diff-functions txyz timestep (runge-kutta) :> rk1
    ! diff-functions txyz rk1 timestamp runge-kutta-2-transform (runge-kutta) :> rk2
    ! diff-functions txyz rk2 timestamp runge-kutta-3-transform (runge-kutta) :> rk3
    ! diff-functions txyz rk3 timestamp runge-kutta-4-transform (runge-kutta) :> rk4
    ! { rk1 rk2 rk3 rk4 } [ v+ ] reduce 1/6 v*n
    ! txyz first timestep - prefix
    ! ; inline

: Vector3>txyz ( Vector3 tlimit -- seq )
    [ { [ x>> ] [ y>> ] [ z>> ] } cleave 3array ] dip
    prefix ;

: (lorenz) ( txyz -- seq )
    { [ dx/dt ] [ dy/dt ] [ dz/dt ] }
    swap
    constant:ts
    runge-kutta ;

: gamma:constant ( -- gamma )
    ! 0.87 ;
    0.1 ;

: alpha:constant ( -- alpha )
    ! 1.1 ;
    0.14 ;
    ! 0.98 ;
    ! 0.1 ;

:: test-dx/dt ( t x y z -- dx )
    ! y ; ! damped chaotic oscillator
    ! y ; ! duffing oscillator
    y z 1 - x dup * + * gamma:constant x * + ; ! rabinovich-fabrikant

:: test-dy/dt ( t x y z -- dy )
    ! y neg x sin - t sin + ;
    ! x sin x z mod t cos * y - + ; ! damped chaotic oscillator
    ! x x 2dup * * - y - t cos + ; ! duffing oscillator
    x 3 z * 1 + x dup * - * gamma:constant y * + ; ! rabinovich-fabrikant

:: test-dz/dt ( t x y z -- dz )
    ! 0 ; 
    ! 0.01 ;
    2 neg z * alpha:constant x y * + * ; ! rabinovich-fabrikant

: (differential) ( txyz -- seq )
    { [ first4 test-dx/dt ] [ first4 test-dy/dt ] [ first4 test-dz/dt ] }
    swap
    0.01 ! constant:ts
    runge-kutta ; inline

: lorenz>Vector3 ( seq -- seq' )
    [ 1 tail first3 Vector3 <struct-boa> ] map ;

: lorenz ( Vector3 tlimit -- seq )
    Vector3>txyz
    [
        dup
        (lorenz)
        [ 1 tail 0 prefix ] dip
        v+
        dup first 0 <=
        [ drop f ]
        [ ]
        if
    ] follow ;

: differential ( Vector3 tlimit -- seq )
    Vector3>txyz
    [
        dup
        (differential)
        [ 1 tail 0 prefix ] dip
        v+
        dup first 0 <=
        [ drop f ]
        [ ]
        if
    ] follow ; inline

: default-lorenz ( -- seq )
    S{ Vector3 f 0 1 21/20 } 150 lorenz ;

: scale-down ( seq -- seq' )
    [ { 1 1/10 1/10 1/10  } v* ] map ;

: <default-lorenz> ( xyz -- xyz seq )
    dup first3 Vector3 <struct-boa> 150 lorenz scale-down lorenz>Vector3 2 clump ;

: <differential> ( xyz -- xyz seq )
    dup first3 Vector3 <struct-boa> 150 differential 
    ! scale-down 
    lorenz>Vector3 2 clump ; inline 

: alt-lorenz ( -- seq )
    S{ Vector3 f 0 0.9 1.04 } 150 lorenz ;

: draw-lorenz ( cam -- )
    [ window-should-close ] 
    [ 
        begin-drawing 
            BLACK clear-background dup
            begin-mode-3d 
                dup update-camera
                default-lorenz lorenz>Vector3 2 clump
                [ first2 LIME draw-line-3d ] each
            end-mode-3d
        end-drawing 
    ]
    until drop close-window ;

: black-b16 ( -- Color )
    0x1d 0x1f 0x21 0xff Color <struct-boa>
    ;

: white-b16 ( -- Color )
    0xe0 0xe0 0xe0 0xff Color <struct-boa>
    ;

: draw-static-lorenz ( cam -- )
    default-lorenz scale-down lorenz>Vector3 2 clump
    swap
    [ window-should-close ] 
    [ 
        begin-drawing 
            black-b16 clear-background dup
            begin-mode-3d 
                dup update-camera
                10 1.0 draw-grid
                [ dup [ first2 white-b16 draw-line-3d ] each ] dip
            end-mode-3d
        end-drawing 
    ]
    until drop close-window 
    drop 
    ;

: draw-dynamic-lorenz ( cam -- )
    { 0 1 21/20 } <default-lorenz>
    rot
    [ window-should-close ] 
    [ 
        begin-drawing 
            black-b16 clear-background dup
            begin-mode-3d 
            dup update-camera
            -rot
            KEY_LEFT enum>number is-key-down [ drop { 0 0 1 } v- <default-lorenz> ] when
            KEY_RIGHT enum>number is-key-down [ drop { 0 0 1 } v+ <default-lorenz> ] when
            KEY_DOWN enum>number is-key-down [ drop { 0 1 0 } v- <default-lorenz> ] when
            KEY_UP enum>number is-key-down [ drop { 0 1 0 } v+ <default-lorenz> ] when
            10 1.0 draw-grid
            dup [ first2 white-b16 draw-line-3d ] each
            rot
            end-mode-3d
        end-drawing 
    ]
    until drop close-window 
    2drop
    ;

: draw-dynamic-diff ( cam -- )
    { 0 1 21/20 } <differential>
    rot
    [ window-should-close ] 
    [ 
        begin-drawing 
            black-b16 clear-background dup
            begin-mode-3d 
            dup update-camera
            -rot
            KEY_LEFT enum>number is-key-down [ drop { 0 0 0.1 } v- <differential> ] when
            KEY_RIGHT enum>number is-key-down [ drop { 0 0 0.1 } v+ <differential> ] when
            KEY_DOWN enum>number is-key-down [ drop { 0 0.1 0 } v- <differential> ] when
            KEY_UP enum>number is-key-down [ drop { 0 0.1 0 } v+ <differential> ] when
            10 1.0 draw-grid
            dup [ first2 white-b16 draw-line-3d ] each
            rot
            end-mode-3d
        end-drawing 
    ]
    until drop close-window 
    2drop
    ; inline

: lorenz-init ( -- )
    640 480 "lorenz test" init-window setup-3d-camera
    60 set-target-fps
    dup CAMERA_FIRST_PERSON enum>number set-camera-mode
    ! draw-lorenz ;
    ! draw-static-lorenz ;
    draw-dynamic-lorenz ; 

: diff-vis ( -- )
    640 480 "differential visualiser" init-window setup-3d-camera
    60 set-target-fps
    dup CAMERA_FIRST_PERSON enum>number set-camera-mode
    draw-dynamic-diff ; inline

: lorenz. ( -- )
    chart new { { -20 20 } { -20 20 } } >>axes line new link-color >>color default-lorenz [ 3 head 1 tail ] map >>data add-gadget gadget. ;

: duffing. ( -- )
    chart new { { -3 3 } { -3 3 } } >>axes line new link-color >>color S{ Vector3 f 0 1 21/20 } 60 differential [ 3 head 1 tail ] map >>data add-gadget gadget. ;

MAIN: lorenz-init
