namespace Chebfun
open System
open MathNet.Numerics.IntegralTransforms

// skip all these layers of abstraction for now:
// chebtech1 < chebtech < smoothfun < onefun
type chebtech2( c: double array , v: double , h: double , i: bool , e: double) = 
    // properties
    let mutable eps = e
    // members
    member this.coeffs = c
    member this.vscale = v
    member this.hscale = h
    member this.ishappy = i
    member this.epslevel with get() = eps and set(e) = eps <- e
    // dependent properties
    member this.length with get() = this.coeffs.Length
    member this.isempty with get() = Array.isEmpty(this.coeffs)

    static member chebpts( n ) =
        match n with 
            // Special case (no points)
            | 0 -> Array.empty<double>
            // Special case (single point)
            | 1 -> [|0.0|]
            // General case
            | _ ->
                // Chebyshev points
                let m = float(n) - 1.0
                Array.map( fun x -> Math.Sin(Math.PI * x / (2.0*m) ) ) [|-m..2.0..m|] // (Use of sine enforces symmetry.)

    member this.points with get() = chebtech2.chebpts(this.length)

    static member barywts( n ) =
        match n with
            // Special case (no points)
            | 0 -> Array.empty<double>
            // Special case (single point)
            | 1 -> [|1.0|]
            // General case
            | _ ->
                // Note v.[n-1] is positive.
                Array.append (Array.create (n-1) 1.0) [|0.5|] 
                |> Array.mapi( fun i x -> 
                    match i with
                    | 0 -> 0.5  * x
                    | i when [|(n-2)..(-2)..0|] |> Array.exists( fun x -> x = i ) -> -1.0
                    | _ -> x ) 

    static member clenshaw( x:double array , c:double array ) = 
        // Clenshaw scheme for scalar-valued functions.
        let mutable bk1 = Array.create x.Length 0.0;
        let mutable bk2 = bk1
        let mutable xb1 = bk1
        let mutable xb2 = bk1
        let n = c.Length - 1;
        for k in [(n)..(-2)..2] do
            xb1 <- Array.map2( fun x y -> x * y ) x bk1
            bk2 <- Array.map2( fun xb1 b2 -> c.[k] + 2.0*xb1 - b2 ) xb1 bk2
            xb2 <- Array.map2( fun x y -> x * y ) x bk2
            bk1 <- Array.map2( fun xb2 b1 -> c.[k-1] + 2.0*xb2 - b1 ) xb2 bk1
        if (n % 2 = 1) then
            xb1 <- Array.map2( fun x y -> x * y ) x bk1
            let tmp1 = Array.map2( fun xb1 b2 -> c.[1] + 2.0*xb1 - b2 ) xb1 bk2
            let tmp2 = bk1
            bk1 <- tmp1
            bk2 <- tmp2
        xb1 <- Array.map2( fun x y -> x * y ) x bk1
        Array.map2( fun xb1 b2 -> c.[0] + xb1 - b2 ) xb1 bk2

    static member feval( f:chebtech2 , x:double array ) =
        match f.isempty with
            | true -> Array.empty<double>
            | false -> chebtech2.clenshaw( x , f.coeffs )

    static member refine( op: double -> double , values: double array ) = 

        // Default refinement function for resampling scheme.
        let refineResampling( op: double -> double , values: double array ) = 

            let mutable giveUp = false
            let mutable n = 0

            if Array.isEmpty(values) then
                // Choose initial n based upon minSamples:
                n <- int(Math.Pow( 2.0, Math.Ceiling( Math.Log(17.0 - 1.0)/Math.Log(2.0) ) )) + 1
            else
                failwith "Resampling with values present is not implemented"

            // n is too large:
            if (n > 65537) then
                if Array.isEmpty(values) then
                    // don't give up if we haven't sampled at least once.
                    n <- 65537
                    giveUp <- false
                else
                    giveUp <- true
             else
                giveUp <- false

            match giveUp with
                | true -> ( values , giveUp )
                | false ->
                    // 2nd-kind Chebyshev grid:
                    let x = chebtech2.chebpts(n)

                    // Evaluate the operator:
                    let newValues = Array.map( op ) x
                    ( newValues , giveUp )

        // Default refinement function for single ('nested') sampling.
        let refineNested( op: double -> double , values: double array ) = 

            let mutable giveUp = false
            let mutable n = 0
            let mutable newValues = Array.empty<double>

            if Array.isEmpty(values) then
                refineResampling( op , values )
            else
                // Compute new n by doubling (we must do this when not resampling).
                n <- 2 * values.Length - 1

                // n is too large:
                if (n > 65537) then
                    giveUp <- true
                 else
                    giveUp <- false

                match giveUp with
                    | true -> ( values , giveUp )
                    | false ->
                        // 2nd-kind Chebyshev grid and take every 2nd element
                        let x = 
                            chebtech2.chebpts(n)
                            |> Array.mapi( fun i x -> (i , x) )
                            |> Array.filter( fun (i, x) -> i % 2 = 1 )
                            |> Array.map snd

                        // Evaluate the operator:
                        let newValues = Array.map( op ) x

                        // Merge stored values:
                        let allValues =
                            Array.append [|values.[0]|] (
                                Array.zip newValues values.[1..values.Length-1]
                                |> Array.collect( fun x -> [|fst(x) ; snd(x)|] ) )

                        ( allValues , giveUp )

        // (Don't yet) decide which refinement to use:
        refineNested( op , values )

    static member extrapolate( values: double array ) = 
        let maskNaN = Array.map( Double.IsNaN ) values
        let maskInf = Array.map( Double.IsInfinity ) values
        let mask = 
            Array.zip maskNaN maskInf 
            |> Array.map( fun x -> fst(x) || snd (x) )
        
        if Array.exists( fun x -> x = true ) mask then
            // Obtain Chebyshev points:
            let n = values.Length
            let x = chebtech2.chebpts(n)

            // The good and the bad
            let xbadTuple , xgoodTuple = 
                x
                |> Array.zip mask
                |> Array.partition( fun (m , x) -> m = true )
            let xbad = xbadTuple |> Array.map snd
            let xgood = xgoodTuple |> Array.map snd

            // if Array.isEmpty(xgood) then
            if not(Array.isEmpty(xbad)) then
                failwith "Too many NaNs/Infs to handle."
            // else
                // // Compute the modified barycentric weights:
                // let mutable w = chebtech2.barywts(n) // Standard weights
                // w <- w 
                // |> Array.zip mask 
                // |> Array.filter( fun (m , x) -> m = false ) 
                // |> Array.map snd // Barycentric weights corresponding to the good points.
                // ...

        ( values , maskNaN , maskInf )

    static member vals2coeffs( values: double array ) = 
        // Get the length ofthe input:
        let n = values.Length

        match n with
            | 0 | 1 -> values // Trivial case (constant)
            | _ -> 
                // Mirror the values (to fake a DCT using an FFT):
                let coeffs = 
                    (Array.append (Array.rev values.[1..values.Length-1]) (Array.rev (Array.rev values).[1..values.Length-1]))
                    |> Array.map( fun x -> Numerics.Complex( x , 0.0 ) )
                // Real-valued case (and truncate and scale the interior coefficients):
                Fourier.Inverse( coeffs , FourierOptions.Matlab )
                (Array.map( fun (x:Numerics.Complex) -> double(x.Real) )  coeffs)
                |> Array.mapi( fun i x -> ( i , x ) )
                |> Array.filter( fun ( i , x ) -> i < n )
                |> Array.map( fun ( i , x ) -> 
                    match i with
                    | 0 -> x
                    | i when i = n-1 -> x
                    | _ -> 2.0 * x )
                // TODO Imaginary-valued case:
                // TODO General case:

    static member coeffs2vals( coeffs: double array ) = 
        // Get the length ofthe input:
        let n = coeffs.Length

        match n with
            | 0 | 1 -> coeffs // Trivial case (constant or empty)
            | _ -> 
                // Scale the coefficients by 1/2 and mirror them (to fake a DCT using an FFT):
                let values = 
                    coeffs
                    |> Array.mapi( fun i x -> ( i , x ) )
                    |> Array.map( fun ( i , x ) -> 
                        match i with
                        | 0 -> x
                        | i when i = n - 2 -> x
                        | _ -> x/2.0 )
                    |> fun x -> Array.append x (Array.rev x.[1..x.Length-2])
                    |> Array.map( fun x -> Numerics.Complex( x , 0.0 ) )
                // Real-valued case (and flip and truncate):
                Fourier.Forward( values , FourierOptions.Matlab )
                (Array.map( fun (x:Numerics.Complex) -> double(x.Real) ) values)
                |> Array.mapi( fun i x -> ( i , x ) )
                |> Array.filter( fun ( i , x ) -> i < n )
                |> Array.rev
                |> Array.map snd
                // TODO Imaginary-valued case:
                // TODO General case:

    static member alias( coeffs:double array , m:int ) =
        let n = coeffs.Length

        match m with
            | m when m > n -> // Pad with zeros
                Array.append coeffs (Array.create (m-n) 0.0)
            // Alias coefficients: (see eq. (4.4) of Trefethen, Approximation Theory and
            // Approximation Practice, SIAM, 2013):
            | 1 -> // Reduce to a single point
                let e =Array.map( fun x -> Math.Pow(-1.0,x) ) [|0.0..Math.Ceiling(double(n)/2.0)|]
                let c = 
                    coeffs 
                    |> Array.mapi( fun i x -> ( i , x ) ) 
                    |> Array.filter( fun ( i , x ) -> i % 2 = 0)
                    |> Array.map snd
                Array.map2( fun x y -> x * y ) e c // and truncate:
                |> Array.mapi( fun i x -> ( i , x ) )
                |> Array.filter( fun ( i , x ) -> i < m )
                |> Array.map snd
            | m when m > n/2 -> // If m > n/2, only single coefficients are aliased, and we can vectorise.
                let j = [|m..(n-1)|]
                let k = Array.map( fun j -> Math.Abs( (j + m - 3 + 1) % (2*m - 2) - m + 2 ) ) j
                coeffs
                |> Array.mapi( fun i x -> 
                    match i with 
                    | i when k|> Array.exists( fun k -> k = i ) 
                        -> x + coeffs.[j.[ k |> Array.findIndex( fun k -> k = i ) ]] 
                    | _ -> x ) // and truncate:
                |> Array.mapi( fun i x -> ( i , x ) )
                |> Array.filter( fun ( i , x ) -> i < m )
                |> Array.map snd
            | _ -> // Otherwise we must do everything in a tight loop. (Which is slower!)
                let j = [|m..(n-1)|]
                let k = Array.map( fun j -> Math.Abs( (j + m - 3 + 1) % (2*m - 2) - m + 2 ) ) j
                coeffs
                |> Array.mapi( fun i x -> 
                    match i with 
                    | i when j|> Array.exists( fun j -> j = i ) 
                        -> x + coeffs.[k.[ j |> Array.findIndex( fun j -> j = i ) ]] 
                    | _ -> x ) // and truncate:
                |> Array.mapi( fun i x -> ( i , x ) )
                |> Array.filter( fun ( i , x ) -> i < m )
                |> Array.map snd   

    static member happinessCheck( f: chebtech2 , op: double -> double , values:double array) = 
        
        let happinessRequirements( values:double array , coeffs:double array , x:double array , vscale:double , hscale:double , epslevel:double ) =
            // Grab the size:
            let n = values.Length

            // We will not allow the estimated rounding errors to be cruder than this value:
            let minPrec = 1e-4 // Worst case precision!

            // Length of tail to test.
            let testLength = Math.Min( n , Math.Max( 5 , int(Math.Round( (double(n)-1.0)/8.0 )) ) )

            // Look at length of tail to loosen tolerance:
            let tailErr = Math.Min( 2.2e-16 * double(testLength) , minPrec )

            // Estimate the condition number of the input function by
            //    ||f(x+eps(x)) - f(x)||_inf / ||f||_inf ~~ (eps(hscale)/vscale)*f'.            
            let dy = Array.map2( fun x y -> y - x ) values.[0..n-2] values.[1..n-1]
            let dx = Array.map2( fun x y -> y - x ) x.[0..n-2] x.[1..n-1]
            let gradEst = // Finite difference approx.
                Array.map2( fun x y -> y/x ) dx dy
                |> Array.map( Math.Abs )
                |> Array.max
            let condEst = // Condition number estimate, using estimate of MATLAB's eps(hscale)
                Math.Min( 2.2e-16 * hscale / vscale * gradEst , minPrec )

            // Choose maximum between prescribed tolerance and estimated rounding error:
            let epslevelOut = Math.Max( Math.Max( epslevel , condEst ) , tailErr )                    

            ( testLength , epslevelOut )
        
        // What does happiness mean to you?
        let classicCheck ( f:chebtech2 , values: double array ) = 
            // Deal with special cases

            // Determine n (the length of the input)
            let n = f.length

            // Grab some preferences:
            let mutable epslevel = 2.2e-16

            // Deal with the trivial case
            match n with
            | n when n < 2 -> ( false , epslevel , n )
            | _ -> // Check the vertical scale
                if f.vscale = 0.0 then
                    // this is the zero function, so we must be happy!
                    ( true , epslevel , 1 )
                elif Double.IsInfinity(f.vscale) then
                    // Inf located. No cutoff.
                    ( false , epslevel , n )
                else
                    // NaNs are not allowed
                    if Array.exists( Double.IsNaN ) f.coeffs then failwith "Function returned NaN when evaluated"

                    // Compute some values if none were given
                    if Array.isEmpty( values ) then failwith "No values were given"
                        // let values = chebtech2.coeffs2vals( f.coeffs )

                    // Check for convergence and chop location

                    // Absolute value of coefficients, relative to vscale:
                    let ac = Array.map( fun (x:double) -> Math.Abs(x)/f.vscale ) f.coeffs
                    
                    // Happiness requirements:
                    let testLength , epslevel = happinessRequirements( values, f.coeffs, f.points, f.vscale, f.hscale, epslevel)
                    
                    let tail = 
                        ac
                        |> Array.mapi( fun i x -> ( i , x ) )
                        |> Array.filter( fun ( i , x ) -> i >= n - testLength )
                        |> Array.map snd
                    if (Array.max tail < epslevel) then // We have converged! Chop tail:
                        // Find last row of coeffs with entry above epslevel
                        let Tloc = 
                            try 
                                n - ( ac
                                |> Array.map( fun x -> x > epslevel )
                                |> Array.rev
                                |> Array.findIndex( fun x -> x) )
                            with :? System.Collections.Generic.KeyNotFoundException -> -1
                        match Tloc with
                            | -1 -> ( true , epslevel , 1 ) // Check for the zero function!
                            | _ -> 
                                // Compute the cumulative max of eps/4 and the tail entries:
                                let t = 0.25 * 2.2e-16
                                let tmp =
                                    ac
                                    |> Array.mapi( fun i x -> ( i , x ) )
                                    |> Array.filter( fun ( i , x ) -> i >= Tloc )
                                    |> Array.map snd
                                let cumMax = // Restrict to coefficients of interest to get cumulative maximum
                                    Array.scan( fun (acc:double) (x:double) -> Math.Max( acc , x ) ) t tmp
                                    |> fun x -> x.[1..x.Length-1] // remove first element

                                // Obtain an estimate for how much accuracy we'd gain 
                                // compared to reducing length ("bang for buck"):
                                let bang = Array.map( fun x -> Math.Log( 1e3*epslevel/x ) ) cumMax
                                let buck = [|(double(n)-1.0)..(-1.0)..double(Tloc)|]
                                let Tbpb = Array.map2( fun x y -> x/y ) bang buck

                                // Compute position at which to chop.
                                let bangForBuck =
                                    Tbpb
                                    |> Array.mapi( fun i x -> ( i , x ) )
                                    |> Array.filter( fun ( i , x ) -> 2 <= i && i <= n-1-Tloc )
                                    |> Array.map snd 
                                    |> Array.max
                                let Tchop = Tbpb |> Array.findIndex( fun x -> x = bangForBuck )

                                ( true , epslevel , n - (Tchop-1) - 2 ) 
                        
                    else // We're unhappy. :(
                        ( false , Array.average tail , 0 )

        let sampleTest( op: double -> double , values:double array , f:chebtech2 ) = 
            // Get the interpolation points:
            let n = f.length
            let x = chebtech2.chebpts(n);

            // Set a tolerance:
            let tol = Math.Max( f.epslevel, 1e3 * 2.2e-16 ) * double(n)

            // Choose a point to evaluate at:
            let xeval = 
                match n with
                | 1 -> 0.61 // Pseudo-random test value
                | _ ->  // Test a point where the (finite difference) gradient of values is largest:
                    let dy = Array.map2( fun x y -> y - x ) values.[0..n-2] values.[1..n-1]
                    let dx = Array.map2( fun x y -> y - x ) x.[0..n-2] x.[1..n-1]
                    let gradEst = // Finite difference approx.
                        Array.map2( fun x y -> y/x ) dx dy
                        |> Array.map( Math.Abs )
                    let index = int( gradEst |> Array.findIndex( fun x -> x = (Array.max gradEst) ) )
                    (x.[index+1] - 1.41*x.[index])/2.41
            let xevals = [| -1.0+1e-12 ; xeval ; 1.0-1e-12 |]

            // Evaluate the CHEBTECH:
            let vFun = chebtech2.feval(f, xevals)

            // Evaluate the op:
            let vOp = Array.map( op ) xevals

            // If the CHEBTECH evaluation differs from the op evaluation, SAMPLETEST failed:
            if Array.forall(fun (x:double) -> Math.Abs(x) < tol ) (Array.map2( fun x y -> ( y - x )/f.vscale ) vFun vOp) then
                true // :)
            else
                false // :(

        // TODO: allow preferences to choose strict, loose, and plateau
        let mutable ishappy , epslevel , cutoff = classicCheck( f , values )

        // Check also that sampleTest is happy:
        if ishappy then
            f.epslevel <- epslevel
            let ishappy = sampleTest( op , values , f )
            if not(ishappy) then
                cutoff <- values.Length

        ( ishappy , epslevel , cutoff )

    static member simplify( f:chebtech2 , tol:double ) = 
        if f.isempty then
            // Deal with empty case:
            ( f.coeffs , f.epslevel )
        elif not(f.ishappy) then
            // Do nothing to an unhappy CHEBTECH:
            ( f.coeffs , f.epslevel )
        else
            // Check for trailing coefficients smaller than 
            // the tolerance relative to F.VSCALE:
            let largeCoeffs = Array.map( fun (x:double) -> (Math.Abs(x) - tol * f.vscale) > 0.0 ) f.coeffs
            let firstNonZeroRow = Array.tryFindIndex( fun x -> x = true ) (Array.rev largeCoeffs)

            // If the whole thing is now zero, leave just one coefficient:
            let coeffs = 
                match firstNonZeroRow with 
                | None -> [|0.0|]
                | _ -> // or remove trailing zeros:
                    f.coeffs
                    |> Array.mapi( fun i x -> ( i , x ) )
                    |> Array.filter( fun ( i , x ) -> i < f.length - firstNonZeroRow.Value )
                    |> Array.map snd

            // Update epslevel:
            let epslevel = Math.Max( f.epslevel , tol )
            
            ( coeffs , epslevel )
                
    // empty constructor
    new() = chebtech2( Array.empty<double> , 0.0 , 1.0 , new bool() , new double() )

    // anon constructor
    new( op: double -> double ) =
        let mutable vscale = 0.0
        let mutable hscale = 1.0 // Different from Chebfun's:
        // "% TODO:  Why do we rescale the hscale like this?
        // data.hscale = data.hscale / diff(data.domain);"
        let mutable ishappy = false
        let mutable giveUp = false
        let mutable epslevel = 2.0e-16

        // adaptive construction
        let mutable values = Array.empty<double>
        let mutable coeffs = Array.empty<double>

        // Loop until ISHAPPY or GIVEUP:
        while (not(ishappy) && not(giveUp)) do
            // Call the appropriate refinement routine:
            let (newValues , newGiveUp) = chebtech2.refine( op , values )
            values <- newValues
            giveUp <- newGiveUp

            match giveUp with
                | true -> ()
                | false ->
                    // Update vertical scale: (Only include sampled finite values)
                    vscale <- Math.Max( vscale , 
                        values
                        |> Array.filter( fun x -> not(Double.IsInfinity(x)) )
                        |> Array.map( Math.Abs )
                        |> Array.max )

                    // Extrapolate out NaNs: (not done)
                    let newValues , maskNaN , maskInf = chebtech2.extrapolate( values )
                    values <- newValues

                    // Compute the Chebyshev coefficients:
                    coeffs <- chebtech2.vals2coeffs( values )
                        
                    // Check for happiness:
                    let f = chebtech2( coeffs , vscale , hscale , ishappy , epslevel )
                    let newIshappy, newEpslevel, cutoff = chebtech2.happinessCheck(f, op, values) 
                    ishappy <- newIshappy
                    epslevel <- newEpslevel

                    match ishappy with
                        | true -> // We're happy! :)
                            // Alias the discarded coefficients
                            coeffs <- chebtech2.alias( coeffs , cutoff )
                        | false -> 
                            // Replace any NaNs or Infs we may have extrapolated:
                            values <- Array.map2( fun v n -> 
                                match n with
                                | false -> v 
                                | true -> Double.NaN ) values maskNaN
                            values <- Array.map2( fun v i -> 
                                match i  with
                                | false -> v 
                                | true -> Double.PositiveInfinity ) values maskInf

        // Update the vscale.
        // Compute the 'true' vscale (as defined in CHEBTECH classdef):
        let vscaleOut = Array.max (Array.map( fun (x:double) -> Math.Abs(x) ) values)

        // Update vertical scale one last time:
        let vscaleGlobal = Math.Max( vscale , vscaleOut )

        // Output the 'true' vscale (i.e., the max of the stored values):
        vscale <- vscaleOut

        // Adjust the epslevel appropriately:
        epslevel <- epslevel * Math.Max(vscaleGlobal , epslevel) / Math.Max(vscaleOut , epslevel)

        // Assign to CHEBTECH object
        let mutable f = chebtech2( coeffs , vscale , hscale , ishappy , epslevel )

        // Output
        if ishappy then
            let newCoeffs , newEpslevel = chebtech2.simplify( f , f.epslevel / 100.0 )
            coeffs <- newCoeffs
            epslevel <- newEpslevel

        chebtech2( coeffs , f.vscale , f.hscale , f.ishappy , epslevel )

type mapping( f : double -> double , d : double -> double , i : double -> double ) = 
    // Note, we must use CAPS otherwise we cannot have 'for' as a member.
    // Forward map:
    member this.For = f
    // Derivative of the map:
    member this.Der = d
    // Inverse of the map:
    member this.Inv = i

    static member linear( dom: double array ) = 
        let a = dom.[0]
        let b = dom.[1]
        let For (y:double) = b*(y + 1.0)/2.0 + a*(1.0 - y)/2.0
        let Der (y:double) = (b - a)/2.0 + 0.0*y
        let Inv (x:double) = (x - a)/(b - a) - (b - x)/(b - a)
        mapping( For , Der , Inv )

    // empty constructor
    new() = mapping( (fun (double) -> (double)) , (fun (double) -> (double)) , (fun (double) -> (double)) )
    
// skip all these layers of abstraction for now:
// bndfun < classicfun < fun
type bndfun(d : double array , m : mapping , o : chebtech2) =
    // members
    member this.domain = d
    member this.mapping = m
    member this.onefun = o
    member this.epslevel with get() = this.onefun.epslevel
    member this.vscale with get() = this.onefun.vscale
    member this.hscale with get() = this.onefun.hscale
    static member createMap( dom: double array) = 
        mapping.linear( dom )
    
    static member feval( f:bndfun , x:double array ) =
        // Map the input:
        let z = Array.map f.mapping.Inv x
        // Evaluate the onefun:
        chebtech2.feval( f.onefun , z )

    // empty constructor
    new() = bndfun( Array.empty<double> , mapping(), chebtech2() )

    // anon-and-domain constructor
    new( op : double -> double , dom: double array ) = 

        // check the domain input
        match dom with
            | d when (not(d.Length.Equals(2)) || (d.[1] - d.[0]) <= 0.0) 
                -> failwith "Domain argument should be a row vector with two entries in increasing order."
            | d when (d |> Array.exists( fun x -> Double.IsInfinity(x) )) 
                -> failwith "Should not encounter unbounded domain in bndfun class."
            | _ -> ()
        
        // Remap the OP to be a function on [-1, 1].
        let map = bndfun.createMap( dom )
        let scaledOp = 
            match dom with
            | [|-1.0;1.0|] -> op
            | _ -> fun x -> op( map.For(x) )
        
        let onefun = chebtech2( scaledOp )
        
        bndfun( dom , map , onefun ) 

type chebfun( f: bndfun array , d: double array ) =
    // members
    member this.domain = d
    member this.funs = f
    member this.vscale 
        with get() = 
            let localVscale( f:bndfun ) = f.vscale
            // Get the maximum of the local vscales:
            Array.max (Array.map localVscale this.funs)
    member this.hscale with get() = this.domain.[1] - this.domain.[0]

    static member feval( f:chebfun , x:double array ) =
        let numFuns = f.funs.Length
        if ( numFuns = 1 ) then
            // Things are simple when there is only a single FUN:
            bndfun.feval( f.funs.[0] , x )
        else
            let mutable out = Array.zeroCreate<double> x.Length
            // For multiple FUNs we must determine which FUN corresponds to each x.
            // Replace the first and last domain entries with +/-inf. (Since we want to
            // use FUN{1} if real(x) < dom(1) and FUN{end} if real(x) > dom(end)).
            let domInf = f.domain
            domInf.[0] <- Double.NegativeInfinity
            domInf.[domInf.Length-1] <- Double.PositiveInfinity
            for k in [| 0 .. numFuns-1 |] do
                let isInThisDom = fun x -> x >= domInf.[k] && x < domInf.[k+1]
                let I = Array.map isInThisDom x
                let theseX = Array.filter isInThisDom x
                // Evaluate the appropriate fun 
                out <- Array.mapi (fun i o -> 
                    match I.[i] with
                    | true -> bndfun.feval( f.funs.[k] , [|x.[i]|] ).[0]
                    | false -> o) out
            out

    // empty contructor
    new() = chebfun( Array.empty<bndfun> , Array.empty<double> )

    // anon-and-domain constructor
    new( op: double -> double , dom: double array ) = 

        let numIntervals = dom.Length-1
        let ends = dom
        let funs = Array.zeroCreate<bndfun> numIntervals

        // Call the FUN constructor
        let getFun( op : double -> double , dom: double array ) =
            // Call the FUN constructor (which is the bndfun
            // consutctor before we support unbndfuns or deltafuns):
            let g = bndfun( op, dom )

            // See if the construction was happy:
            let ishappy = g.onefun.ishappy
            ( g, ishappy )

        for k in [| 0 .. numIntervals-1 |] do
            let endsk = ends.[k..k+1]
            let funk , ishappy = getFun( op, endsk)
            funs.[k] <- funk
            if not(ishappy) then
                failwith ("Function not resolved using " + funk.onefun.length.ToString() + " pts.")
        chebfun( funs , dom )

    // anon-only constructor
    new( op : double -> double ) = chebfun( op, [|-1.0;1.0|] )