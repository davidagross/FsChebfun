namespace Chebfun

// skip all these layers of abstraction for now:
// chebtech1 < chebtech < smoothfun < onefun
type chebtech1( c: double list, v: double, h: double, i: bool, e: double) = 
    member this.coeffs = c
    member this.vscale = v
    member this.hscale = h
    member this.ishappy = i
    member this.epslevel = e
    new() = chebtech1( List.empty<double> , 0.0, 1.0, new bool(), new double() )

type mapping( f : 'a -> 'a, d : 'a -> 'a, i : 'a -> 'a ) = 
    // Note, we must use CAPS otherwise we cannot have 'for' as a member.
    // Forward map:
    member this.For = f
    // Derivative of the map:
    member this.Der = d
    // Inverse of the map:
    member this.Inv = i
    new() = mapping( (fun (double) -> (double)) , (fun (double) -> (double)) , (fun (double) -> (double)) )
    

// skip all these layers of abstraction for now:
// bndfun < classicfun < fun
type bndfun(d : double list , m : mapping , o : chebtech1) =
    member this.domain = d
    member this.mapping = m
    member this.onefun = o
    new() = bndfun( List.empty<double> , mapping(), chebtech1() )

type chebfun( d: double list , b: bndfun list) =
    member this.domain = d
    member this.funs = b
    new() = chebfun( List.empty<double> , List.empty<bndfun> )
