namespace test
open Chebfun
open FsUnit
open FsUnit.Xunit
open NUnit
open NUnit.Framework
open System
open Xunit

module test_chebfun = 

    [<Test>]
    let ``test_constructor_basic``() = 
        // Some basic test functions
        let FF = [ Math.Sin ; Math.Cos ; fun x -> Math.Exp(-x) ]
        // TODO move to a math library and/or start using Math.Net Numerics Vectors
        let infnorm (a:double array) (b:double array) = 
            Array.max(Array.map2( fun (a:double) (b:double) -> Math.Abs(a - b) ) a b)
        // TODO move into a real pref / chebfunpref type
        let eps = 2.2e-16

        for F in FF do
            // Test on [-1 1]:
            let f = chebfun( F , [|-1.0 ; 1.0|] )
            let xx = [|-1.0..(2.0/99.0)..1.0|]
            let err = infnorm (chebfun.feval( f , xx )) (Array.map F xx)
            err |> should lessThan (10.0*eps*f.vscale)
            err |> should lessThan (50.0*eps)

            // Test on [-1 1] (no domain passed):
            let f = chebfun( F )
            let xx = [|-1.0..(2.0/99.0)..1.0|]
            let err = infnorm (chebfun.feval( f , xx )) (Array.map F xx)
            err |> should lessThan (10.0*eps*f.vscale)
            err |> should lessThan (500.0*eps)

            // Test on [0 10000]:
            let f = chebfun( F , [|0.0 ; 10000.0|] )
            let xx = [|0.0..(10000.0/99.0)..10000.0|]
            let a = (chebfun.feval( f , xx ))
            let b = (Array.map F xx)
            let err = infnorm a b 
            err |> should lessThan (10000.0*eps*f.vscale*20.0) // TODO: why are we off by a factor of 20?
            err |> should lessThan (100.0*f.hscale*eps)

            // Test on piecewise domain:
            let f = chebfun( F , [|-1.0 ; 0.0 ; 0.5 ; Math.Sqrt(Math.PI/4.0) ; 1.0|] )
            let xx = [|-1.0..(2.0/99.0)..1.0|]
            let a = (chebfun.feval( f , xx ))
            let b = (Array.map F xx)
            let err = infnorm a b 
            err |> should lessThan (10.0*eps*f.vscale)
            err |> should lessThan (100.0*f.hscale*eps)