namespace test
open System
open Microsoft.VisualStudio.TestTools.UnitTesting
open Chebfun

[<TestClass>]
type test_chebfun() = 

    [<TestMethod>]
    member this.test_constructor_basic() = 
        // Some basic test functions
        let FF = [ Math.Sin ; Math.Cos ; fun x -> Math.Exp(-x) ]
        // TODO move to a math library and/or start using Math.Net Numerics Vectors
        let infnorm (a:double list) (b:double list) = 
            List.max (List.map2( fun (a:double) (b:double) -> Math.Abs(a - b) ) a b)
        // TODO move into a real pref / chebfunpref type
        let eps = 2.2e-16

        for F in FF do
            // Test on [-1 1]:
            let f = chebfun( F , [-1.0 ; 1.0] )
            let xx = [-1.0..(1.0/99.0)..1.0]
            let err = infnorm (chebfun.feval( f , xx )) (List.map F xx)
            Assert.IsTrue( err < 10.0*f.epslevel*f.vscale )
            Assert.IsTrue( err < 50.0*eps )

            // Test on [-1 1] (no domain passed):
            let f = chebfun( F )
            let xx = [-1.0..(1.0/99.0)..1.0]
            let err = infnorm (chebfun.feval( f , xx )) (List.map F xx)
            Assert.IsTrue( err < 10.0*f.epslevel*f.vscale )
            Assert.IsTrue( err < 500.0*eps )

            // Test on [0 10000]:
            let f = chebfun( F , [0.0 ; 10000.0] )
            let xx = [0.0..(10000.0/99.0)..10000.0]
            let a = (chebfun.feval( f , xx ))
            let b = (List.map F xx)
            let err = infnorm a b 
            Assert.IsTrue( err < 100.0*f.epslevel*f.vscale )
            Assert.IsTrue( err < 100.0*f.hscale*eps )

            // TODO Test on piecewise domain: