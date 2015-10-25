About
=====

FsChebfun is a port of [Chebfun][1], the open-source software system for numerical computing with functions, from MATLAB to F#.  Its current abilities are minimal, but are intended to grow to include, at least, inversion of monotonic functions, differentiation, and integration, allowing for some modest but still powerful statistical analyses.

Structure
=========

For now, a test-driven development paradigm has been used to add functionality to this project.  Currently, we have a passing `test_constructor_basic` from Chebfun by implementing the following types and their constructors in `chebfun.fs`.

* `chebtech2`
  - `chebpts`
  - `barywts`
  - `clenshaw`
  - `feval`
  - `refine`
  - `extrapolate`
  - `vals2coeffs`
  - `coeffs2vals`
  - `alias`
  - `happinessCheck`
  - `simplify`
* `mapping`
  - `linear`
* `bndfun`
  - `createMap`
  - `feval`
* `chebfun`
  - `feval`

Detailing the available codebase here should be discontinued once the structure of the source and its documentation is enhanced.

Usage
=====

The test methods in `test_chebfun.fs` detail some of the usage of the `chebfun` constructor.

``` f#
let f = chebfun( fun x -> Math.Exp(-x) , [|-1.0 ; 1.0|])
```

While function composition, max, plotting, and root-finding will not be primary priorities for this project, great power should still derived from `chebfun` creation from MATLAB-style "anonymous" function in F#, as above.  Creating a PDF from values on equispaced points and inverting their CDF is a long-term goal of this project.

License
=======

See `LICENSE.md` for FsChebfun's licensing information.

[1]: https://github.com/chebfun/chebfun
