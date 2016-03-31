// Learn more about F# at http://fsharp.net. See the 'F# Tutorial' project
// for more guidance on F# programming.

#load "../packages/MathNet.Numerics.FSharp.3.5.0/MathNet.Numerics.fsx"
#load "chebfun.fs"
// #r "bin/Debug/chebfun.dll"
open Chebfun
open System

let x = new chebfun( (fun x -> x*x) )

chebfun.feval(x , [|Math.PI/10.0|]) |> Array.map Math.Sqrt