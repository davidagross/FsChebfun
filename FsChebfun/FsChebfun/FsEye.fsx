//reference FsEye.dll (e.g. "Send to F# Interactive" on the project reference availabe in some IDE versions)
#r @"..\packages\FsEye.2.1.0\lib\net40\FsEye.dll"

//Execute the following two lines of code to bring the eye singleton into scope and bind it to FSI
open Swensen.FsEye.Fsi //bring the eye singleton into scope
fsi.AddPrintTransformer eye.Listener //attached the listener