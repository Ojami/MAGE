%% About defineqfc.mlx
% This file defines the MATLAB interface to the library |qfc|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for &lt;SHAPE&gt;, &lt;DIRECTION&gt;, etc. For more
% information, see <matlab:helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') Define MATLAB Interface for C++ Library>.



%% Setup. Do not edit this section.
function libDef = defineqfc()
libDef = clibgen.LibraryDefinition("qfcData.xml");
%% OutputFolder and Libraries 
libDef.OutputFolder = "C:\Users\xjamov\Desktop\test\qfc";
libDef.Libraries = "";

%% C++ function |exp1| with MATLAB name |clib.qfc.exp1|
% C++ Signature: double exp1(double x)
exp1Definition = addFunction(libDef, ...
    "double exp1(double x)", ...
    "MATLABName", "clib.qfc.exp1", ...
    "Description", "clib.qfc.exp1    Representation of C++ function exp1."); % Modify help description values as needed.
defineArgument(exp1Definition, "x", "double");
defineOutput(exp1Definition, "RetVal", "double");
validate(exp1Definition);

%% C++ function |counter| with MATLAB name |clib.qfc.counter|
% C++ Signature: void counter()
counterDefinition = addFunction(libDef, ...
    "void counter()", ...
    "MATLABName", "clib.qfc.counter", ...
    "Description", "clib.qfc.counter    Representation of C++ function counter."); % Modify help description values as needed.
validate(counterDefinition);

%% C++ function |square| with MATLAB name |clib.qfc.square|
% C++ Signature: double square(double x)
squareDefinition = addFunction(libDef, ...
    "double square(double x)", ...
    "MATLABName", "clib.qfc.square", ...
    "Description", "clib.qfc.square    Representation of C++ function square."); % Modify help description values as needed.
defineArgument(squareDefinition, "x", "double");
defineOutput(squareDefinition, "RetVal", "double");
validate(squareDefinition);

%% C++ function |cube| with MATLAB name |clib.qfc.cube|
% C++ Signature: double cube(double x)
cubeDefinition = addFunction(libDef, ...
    "double cube(double x)", ...
    "MATLABName", "clib.qfc.cube", ...
    "Description", "clib.qfc.cube    Representation of C++ function cube."); % Modify help description values as needed.
defineArgument(cubeDefinition, "x", "double");
defineOutput(cubeDefinition, "RetVal", "double");
validate(cubeDefinition);

%% C++ function |log1| with MATLAB name |clib.qfc.log1|
% C++ Signature: double log1(double x,BOOL first)
log1Definition = addFunction(libDef, ...
    "double log1(double x,BOOL first)", ...
    "MATLABName", "clib.qfc.log1", ...
    "Description", "clib.qfc.log1    Representation of C++ function log1."); % Modify help description values as needed.
defineArgument(log1Definition, "x", "double");
defineArgument(log1Definition, "first", "int32");
defineOutput(log1Definition, "RetVal", "double");
validate(log1Definition);

%% C++ function |order| with MATLAB name |clib.qfc.order|
% C++ Signature: void order()
orderDefinition = addFunction(libDef, ...
    "void order()", ...
    "MATLABName", "clib.qfc.order", ...
    "Description", "clib.qfc.order    Representation of C++ function order."); % Modify help description values as needed.
validate(orderDefinition);

%% C++ function |errbd| with MATLAB name |clib.qfc.errbd|
% C++ Signature: double errbd(double u,double * cx)
%errbdDefinition = addFunction(libDef, ...
%    "double errbd(double u,double * cx)", ...
%    "MATLABName", "clib.qfc.errbd", ...
%    "Description", "clib.qfc.errbd    Representation of C++ function errbd."); % Modify help description values as needed.
%defineArgument(errbdDefinition, "u", "double");
%defineArgument(errbdDefinition, "cx", "clib.array.qfc.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Double", or "double"
%defineOutput(errbdDefinition, "RetVal", "double");
%validate(errbdDefinition);

%% C++ function |ctff| with MATLAB name |clib.qfc.ctff|
% C++ Signature: double ctff(double accx,double * upn)
%ctffDefinition = addFunction(libDef, ...
%    "double ctff(double accx,double * upn)", ...
%    "MATLABName", "clib.qfc.ctff", ...
%    "Description", "clib.qfc.ctff    Representation of C++ function ctff."); % Modify help description values as needed.
%defineArgument(ctffDefinition, "accx", "double");
%defineArgument(ctffDefinition, "upn", "clib.array.qfc.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Double", or "double"
%defineOutput(ctffDefinition, "RetVal", "double");
%validate(ctffDefinition);

%% C++ function |truncation| with MATLAB name |clib.qfc.truncation|
% C++ Signature: double truncation(double u,double tausq)
truncationDefinition = addFunction(libDef, ...
    "double truncation(double u,double tausq)", ...
    "MATLABName", "clib.qfc.truncation", ...
    "Description", "clib.qfc.truncation    Representation of C++ function truncation."); % Modify help description values as needed.
defineArgument(truncationDefinition, "u", "double");
defineArgument(truncationDefinition, "tausq", "double");
defineOutput(truncationDefinition, "RetVal", "double");
validate(truncationDefinition);

%% C++ function |findu| with MATLAB name |clib.qfc.findu|
% C++ Signature: void findu(double * utx,double accx)
%finduDefinition = addFunction(libDef, ...
%    "void findu(double * utx,double accx)", ...
%    "MATLABName", "clib.qfc.findu", ...
%    "Description", "clib.qfc.findu    Representation of C++ function findu."); % Modify help description values as needed.
%defineArgument(finduDefinition, "utx", "clib.array.qfc.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Double", or "double"
%defineArgument(finduDefinition, "accx", "double");
%validate(finduDefinition);

%% C++ function |integrate| with MATLAB name |clib.qfc.integrate|
% C++ Signature: void integrate(int nterm,double interv,double tausq,BOOL mainx)
integrateDefinition = addFunction(libDef, ...
    "void integrate(int nterm,double interv,double tausq,BOOL mainx)", ...
    "MATLABName", "clib.qfc.integrate", ...
    "Description", "clib.qfc.integrate    Representation of C++ function integrate."); % Modify help description values as needed.
defineArgument(integrateDefinition, "nterm", "int32");
defineArgument(integrateDefinition, "interv", "double");
defineArgument(integrateDefinition, "tausq", "double");
defineArgument(integrateDefinition, "mainx", "int32");
validate(integrateDefinition);

%% C++ function |cfe| with MATLAB name |clib.qfc.cfe|
% C++ Signature: double cfe(double x)
cfeDefinition = addFunction(libDef, ...
    "double cfe(double x)", ...
    "MATLABName", "clib.qfc.cfe", ...
    "Description", "clib.qfc.cfe    Representation of C++ function cfe."); % Modify help description values as needed.
defineArgument(cfeDefinition, "x", "double");
defineOutput(cfeDefinition, "RetVal", "double");
validate(cfeDefinition);

%% C++ function |qf| with MATLAB name |clib.qfc.qf|
% C++ Signature: double qf(double * lb1,double * nc1,int * n1,size_t r1,double sigma,double c1,int lim1,double acc,double * trace,int * ifault)
%qfDefinition = addFunction(libDef, ...
%    "double qf(double * lb1,double * nc1,int * n1,size_t r1,double sigma,double c1,int lim1,double acc,double * trace,int * ifault)", ...
%    "MATLABName", "clib.qfc.qf", ...
%    "Description", "clib.qfc.qf    Representation of C++ function qf."); % Modify help description values as needed.
%defineArgument(qfDefinition, "lb1", "clib.array.qfc.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Double", or "double"
%defineArgument(qfDefinition, "nc1", "clib.array.qfc.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Double", or "double"
%defineArgument(qfDefinition, "n1", "clib.array.qfc.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Int", or "int32"
%defineArgument(qfDefinition, "r1", "uint64");
%defineArgument(qfDefinition, "sigma", "double");
%defineArgument(qfDefinition, "c1", "double");
%defineArgument(qfDefinition, "lim1", "int32");
%defineArgument(qfDefinition, "acc", "double");
%defineArgument(qfDefinition, "trace", "clib.array.qfc.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Double", or "double"
%defineArgument(qfDefinition, "ifault", "clib.array.qfc.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.qfc.Int", or "int32"
%defineOutput(qfDefinition, "RetVal", "double");
%validate(qfDefinition);

%% Validate the library definition
validate(libDef);

end
