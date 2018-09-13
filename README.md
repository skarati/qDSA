qDSA using Kummer Line:
-----------------------
This repository contains implementations of qDSA algorithms based 
on the Kummer Line from https://github.com/skarati/KummerLineV02.

There are three qDSA algorithms:
   1. qDSA_KL2519  : qDSA based on KL2519,
   2. qDSA_KL25519 : qDSA based on KL25519, and
   3. qDSA_KL2663  : qDSA based on KL2663.

Each directory contains the file qDSA.h which includes the Key Generation,
Signing and Verification algorithms. The scalar.c and scalar.h are modified
versions of the files which are available at the hompage of Joost Renes. The
test.c file used to check and measure the software. The file run.sh contains
the compilation flags used for the code.

Implementor:
------------
Sabyasachi Karati
