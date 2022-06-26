# murk2

## Introduction
murk2 is a C++ toolkit for researching attacks on modern systems, particularly cryptographic ones.

This framework aims to satisfy the following objectives (in order of importance):

1. Useful: it must comprehensively cover real-world attacks on modern systems.
2. Reliable: it needs to give usable output, so that it can form part of a research paper or writeup.
3. Reproducable: any attack must be reproducable by that exact version, so that it can be cited (later versions can break this).
4. Usable: while the API will not be the cleanest, it should at least be navigable by a relatively new C++ programmer.
5. Modular: there should exist a way of stringing together theoretically compatible components.
6. Performant: it needs to be actually usable in a research setting, not just some "guaranteed to return" function.

## Motivation
While there is a lot of very high quality cryptanalysis out there, it is often stuck in the form of vague pseudocode,
or disparate unmaintained git repos. This hinders its applicability to real-world engagements or further research, 
which in turn reduces the security of the systems we use. 
murk2 aims to properly unify these attacks into one framework, where complex procedures can be modified and examined by researchers,
as well as directly applied to real-world scenarios by pentesters.

Similarly, many public CTF and pentest writeups describe complex bash pipelines, without giving reproducible results. 
The original murk (for all it's *many* faults) provided a clean and (reasonably) efficient API, as well as modular components, 
that allowed for the creation of workable complex attacks.
murk2 aims to keep this idea alive, by exposing interfaces to apply a theoretical attack to a real-world scenario.

Note that the aim is *not* to create another Metasploit, as Metasploit already does that pretty well. 
Instead, murk2 provides attack primitives, rather than ready-made attacks for existing software.
