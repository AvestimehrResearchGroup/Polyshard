Codinglib is mostly deprecated
======

**NOTE:** I'm deprecating Codinglib. The [ACTIS project][actis] for improving native coding theory support in Sage has now successfully completed, and many of the base elements of Codinglib are now better supported in the most recent Sage release.

Some elements of Codinglib are not yet ported into Sage, most notably the
alternative implementations of Guruswami--Sudan list-decoding and the
implementation of Power decoding. Note however that native Sage now supports one
reasonably fast implementation of the Guruswami--Sudan algorithm.

Codinglib is kept (mostly frozen) for reference and for compatibility with old
implementations. At time of writing (jan 2017) Codinglib still works in the
newest release of Sage; it's components are just largely incompatible with the
native Sage ones. Several of my past papers have proof-of-concept
implementations in separate repositories which use subroutines from Codinglib; I
will likely not update these.

I might later start a successor project of Codinglib, building on the native
coding theory structure in Sage. If this materialises, I will make this an
standard Sage package which should ease installation considerably.

Johan Rosenkilde (né Nielsen)


Previous Introduction
=========

This is **Codinglib**, a library for experimenting with algebraic coding theory.  
Copyright Johan S. H. Rosenkilde  
This is free software; every part of this released under the GNU
Public License version 3 or later (at your option).

The library is an extension to the Open Source computer algebra system Sage. I
have developed and maintain Codinglib for my own research on algebraic coding
theory. Most of the code is quite specific for this area.

Codinglib is partially self-documenting by the use doc-strings for every
function; however, it probably takes some effort to properly get to know it. It
is also quite idiosyncratic to my specific needs and workflow. Finally, in any
new version, I might change the name, semantics or calling convention of any
existing functions.

However, if you find some of its functionality interesting, I will gladly answer
specific questions to its usage and source code. You are also highly encouraged
to continue working on and improve any of the code; if you do, I would very much
like to know of your work. Also, I would be more than happy if any of my code
ended in the actual Sage distribution; I would work on this myself if I had more
time.


Overview
======
The following is a short description of each of the modules in Codinglib:

-  **bma**: The Berlekamp-Massey algorithm.

-  **channel**: Abstract class for representing channels

-  **chinese_remainder**: Chinese remainder codes: poly-alphabetic codes over `ZZ/pZZ` for different
	primes `p`, which are the `ZZ`-equivalent of Reed-Solomon codes. Note that this
	code implementation does not inherit from BlockCodeAbstract since it is
	poly-alphabetic. Also implementation of InterleavedChineseRemainderCode

-  **code**: Abstract and concrete class for representing linear block codes

-  **codeTesting**: Helper functions for testing error correcting code functions

-  **counting**: Wrapper-algebra for counting no. of operations in fields

-  **cyclic**: Cyclic linear codes.

-  **ea**: The Extended Euclidean algorithm with stopping conditions.
	Note that using modules over F[x] and weak Popov form is a more flexible
	approach.

-  **generate_readme**: Generates the package overview in the readme

-  **goppa**: Class for representing Goppa codes with corresponding utility
	functions

-  **graphcode**: Class for representing bipartite graph codes with corresponding utility
	functions

-  **gs**: The Guruswami-Sudan algorithm for Reed-Solomon codes: implementation and utilities

-  **gsKV**: Utilities for the variant of the Guruswami-Sudan algorithm using the small
	field Kotter-Vardy multiplicity assignment.

-  **hermitian_code**: Class for representing and decoding Hermitian codes.
	
	The codes follow the definition of [jsrn], which is a slightly larger class of
	codes than in [Dec].
	The decoding algorithms are described in [Dec].
	
	- [jsrn]: Johan S. R. Nielsen, List decoding of algebraic codes, 2013 PhD thesis
	          at Technical University of Denmark.
	- [Dec]:  Johan S. R. Nielsen, Peter Beelen, Sub-quadratic Decoding of One-Point
	          Hermitian Codes, in review for IEEE Transactions of Information Theory
	

-  **key2d**: 2D key equations and solving them.
	
	The name 2D Key Equations was introduced in [jsrn]. They are in Computer Algebra
	sometimes known as Simultaneous Hermitian Pade approximations.
	
	The notation here follows mainly that of [jsrn], but uses an improved handling
	of weights described in [Dec].
	
	- [jsrn]: Johan S. R. Nielsen, List decoding of algebraic codes, 2013 PhD thesis
	          at Technical University of Denmark.
	- [Dec]:  Johan S. R. Nielsen, Peter Beelen, Sub-quadratic Decoding of One-Point
	          Hermitian Codes, in review for IEEE Transactions of Information Theory
	

-  **listdecoding**: Utility functions related to list decoding

-  **module**: Functionality for free F[x]-modules, in particular for calculating
	bases in weak Popov form. The weak Popov form is slightly stronger
	than row reduced form, and it is also a type of reduced Groebner bases
	over some simple module monomial orderings.
	
	- [Dec]:  Johan S. R. Nielsen, Peter Beelen, Sub-quadratic Decoding of One-Point
	          Hermitian Codes, in review for IEEE Transactions of Information Theory
	

-  **parallel**: Module for a simple thread pool, along with invoking workers over SSH.
	

-  **plot_util**: Utility functions for plotting. Contains an (now old) copy of the
	implementation of points for plotting, in order to allow custom marker symbols
	in list plots.
	
	Ideally, this should be reported and fixed upstream in Sage rather than be here.
	

-  **polynomial_util**: Generally usable polynomial functions

-  **ratinter**: Rational interpolation using the Wu algorithm

-  **remote**: Functions for importing Codinglib remotely, loaded on-the-fly
	directly from the repository. This is useful for loading Codinglib on
	a public Sage server.

-  **rootfinding**: F[x] Root-finding in polynomials over F[x,y].
	

-  **rs**: Class for representing GRS codes with corresponding utility functions

-  **rs_mindist**: Implementation of decoding algorithms Reed-Solomon up to half the
	minimum distance

-  **rs_power**: Functionality for Power decoding of Reed-Solomon codes.
	Power decoding was called "Decoding by virtual extension to an interleaved
	code" in the initial description of [SSB].
	
	In the initial algorithm of [SSB], the decoding radius is bounded by (roughly)
	the Sudan radius. In [N15], it was shown how to incorporate multiplicities to
	decode up to the Johnson radius. The implementation here supports that.
	
	REFERENCES:
	 - Schmidt, G., V.R. Sidorenko, and M. Bossert. "Syndrome Decoding of
	   Reed-Solomon Codes Beyond Half the Minimum Distance Based on Shift-Register
	   Synthesis." IEEE Transactions on Information Theory 56, no. 10 (2010):
	   5245-52. doi:10.1109/TIT.2010.2060130.
	
	 - Nielsen, Johan S. R. "Power Decoding Reed-Solomon Up to the Johnson
	   Radius.". Submitted to IEEE Transactions of Information Theory, 2015.
	

-  **subcode**: Direct construction of subfield subcodes by expansion of parity check matrices

-  **util**: Various utility functions used in codinglib

-  **wu**: Wu list decoder for Reed-Solomon and Binary Goppa codes.
	Parameter choices. See also ratinter for the core rational
	interpolation functionality.


Furthermore, the library contains a number of `sheet` files: these are
similar to the notebook sheets, in that they contain snippets of code
demonstrating or testing the core functionality of Codinglib. Most of
these sheets are indeed intended as tests.

For instance, `rs_decoding.sheet` demonstrates all decoding algorithms
for Reed-Solomon codes included in Codinglib.

For Emacs users, sheet files can elegantly be handled when using
[sage-mode] by using the `sage-blocks` functionality.



Usage
======
Sage does not currently work well with non-spkg library code, so the importing
of Codinglib is a bit peculiar.

Codinglib is written using the Sage language extensions to Python but at the same
time is highly dependent on itself. To make this work in Sage currently, one
needs to preprocess Codinglib's .sage files and load all the generated .py files
at once as a Python package. In particular, using load/attach on the
.sage files from a Sage prompt will usually not work.

Method 1: Downloading to personal computer
-----------
-   Download the whole source (mirror git repo) and put in a folder `codinglib`.
    Add the parent folder to the SAGE_PATH environment variable.

    For example, you have created a folder

        \home\foo\bar\codinglib

    on your system, which contains the `.sage` files and the rest of Codinglib.
    You then add 

        \home\foo\bar

    to your `SAGE_PATH`.

-   Navigate to the `codinglib` folder and type the command

           make

    This will preprocess all the `.sage` files to reduce Sage syntax to standard
    Python.

-   From the Sage prompt, the library can be loaded with the command

        import codinglib 


-   From the Sage Notebook, if running on a local Sage server, add the following
    to one of the first cells of each sheet in which you intend to use the
    library:

        #auto
        import codinglib 
        all = __builtins__.all

    The last line is to restore the Python standard library function `all`. 

    In either of the above, you can of course replace `import codinglib` with the
    following to include commands directly in the namespace

        from codinglib import *



Method 2: On SageMathCloud
-----------
The SageMathCloud has excellent support for making modules within a project
which can be run from any worksheet in that project.

To add Codinglib to a project, you should make a copy of Codinglib's source in a
folder `codinglib` within your project. The most direct way to do this is:

- Create a folder `codinglib` in the project on SageMathCloud and go to it's Add
    Files page.

- Download Codinglib's source to your own machine.

- Open a file explorer and drag-and-drop all the downloaded files onto the Drop
Files-area. files into the SageMathCloud's add-page.

- In SageMathCloud, go to the `codinglib`-folder's page. In the upper right-hand
    corner is a "Terminal command" text area. Enter

        make build

    This will create `.py` files for each of the `.sage` files.

    SageMathCloud sometimes complains that the command timed out after 15
    seconds and was halted. This might mean nothing, or it might mean that some
    of the `.py`-files were not properly generated. In that case, simply run
    `make build` again.
    

- You should now be able to run

        from codinglib import *

    from any worksheet in this project.


If you later need to update Codinglib, you simply need to overwrite the `.sage.`
files with the new ones and rerun the `make` command.


Method 3: On a shared Sage Notebook sheet without write-access to its server
-----------
This is much more difficult as the Sage Notebook is not really geared to
support non-spkg add-on Sage libraries. The following method basically
dynamically downloads, preprocesses and imports Codinglib whenever one starts the
Notebook sheet. Thus, it is very slow to start up and a bit tedious. However, it
works without assuming any special rights on the Sage server. Once loaded, it
will of course run in the usual speed of the server.

In one of the first cells of the Notebook sheet, add the following code

    changeset = "ddd6636b2da0"
    baseurl ="https://bitbucket.org/jsrn/codinglib/raw/" + changeset + "/"

    # Manually load the remote functionality
    import imp
    from sage.misc.remote_file import get_remote_file
    remoteFile = get_remote_file(baseurl + "remote.py", verbose=False)
    remote = imp.new_module("remote")
    exec open(remoteFile).read() in remote.__dict__

    # Load the Sage library files and import them into the global namespace
    imports = remote.retrieve_codinglib(baseurl)
    for stm in imports:
        exec(stm)

Replace the string in changeset with the Git revision no. of the wanted version
of Codinglib, as provided by BitBucket. For the most recent version, visit

    https://bitbucket.org/jsrn/codinglib

and copy the revision no. of the top commit.

Executing this cell takes 10-20 seconds depending on your server's internet
connection.



Regards,  
Johan S. H. Rosenkilde__
jsrn@jsrn.dk  
[jsrn.dk]

[jsrn.dk]: http://jsrn.dk
[sage-mode]: https://bitbucket.org/gvol/sage-mode
[actis]: https://bitbucket.org/lucasdavid/sage_coding_project/wiki/Home
