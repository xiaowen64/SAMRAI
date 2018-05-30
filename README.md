# SAMRAI

SAMRAI (Structured Adaptive Mesh Refinement Application Infrastructure) is an
object-oriented C++ software library that enables exploration of numerical,
algorithmic, parallel computing, and software issues associated with applying
structured adaptive mesh refinement (SAMR) technology in large-scale parallel
application development. SAMRAI provides software tools for developing SAMR
applications that involve coupled physics models, sophisticated numerical
solution methods, and which require high-performance parallel computing
hardware. SAMRAI enables integration of SAMR technology into existing codes and
simplifies the exploration of SAMR methods in new application domains. 

## New Release

SAMRAI v. 2.13.0 is now available! The major change is that we no longer use
the Boost library. This affects many interfaces in the API that used
boost::shared_ptr, which has been replaced by std::shared_ptr.

## Get Involved

SAMRAI is an open source project, and questions, discussion and contributions
are welcome!

### Mailing List

To get in touch with all the SAMRAI developers, please email samrai@llnl.gov

### Contributions

Contributing to SAMRAI should be easy! We are managing contributions through
pull requents here on GitHub. When you create your pull request, please make
`develop` the target branch.

Your PR must pass all of SAMRAI's unit tests, which are enforced using Travis
CI. For information on how to run these tests locally, please see our
[contribution guidelines](CONTRIBUTING.md)

Our development process roughly follows the
[gitflow](http://nvie.com/posts/a-successful-git-branching-model/) model. The
`master` branch points to the last stable, tagged release. The `develop` branch
contains the latest changes. New features should be created in `feature/<name>`
branches and be based on `develop`.

## Citing SAMRAI

We maintain a list of publications
[here](https://computation.llnl.gov/projects/samrai/publications).

## Release

Copyright (c) 1997-2018, Lawrence Livermore National Security, LLC.

Produced at the Lawrence Livermore National Laboratory.

All rights reserved.

Released under LGPL v2.1

For release details and restrictions, please read the LICENSE file. It is also
linked here: [LICENSE](./LICENSE)

`LLNL-CODE-434871`
