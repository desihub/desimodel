=================
Testing desimodel
=================

Introduction
------------

Testing desimodel is a bit tricky since most of the code involves reading
data files that are not included with the git product.  In addition, the data
files can be rather large, both individually and as a package. This document
describes how to create lightweight test branches of the desimodel *data*
that can be used for rapid testing.
