Athena
--------

Athena is a read cloud assembler for large low-copy repeat families within
mammalian genomes. 

This repository and setup material is currently under construction.  It
will soon be updated with easy to use instructions.

Installation
============

Athena depends on the `phaser <https://github.com/abishara/phaser` library, which must first be
installed.

The following prerequisites must also be installed:

* `idba_ud <https://github.com/abishara/idba/releases/tag/1.1.3a1>`_ -- please use **this** version, which is modified both to handle longer short-read lengths and to locally assemble subsampled barcoded reads clouds.  Ensure all compiled binaries, including ``idba_subasm``, are in your ``$PATH``


Troubleshooting
===============

The ``athena-meta`` command may be run multiple times to resume the pipeline.

If an error arises, the output from ``athena-meta`` or the log files may
be informative.

**ShortSequence: Sequence is too long.** If you get this error during
assembly, please make sure you are using `the right fork of idba_ud
<https://github.com/abishara/idba/releases/tag/1.1.3a1>`_.

Please submit issues on the `github page for Athena
<https://github.com/abishara/athena_meta/issues>`_.

