Installation
============

You can install GeneBe Utils using pip:

.. code-block:: bash

   pip install genebe

If you wish to install a faster `mmh3` implementation or use the option of annotating VCF files, install using:

.. code-block:: bash

   pip install genebe[cpp]

or install modules separately:

.. code-block:: bash

   pip install cyvcf2
   pip install mmh3

This step will require build tools installed on your computer.

Docker
------

There is a dockerized version of this package, available at https://hub.docker.com/r/genebe/pygenebe.

Usage example, reading from file `input.vcf` and writing output to `stdout`:

.. code-block:: bash

   docker run -v input.vcf:/tmp/input.vcf --rm genebe/pygenebe:0.0.14 genebe annotate --input /tmp/input.vcf --output /dev/stdout

Limits
------

If you wish to annotate thousands of variants, please log in to https://genebe.net, generate an API Key, and provide it using `--username` and `--api-key`.

The number of daily requests from a single IP is restricted to prevent abuse and excessive resource consumption on our server. Account holders with an API Key enjoy significantly higher limits (in the tens of thou
