==========
User guide
==========

Installation
------------

The framol package can currently only be installed from source.
The steps are:

.. code-block::

	git clone https://github.com/saraidery/fragment-molecule
	cd fragment-molecule
	pip install --editable .


Basic functionality
-------------------

The basic functionality of framol is the fragmentation of a molecular system.
This is done using the ``MolecularFragmenter``.

For example, we can fragment a DNA strand (given in the file ``dna_strand.xyz``), with fragments that do not exceed 50 atoms, by the following lines:

.. code-block:: python

	from framol import MolecularFragmenter

	f = MolecularFragmenter(50, "dna_strand.xyz")
	f.fragment()

We can visualize the result of the fragmentation by a call to the ``plot_fragments`` method. Either using the default (`CPK <https://en.wikipedia.org/wiki/CPK_coloring>`_) colors

.. code-block:: python

	f.plot_fragments()

or by giving each fragment a random color

.. code-block:: python

	f.plot_fragments(color="random")

The fragment geometries can be stored to XYZ-files by calling the ``store_fragments`` method:

.. code-block:: python

	f.store_fragments("dna_strand")

The files ``dna_strand_fragment_i.xyz`` are then generated for :math:`i=0:n - 1` where :math:`n`
is the number of fragments.


Properties of the fragmented system
-----------------------------------

The size of the molecular system to be fragmented is ``f.size`` and the size of the individual fragments is ``f.fragment_sizes``. The number of fragments is ``f.n_fragments`` and the number of capped bonds---i.e., how many bonds there are between different fragments---is ``f.n_capped_bonds``.

.. note::

	The MolecularFragmenter has no fragments before its ``fragment`` method is called.


Advanced functionality
----------------------

The framol package can add hydrogen atoms at the ends of capped bonds through:

.. code-block:: python

	f.add_H_to_capped_bonds()
	f.store_fragments("dna_strand")

This will result in addition of H atoms along the bond that was capped.


.. code-block:: python

	i = f.find_central_fragment()
	f.swap_fragments(0, i)
	f.group_fragments_by_size()
	f.store_fragments("dna_strand")
