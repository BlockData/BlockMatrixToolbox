# BlockMatrixToolbox
Implementation of multi-blocks data structure (matrices) for Matlab

This toolbox for Matlab allows creation and manipulation of multi-block matrices, also called
"block-matrices".
Block-matrices are matrices endowed with a "block-dimensions", that defines a decomposition of 
the matrix data into several blocks. 

The motivation is to manage efficiently concatenation of several data sets while preserving the 
original structure of the experiment.

# Overview

The toolbox can be downloaded by cloning the project from the following URL:
[http://github.com/BlockData/BlockMatrixToolbox](http://github.com/BlockData/BlockMatrixToolbox). Select the "release" tab to download archive, or simply cone the projet on your local computer.
Then, simply add the path to the "BlockMatrixToolbox/multiBlock" directory to the path.

## Creation of Block-Matrices

Let us start with a rather simple matrix:

    >> data = reshape(1:28, [7 4])';
    >> data
    data =
         1     2     3     4     5     6     7
         8     9    10    11    12    13    14
        15    16    17    18    19    20    21
        22    23    24    25    26    27    28

A Block-Dimension can be assigned to this matrix, using the **BlockMatrix** class:

    >> X = BlockMatrix(data, [2 2], [2 3 2]);
 
As an alternative syntax, the **BlockDimension** object can be explicitely created and given as argument:

    >> BD = BlockDimension({[2 2], [2 3 2]});
    >> X = BlockMatrix(data, BD);

The content of a BlockMatrix can be shown by using the `disp` function, or simply by omitting the 
final semi colon:

	>> X
    X = 
    BlockMatrix object with 4 rows and 7 columns
      row dims: 2 2
      col dims: 2 3 2
            1           2           3           4           5           6           7   
            8           9          10          11          12          13          14   
           15          16          17          18          19          20          21   
           22          23          24          25          26          27          28   

## Operations on Block-Matrices

Block-Matrices can be transposed, and multiplied. The block-structure of the resulting 
block-matrices are preserved.

	>> X * X'
	ans =
	BlockMatrix object with 4 rows and 4 columns
    row dims: 2 2
    col dims: 2 2
          140         336         532         728   
          336         875        1414        1953   
          532        1414        2296        3178   
          728        1953        3178        4403   

The content of Block-matrices is accessible by different ways. One possibility is to use the 
`getBlock` and `setBlock` methods.

	>> getBlock(X, 1, 2)
	ans =
		 3     4     5
		10    11    12

	>> setBlock(X, 1, 2, ones(2, 3));
	>> X
	X = 
	BlockMatrix object with 4 rows and 7 columns
	  row dims: 2 2
	  col dims: 2 3 2
			1           2           1           1           1           6           7   
			8           9           1           1           1          13          14   
		   15          16          17          18          19          20          21   
		   22          23          24          25          26          27          28   
		
A similar result is obtained by using braces indexing:

	>> X{2, 3}
	ans =
		20    21
		27	  28

		
Elements of a block matrix can be accessed using parens indexing

	>> X(3, 2)
	ans =
		16

## Block-Diagonal matrices

The toolbox provides support for **Block-Diagonal matrices**. Block-Diagonal matrices are block-matrices
with the particularity that the blocks that are not on the diagonal contains only zero values.
An example of block-diagonal matrix is as follow:

	>> BD = BlockDiagonal({[1 2 3; 4 5 6], [7 8;9 10], [11 12]})
	BD =
	BlockMatrix object with 5 rows and 7 columns
	  row dims: 2 2 1
	  col dims: 3 2 2
			1           2           3           0           0           0           0   
			4           5           6           0           0           0           0   
			0           0           0           7           8           0           0   
			0           0           0           9          10           0           0   
			0           0           0           0           0          11          12   

BlockDiagonal matrices can be multiplied as common block-matrices:

	>> BD * BD'
	ans = 
	BlockMatrix object with 5 rows and 5 columns
	  row dims: 2 2 1
	  col dims: 2 2 1
		   14          32           0           0           0   
		   32          77           0           0           0   
			0           0         113         143           0   
			0           0         143         181           0   
			0           0           0           0         265   

# Credits

**BlockMatrixToolbox** is developped jointly by INRA and ONIRIS. The development started at the occasion of 
the "CouplImSpec" project, involving also the SOLEIL synchrotron, and during the AI-Fruit project
of Region "Pays-de-la-Loire"

Authors:

* Mohamed Hanafi (oniris, mohamed.hanafi, at oniris-nantes.fr)
* David Legland (INRA, BIBS Platform, david.legland, at nantes.inra.fr)

See also:

* [http://github.com/gastonstat/blockberry](http://github.com/gastonstat/blockberry)
The R implementation of multi-block matrices that inspired this package
