classdef BlockDiagonalTest < matlab.unittest.TestCase
%BLOCKDIAGONALTEST Test Suite for BlockDiagonal class
%
%   Class BlockDiagonalTest
%
%   Example
%   BlockDiagonalTest
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-05-04,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.



%% Methods
methods (Test)
    function test_BlockDiagonal(testCase)
        % test for constructor

        % create the BlockMatrix object
        BD = BlockDiagonal({rand(2, 3), rand(2, 2), rand(1, 2)});
        
        siz = [5 7];
        testCase.verifyEqual(siz, size(BD));
    end
    
    function test_BlockDiagonal_Copy(testCase)
        % test for constructor

        % create the BlockMatrix object
        BD0 = BlockDiagonal({rand(2, 3), rand(2, 2), rand(1, 2)});
        BD = BlockDiagonal(BD0);
        
        siz = [5 7];
        testCase.verifyEqual(siz, size(BD));
    end
end

%% Test Functions for advanced computations
methods (Test)
    function test_blockNorm(testCase)
        diags = {[1 2;3 4], [5 6 7;8 9 10], [11;12]};
        BM = BlockDiagonal(diags);
        
        BM2 = blockNorm(BM);

        % result should be an instance of BlockDiagonal
        testCase.assertTrue(isa(BM2, 'BlockDiagonal'));
        % block size of result should match block size of input
        testCase.verifyEqual(blockSize(BM), blockSize(BM2));
        % result should be a scalar block matrix (all blocks are 1-by-1).
        testCase.verifyTrue(isScalarBlock(BM2));
    end
    
    function test_fapply(testCase)
        block1 = [1 2 3;4 5 6];
        block2 = [7 8 ; 9 10];
        block3 = [11 ; 12];
        BD = BlockDiagonal({block1, block2, block3});
        
        BD2 = fapply(@sqrt, BD);
        
        testCase.verifyEqual(blockSize(BD), blockSize(BD2));
        for i = 1:3
            testCase.verifyEqual(sqrt(diagonalBlock(BD, i)), diagonalBlock(BD2, i));
        end
    end
end


%% Test Functions for data access and manipulation
methods (Test)
    function test_diagonalBlock(testCase)
        block1 = [1 2 3;4 5 6];
        block2 = [7 8 ; 9 10];
        block3 = [11 ; 12];
        BD = BlockDiagonal({block1, block2, block3});
        
        testCase.verifyEqual(block1, diagonalBlock(BD, 1));
        testCase.verifyEqual(block2, diagonalBlock(BD, 2));
        testCase.verifyEqual(block3, diagonalBlock(BD, 3));
    end
        
    function test_diagonalBlocks(testCase)
        block1 = [1 2 3;4 5 6];
        block2 = [7 8 ; 9 10];
        block3 = [11 ; 12];
        BD = BlockDiagonal({block1, block2, block3});
        
        blocks = diagonalBlocks(BD);
        testCase.verifyEqual(block1, blocks{1});
        testCase.verifyEqual(block2, blocks{2});
        testCase.verifyEqual(block3, blocks{3});
    end
    
end % end methods

end % end classdef

