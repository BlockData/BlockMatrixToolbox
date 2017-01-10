classdef BlockMatrixTest < matlab.unittest.TestCase
%BLOCKMATRIXTEST Test Suite for BlockMatrix class
%
%   Class BlockMatrixTest
%
%   Example
%   run(BlockMatrixTest);
%
%   See also
%   BlockMatrix

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-02-27,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.

%% Test Functions for constructor
methods (Test)
    
    function test_BlockMatrix(testCase)
        % test for constructor

        % create the BlockMatrix object
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix(data, parts);

        testCase.verifyEqual(isempty(BM), false);
    end
end


%% Test Functions for advanced computations
methods (Test)
    function test_blockNorm(testCase)
        data = reshape(1:28, [7 4])';
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);
        
        BM2 = blockNorm(BM);

        % result should be an instance of BlockMatrix
        testCase.assertTrue(isa(BM2, 'BlockMatrix'));
        % block size of result should match block size of input
        testCase.verifyEqual(blockSize(BM), blockSize(BM2));
        % result should be a scalar block matrix (all blocks are 1-by-1).
        testCase.verifyTrue(isScalarBlock(BM2));
    end
    
    function test_fapply(testCase)
        data = reshape(1:28, [7 4])';
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);
        
        BM2 = fapply(@sqrt, BM);
        
        testCase.verifyEqual(blockSize(BM), blockSize(BM2));
        testCase.verifyEqual(sqrt(data), BM2.data);
    end
end


%% Test Functions for static methods
methods (Test)
    function test_oneBlock_matrix(testCase)
        data = [1 2 3;4 5 6];
        BM = BlockMatrix.oneBlock(data);

        testCase.verifyEqual(data, getMatrix(BM));
        testCase.verifyEqual(1, blockSize(BM, 1));
        testCase.verifyEqual(1, blockSize(BM, 2));
    end
    
    function test_oneBlock_blockMatrix(testCase)
        data = reshape((1:28)', [7 4])';
        BM0 = BlockMatrix.create(data, [2 2], [2 3 2]);
        
        BM = BlockMatrix.oneBlock(BM0);

        testCase.verifyEqual(data, getMatrix(BM));
        testCase.verifyEqual(1, blockSize(BM, 1));
        testCase.verifyEqual(1, blockSize(BM, 2));
    end
    
    function test_scalarBlock_matrix(testCase)
        data = reshape((1:28)', [7 4])';
        BM = BlockMatrix.scalarBlock(data);

        testCase.verifyEqual(data, getMatrix(BM));
        testCase.verifyEqual(4, blockSize(BM, 1));
        testCase.verifyEqual(7, blockSize(BM, 2));
    end
    
    function test_scalarBlock_blockMatrix(testCase)
        data = reshape((1:28)', [7 4])';
        BM0 = BlockMatrix.create(data, [2 2], [2 3 2]);
        
        BM = BlockMatrix.scalarBlock(BM0);

        testCase.verifyEqual(data, getMatrix(BM));
        testCase.verifyEqual(4, blockSize(BM, 1));
        testCase.verifyEqual(7, blockSize(BM, 2));
    end
    
    function test_uniformBlock_matrix(testCase)
        data = reshape((1:24)', [6 4])';
        
        BM = BlockMatrix.uniformBlocks(data, [2 2]);
        
        testCase.verifyEqual(data, getMatrix(BM));
        testCase.verifyEqual(2, blockSize(BM, 1));
        testCase.verifyEqual(3, blockSize(BM, 2));
    end
    
    function test_uniformBlock_blockMatrix(testCase)
        data = reshape((1:24)', [6 4])';
        BM0 = BlockMatrix.create(data, [1 2 1], [2 1 3]);
        
        BM = BlockMatrix.uniformBlocks(BM0, [2 2]);
        
        testCase.verifyEqual(data, getMatrix(BM));
        testCase.verifyEqual(2, blockSize(BM, 1));
        testCase.verifyEqual(3, blockSize(BM, 2));
    end
    
end


%% Test Functions for block products
methods (Test)
    function test_blockProduct_ss(testCase)
        data = reshape(1:16, [8 2]);
        scalar = 3;
        
        A = BlockMatrix(data, {[4 4], [1 1]});
        lambda = BlockMatrix.oneBlock(scalar);
        AA = blockProduct_ss(lambda, A);
        
        testCase.verifyEqual(blockSize(A), blockSize(AA));
        testCase.verifyEqual(data * scalar, getMatrix(AA));
    end
    
    function test_blockProduct_sh(testCase)
        dataA = [1 2 3;3 2 1];
        dataB = reshape(1:36, [4 9]);
        
        A = BlockMatrix.oneBlock(dataA);
        B = BlockMatrix.create(dataB, {[2 2], [3 3 3]});
        X = blockProduct_sh(A, B);
        
        testCase.verifyEqual(blockSize(B), blockSize(X));
        exp = dataB .* repmat(dataA, 2, 3);
        testCase.verifyEqual(exp, getMatrix(X));
    end
    
    function test_blockProduct_su(testCase)
        blockA = [1 2 3;3 2 1];
        A = BlockMatrix.oneBlock(blockA);
        B = BlockMatrix.create(reshape(1:36, [6 6]), {[3 3], [2 2 2]});
        X = blockProduct_su(A, B);
        
        % verify the dimensions of the result BlockMatrix
        testCase.verifyEqual([4 6], size(X));
        testCase.verifyEqual(blockSize(B), blockSize(X));
        testCase.verifyEqual(IntegerPartition([2 2]), blockDimensions(X, 1));
        testCase.verifyEqual(IntegerPartition([2 2 2]), blockDimensions(X, 2));
        
        % verify that each block correspond to usual product of input
        % blocks
        for iBlock = 1:blockSize(B, 1)
            for jBlock = 1:blockSize(B, 2)
                block = blockA * B{iBlock, jBlock};
                testCase.verifyEqual(block, X{iBlock, jBlock});
            end
        end            
    end
    
    function test_blockProduct_sk(testCase)
        A = BlockMatrix.oneBlock(reshape(1:4, [2 2]));
        B = BlockMatrix.create(reshape(1:36, [6 6]), {[3 3], [2 2 2]});
        X = blockProduct_sk(A, B);
        
        testCase.verifyEqual(size(A) .* size(B), size(X));
        testCase.verifyEqual(blockSize(B), blockSize(X));        
    end
    
    function test_blockProduct_hs(testCase)
        A = BlockMatrix.scalarBlock(reshape(1:4, [2 2]));
        B = BlockMatrix.create(reshape(1:16, [4 4]), {[2 2], [2 2]});
        X = blockProduct_hs(A, B);
 
        testCase.verifyEqual(size(B), size(X));
        testCase.verifyEqual(blockSize(B), blockSize(X));
        
        for iBlock = 1:blockSize(B, 1)
            for jBlock = 1:blockSize(B, 2)
                block = A(iBlock, jBlock) * B{iBlock, jBlock};
                testCase.verifyEqual(block, X{iBlock, jBlock});
            end
        end            
    end
    
    function test_blockProduct_hh(testCase)
        dataA = reshape(1:36, [6 6]);
        dims = BlockDimensions({[3 3], [2 2 2]});
        A = BlockMatrix.create(dataA, dims);
        dataB = magic(6);
        B = BlockMatrix.create(dataB, dims);
        X = blockProduct_hh(A, B);

        testCase.verifyEqual(size(B), size(X));
        testCase.verifyEqual(blockSize(B), blockSize(X));
        exp = dataA .* dataB;
        testCase.verifyEqual(exp, getMatrix(X));
    end
    
    function test_blockProduct_hu(testCase)
        A = BlockMatrix.create(reshape(1:20, [5 4]), {[3 2], [2 2]});
        B = BlockMatrix.create(reshape(1:12, [4 3]), {[2 2], [2 1]});
        X = blockProduct_hu(A, B);
        
        % verify validity of sizes
        testCase.verifyEqual(size(A, 1), size(X, 1));
        testCase.verifyEqual(size(B, 2), size(X, 2));
        
        % verify validity of block-sizes
        dimsA = blockDimensions(A);
        dimsB = blockDimensions(B);
        dimsX = blockDimensions(X);
        testCase.verifyEqual(dimsA{1}, dimsX{1});
        testCase.verifyEqual(dimsB{2}, dimsX{2});
        
        for iBlock = 1:blockSize(A, 1)
            for jBlock = 1:blockSize(A, 2)
                block = A{iBlock, jBlock} * B{iBlock, jBlock};
                testCase.verifyEqual(block, X{iBlock, jBlock});
            end
        end
    end
    
    function test_blockProduct_hk(testCase)
        A = BlockMatrix.create(reshape(1:12, [3 4])', [2 2], [1 1 1]);
        B = BlockMatrix.create(reshape(1:36, [4 9]), [2 2], [3 3 3]);
        X = blockProduct_hk(A, B);
       
        testCase.verifyEqual(blockSize(A, 1), blockSize(X, 1));
        testCase.verifyEqual(blockSize(B, 2), blockSize(X, 2));
        
        dimsA = blockDimensions(A);
        dimsB = blockDimensions(B);
        dimsX = blockDimensions(X);
        testCase.verifyEqual(dimsA{1} .* dimsB{1}, dimsX{1});
        testCase.verifyEqual(dimsA{2} .* dimsB{2}, dimsX{2});
        
    end
    
    function test_blockProduct_us(testCase)
        A = BlockMatrix.scalarBlock([1 2 3; 3 2 1]);
        dataB = reshape(1:36, [4 9])';
        B = BlockMatrix.create(dataB, {[3 3 3], [2 2]});
        X = blockProduct_us(A, B);
        
        testCase.verifyEqual(blockSize(A, 1), blockSize(X, 1));
        testCase.verifyEqual(blockSize(B, 2), blockSize(X, 2));
    end
    
    function test_blockProduct_uh(testCase)
        A = BlockMatrix.create(reshape(1:36, [9 4])', {[2 2], [3 3 3]});
        B = BlockMatrix.create(reshape(1:36, [6 6]), {[2 2 2], [3 3]});
        X = blockProduct_uh(A, B);
        
        testCase.verifyEqual(blockSize(A, 1), blockSize(X, 1));
        testCase.verifyEqual(blockSize(B, 2), blockSize(X, 2));
    end
    
    function test_blockProduct_uu(testCase)
        A = BlockMatrix.create(reshape(1:36, [6 6]), {[2 2 2], [3 3]});
        B = BlockMatrix.create(reshape(1:36, [6 6]), {[3 3], [2 1 2 1]});
        X = blockProduct_uu(A, B);
        
        testCase.verifyEqual(blockSize(A, 1), blockSize(X, 1));
        testCase.verifyEqual(blockSize(B, 2), blockSize(X, 2));
    end
    
    function test_blockProduct_uk(testCase)
        A = BlockMatrix.create(reshape(1:12, [3 4]), 3, [1 1 1 1 ]);
        B = BlockMatrix.create(reshape(1:12, [4 3]), [1 1 1 1], [1 2]);
        X = blockProduct_uk(A, B);
        
        testCase.verifyEqual(blockSize(A, 1), blockSize(X, 1));
        testCase.verifyEqual(blockSize(B, 2), blockSize(X, 2));
    end
    
    function test_blockProduct_ks(testCase)
        A = BlockMatrix.scalarBlock(reshape(1:4, [2 2]));
        B = BlockMatrix.oneBlock(reshape(1:9, [3 3]));
        X = blockProduct_ks(A, B);
       
        testCase.verifyEqual(blockSize(A), blockSize(X));
        testCase.verifyEqual(size(A) .* size(B), size(X));
    end
    
    function test_blockProduct_kh(testCase)
        A = BlockMatrix.create(reshape(1:36, [9 4]), {[3 3 3], [2 2]});
        B = BlockMatrix.oneBlock(2*ones(3,2));
        X = blockProduct_kh(A, B);
        
        testCase.verifyEqual(blockSize(A), blockSize(X));
    end
    
   function test_blockProduct_ku(testCase)
       A = BlockMatrix.create(reshape(1:36, [6 6]), {[2 2 2], [3 3]});
       B = BlockMatrix.oneBlock(ones(3,2));
       X = blockProduct_ku(A, B);
       
       testCase.verifyEqual(blockSize(A), blockSize(X));
   end
end % end of blockProduct test cases


%% Test Functions for basic array manipulation
methods (Test)
    function test_size(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);
        
        siz = size(BM);
        testCase.verifyEqual(siz(1), 4);
        testCase.verifyEqual(siz(2), 7);
    end
    
    function test_size_two_outputs(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);
        
        [siz1, siz2] = size(BM);
        testCase.verifyEqual(siz1, 4);
        testCase.verifyEqual(siz2, 7);
    end
    
    function test_blockSize(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);
        
        siz = blockSize(BM);
        testCase.verifyEqual(siz(1), 2);
        testCase.verifyEqual(siz(2), 3);
    end
    
    function test_blockSize_two_outputs(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);
        
        [siz1, siz2] = blockSize(BM);
        testCase.verifyEqual(siz1, 2);
        testCase.verifyEqual(siz2, 3);
    end
    
end


%% Test Functions for testing block matrix nature
methods (Test)
    function test_isOneBlock_true(testCase)
       data = reshape(1:28, [7 4])';
       BM = BlockMatrix.create(data, {4, 7});
       testCase.verifyTrue(isOneBlock(BM));
    end
    
    function test_isOneBlock_false(testCase)
        data = reshape(1:28, [7 4])';
        BM = BlockMatrix.create(data, {[2 2], [2 3 2]});
        testCase.verifyFalse(isOneBlock(BM));
    end
    
    function test_isScalarBlock_true(testCase)
        data = reshape(1:28, [7 4])';
        BM = BlockMatrix.create(data, {ones(1,4), ones(1,7)});
        testCase.verifyTrue(isScalarBlock(BM));
    end
    
    function test_isScalarBlock_false(testCase)
        data = reshape(1:28, [7 4])';
        BM = BlockMatrix.create(data, {[2 2], [2 3 2]});
        testCase.verifyFalse(isScalarBlock(BM));
    end
    
    function test_isUniformBlock_true(testCase)
       data = reshape(1:24, [6 4])';
       BM = BlockMatrix.create(data, {[2 2], [3 3]});
       testCase.verifyTrue(isUniformBlock(BM));
    end
    
    function test_isUniformBlock_false(testCase)
       data = reshape(1:24, [6 4])';
       BM = BlockMatrix.create(data, {[2 2], [4 2]});
       testCase.verifyFalse(isUniformBlock(BM));
    end

    function test_isVectorBlock_true(testCase)
       data = reshape(1:24, [6 4])';
       BM = BlockMatrix.create(data, {4, [2 2 2]});
       testCase.verifyTrue(isVectorBlock(BM));
    end
    
    function test_isVectorBlock_false(testCase)
       data = reshape(1:24, [6 4])';
       BM = BlockMatrix.create(data, {[2 2], [4 2]});
       testCase.verifyFalse(isVectorBlock(BM));
    end

end


%% Test Functions arithmetic operations
methods (Test)
    function test_log_exp(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);

        res = log(exp(BM));
        
        testCase.verifyEqual(data, getMatrix(res));
    end
end

%% Test Functions for basic array manipulation
methods (Test)

    function test_transpose(testCase)
        % test the transpose method
        
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);
        
        BM2 = BM';
        
        dim1 = size(BM2, 1);
        dim2 = size(BM2, 2);
        testCase.verifyEqual(dim1, 7);
        testCase.verifyEqual(dim2, 4);
    end
    
    function test_cat_dir1(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);

        BM2 = cat(1, BM, BM);

        dim1 = size(BM2, 1);
        dim2 = size(BM2, 2);
        testCase.verifyEqual(dim1, 8);
        testCase.verifyEqual(dim2, 7);
    end
    
    function test_cat_dir2(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);

        BM2 = cat(2, BM, BM);

        dim1 = size(BM2, 1);
        dim2 = size(BM2, 2);
        testCase.verifyEqual(dim1, 4);
        testCase.verifyEqual(dim2, 14);
    end
    
    function test_horzcat(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);

        BM2 = [BM BM];

        dim1 = size(BM2, 1);
        dim2 = size(BM2, 2);
        testCase.verifyEqual(dim1, 4);
        testCase.verifyEqual(dim2, 14);
    end
    
    function test_vertcat(testCase)
        data = reshape(1:28, [4 7]);
        parts = {[2 2], [2 3 2]};
        BM = BlockMatrix.create(data, parts);

        BM2 = [BM ; BM];

        dim1 = size(BM2, 1);
        dim2 = size(BM2, 2);
        testCase.verifyEqual(dim1, 8);
        testCase.verifyEqual(dim2, 7);
    end
end


%% Test Functions for overloadig arithmetic operations
methods (Test)
    
    function test_times(testCase)
        data1 = reshape(1:28, [4 7]);
        parts1 = {[2 2], [2 3 2]};
        BM1 = BlockMatrix.create(data1, parts1);
        
        data2 = reshape(1:35, [7 5]);
        parts2 = {[2 3 2], [2 1 2]};
        BM2 = BlockMatrix.create(data2, parts2);
        
        data3 = data1 * data2;
        BM3 = BM1 * BM2;
        testCase.verifyEqual(data3, BM3.data, 'AbsTol', .1);
    end

end

%% Test Functions for array indexing
methods (Test)
    
    function test_subsref_parens(testCase)
        BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        res = BM(2, 3);
        testCase.verifyEqual(10, res, 'AbsTol', .1);
    end
    
    function test_subsref_braces(testCase)
        BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        res = BM{2, 3};
        exp = BlockMatrix.oneBlock([20 21; 27 28]);
        testCase.verifyEqual(exp, res, 'AbsTol', .1);
    end

    function test_subsref_braces_recurse(testCase)
        BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        res = BM{2, 2}(1, 2);
        exp = 18;
        testCase.verifyEqual(exp, res, 'AbsTol', .1);
    end

    function test_subsasgn_parens(testCase)
        BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        BM(2, 3) = 4;
        testCase.verifyEqual(4, BM(2, 3));
    end
    
    function test_subsasgn_braces(testCase)
        BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 2 3]);
        BM{2, 3} = [1 2 3;4 5 6];
        testCase.verifyEqual(4, BM(4, 5));
    end

    function test_subsasgn_chained_scalarToScalar(testCase)
        BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        BM{1,2}(2,3) = 10;
        val = BM(2, 5);
        testCase.verifyEqual(10, val, 'AbsTol', .1);
    end
    
    function test_subsasgn_chained_matrixToMatrix(testCase)
        BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        BM{1,2}(1:2, 1:3) = [10 11 12;21 22 23];
        val = BM(2, 5);
        testCase.verifyEqual(23, val, 'AbsTol', .1);
    end
    
end

end % end classdef

