package cc.redberry.physics.feyncalc;

import cc.redberry.core.indexmapping.IndexMapping;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.indices.Indices;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.TensorLastIterator;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.utils.ArraysUtils;
import cc.redberry.physics.feyncalc.context.FCC;

import java.util.Arrays;

import static cc.redberry.core.indices.IndicesUtils.inverseIndexState;
import static cc.redberry.core.tensor.FullContractionsStructure.getFromIndexId;
import static cc.redberry.core.tensor.FullContractionsStructure.getToTensorIndex;
import static cc.redberry.physics.feyncalc.MatrixUtils.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class InverseProductOfMatrices implements Transformation {
    final IndexType matrixType;

    public InverseProductOfMatrices(IndexType matrixType) {
        this.matrixType = matrixType;
    }

    @Override
    public Tensor transform(Tensor t) {
        return inverseProductsOfMatrices(t, matrixType);
    }

    public static Tensor inverseProductsOfMatrices(Tensor t, IndexType matrixType) {
        TensorLastIterator iterator = new TensorLastIterator(t);
        Tensor c;
        while ((c = iterator.next()) != null)
            if (c instanceof Product)
                iterator.set(inverseProductOfMatrices((Product) c, matrixType));

        return iterator.result();
    }


    public static Product inverseProductOfMatrices(Product p, IndexType matrixType) {
        boolean needToTransform = false, containsVectors = false;

        for (Tensor m : p)
            if (isMatrix(m, matrixType)) {
                needToTransform = true;
                break;
            } else if (isVectorOrCovector(m, matrixType))
                containsVectors = true;

        if (!needToTransform)
            return p;
        Indices freeIndices = p.getIndices().getFree().getOfType(matrixType);
        if (!containsVectors && freeIndices.size() == 0)
            return p;
        int[] matrixIndices;
        if (!containsVectors)
            matrixIndices = freeIndices.getAllIndices().copy();
        else {
            matrixIndices = new int[2];
            for (Tensor m : p) {
                if (isCovector(m, matrixType))
                    matrixIndices[0] = inverseIndexState(m.getIndices().get(matrixType, 0));
                if (isVector(m, matrixType))
                    matrixIndices[1] = inverseIndexState(m.getIndices().get(matrixType, 0));
            }
        }
        //indices are sorted, so first is upper
        int firstUpper = matrixIndices[0];

        ProductContent content = p.getContent();
        FullContractionsStructure fs = content.getFullContractionsStructure();

        final Tensor[] data = content.getDataCopy();
        //search the first upper matrix index
        int currentIndex = 0;
        Tensor temp;
        for (; currentIndex < data.length; ++currentIndex) {
            temp = data[currentIndex];
            if (!isMatrix(temp, matrixType))
                continue;
            if (temp.getIndices().get(matrixType, 0) == firstUpper)
                break;
        }

        assert currentIndex != content.size();

        //this is the first matrix in the initial product
        //and the last in the result
        Tensor previousMatrix = data[currentIndex], currentMatrix, nextMatrix;
        //index mapping
        int[][] fromTo = new int[2][2];

        long fromContractions[];
        int nextIndex = -1;
        currentMatrix = previousMatrix;
        while (true) {//traversing through the product graph

            nextMatrix = null;
            fromContractions = fs.contractions[currentIndex];
            for (long from : fromContractions) {
                if (getFromIndexId(from) != FCC.getGammaMatricesIndicesIds().get(1))
                    continue;
                nextIndex = getToTensorIndex(from);
                if (nextIndex == -1)
                    continue;
                nextMatrix = data[nextIndex];
                if (nextMatrix != previousMatrix && isMatrix(nextMatrix, matrixType))
                    break;
                else nextMatrix = null;
            }
            fromTo[0] = currentMatrix.getIndices().getOfType(matrixType).getAllIndices().copy();

            if (currentMatrix == previousMatrix)
                //rewriting first matrix
                fromTo[1] = renameMatrixIndices(
                        currentMatrix,
                        matrixType,
                        new int[]{inverseIndexState(fromTo[0][1]), matrixIndices[1]});
            else if (nextMatrix == null)
                //rewriting the last matrix
                fromTo[1] = renameMatrixIndices(
                        currentMatrix,
                        matrixType,
                        new int[]{matrixIndices[0], inverseIndexState(fromTo[0][0])});
            else
                //rewriting current matrix
                fromTo[1] = renameMatrixIndices(
                        currentMatrix,
                        matrixType,
                        new int[]{inverseIndexState(fromTo[0][1]), inverseIndexState(fromTo[0][0])});

            currentMatrix =
                    ApplyIndexMapping.applyIndexMapping(currentMatrix,
                            renameMatrixIndices(currentMatrix, matrixType, fromTo[0]),
                            fromTo[1], new int[0]);
            data[currentIndex] = currentMatrix;

            if (nextMatrix == null)
                break;

            currentIndex = nextIndex;
            previousMatrix = currentMatrix;
            currentMatrix = nextMatrix;
        }
        return (Product) Tensors.multiply(ArraysUtils.addAll(data, p.getIndexlessSubProduct()));
    }

    private static int[] renameMatrixIndices(Tensor tensor, IndexType matrixType, int[] newMatrixIndices) {
        Indices free = tensor.getIndices().getFree();
        if (free.size() == newMatrixIndices.length)
            return newMatrixIndices;
        int[] from = free.getOfType(matrixType).getAllIndices().copy();
        assert from.length == newMatrixIndices.length;
        ArraysUtils.quickSort(from, newMatrixIndices);
        return free.applyIndexMapping(new IndexMapper(from, newMatrixIndices)).getAllIndices().copy();
    }

    private static final class IndexMapper implements IndexMapping {
        final int[] from, to;

        private IndexMapper(int[] from, int[] to) {
            this.from = from;
            this.to = to;
        }

        @Override
        public int map(int index) {
            int i;
            if ((i = Arrays.binarySearch(from, index)) >= 0)
                return to[i];
            return index;
        }
    }

}
