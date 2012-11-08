package cc.redberry.physics.feyncalc;

import cc.redberry.core.indices.IndexType;
import cc.redberry.core.indices.IndicesUtils;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.tensor.Tensor;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class MatrixUtils {
    private MatrixUtils() {
    }

    public static boolean isMatrixOrVectorOrCovector(Tensor t, IndexType matrixType) {
        int size;
        return t instanceof SimpleTensor
                && ((size = t.getIndices().size(matrixType)) == 2
                || size == 1);
    }

    public static boolean isVectorOrCovector(Tensor t, IndexType matrixType) {
        return t instanceof SimpleTensor
                && t.getIndices().size(matrixType) == 1;
    }

    public static boolean isMatrix(Tensor t, IndexType matrixType) {
        return t instanceof SimpleTensor
                && t.getIndices().size(matrixType) == 2;
    }

    public static boolean isVector(Tensor t, IndexType matrixType) {
        return t instanceof SimpleTensor
                && t.getIndices().size(matrixType) == 1
                && IndicesUtils.getState(t.getIndices().get(matrixType, 0));
    }

    public static boolean isCovector(Tensor t, IndexType matrixType) {
        return t instanceof SimpleTensor
                && t.getIndices().size(matrixType) == 1
                && !IndicesUtils.getState(t.getIndices().get(matrixType, 0));
    }
}
