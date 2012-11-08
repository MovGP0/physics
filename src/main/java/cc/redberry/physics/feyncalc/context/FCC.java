package cc.redberry.physics.feyncalc.context;

import cc.redberry.core.indices.IndexType;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.utils.IntArray;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class FCC {
    private FCC() {
    }

    private static FeynCalcContext current() {
        return FeynCalcContext.current();
    }

    public static boolean isGammaMatrix(Tensor t) {
        return current().isGammaMatrix(t);
    }

    public static IndexType getDiracMatrixIndexType() {
        return current().getDiracMatrixIndexType();
    }

    public static IndexType getLorentzIndexType() {
        return current().getLorentzIndexType();
    }

    public static IntArray getGammaMatricesIndicesIds() {
        return current().getGammaMatrixIndicesIds();
    }

    public static String getGammaMatrixStringName() {
        return current().getGammaMatrixStringName();
    }
}
