package cc.redberry.physics.feyncalc.context;

import cc.redberry.core.indices.IndexType;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class DefaultFeynCalcSettings {
    private DefaultFeynCalcSettings() {
    }

    public static FeynCalcSettings create() {
        String defaultDiracGamma = "G";
        IndexType defaultLorentzType = IndexType.LatinLower;
        IndexType defaultDiracMatrixType = IndexType.LatinLower1;
        return new FeynCalcSettings(defaultDiracGamma, defaultLorentzType, defaultDiracMatrixType);
    }
}
