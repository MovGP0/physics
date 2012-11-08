package cc.redberry.physics.feyncalc.context;

import cc.redberry.core.indices.IndexType;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class FeynCalcSettings {
    public String DiracGammaName;
    public IndexType LorentzType;
    public IndexType DiracMatrixType;

    public FeynCalcSettings(String diracGammaName, IndexType lorentzType, IndexType diracMatrixType) {
        DiracGammaName = diracGammaName;
        LorentzType = lorentzType;
        DiracMatrixType = diracMatrixType;
    }

}
