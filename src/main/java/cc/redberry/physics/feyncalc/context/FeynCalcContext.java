package cc.redberry.physics.feyncalc.context;

import cc.redberry.core.context.CC;
import cc.redberry.core.context.NameDescriptor;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.indices.IndicesTypeStructure;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.utils.IntArray;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class FeynCalcContext {
    private final NameDescriptor gammaMatrixName;
    private final String gammaMatrixStringName;
    private final IntArray gammaMatrixIndicesIds;
    private final IndexType DiracMatrixIndexType;
    private final IndexType LorentzIndexType;


    public FeynCalcContext(FeynCalcSettings settings) {
        IndicesTypeStructure gammaMatrixTypesStruct = new IndicesTypeStructure(
                new byte[]{settings.LorentzType.getType(),
                        settings.DiracMatrixType.getType()},
                new int[]{1, 2});
        gammaMatrixName = CC.getNameManager().mapNameDescriptor(settings.DiracGammaName, gammaMatrixTypesStruct);
        DiracMatrixIndexType = settings.DiracMatrixType;
        LorentzIndexType = settings.LorentzType;
        this.gammaMatrixStringName = settings.DiracGammaName;
        int[] gammaMatrixIndicesIds = new int[2];
        if (LorentzIndexType.getType() < DiracMatrixIndexType.getType()) {
            gammaMatrixIndicesIds[0] = 1;
            gammaMatrixIndicesIds[1] = 2;
        } else {
            gammaMatrixIndicesIds[0] = 0;
            gammaMatrixIndicesIds[2] = 1;
        }
        this.gammaMatrixIndicesIds = new IntArray(gammaMatrixIndicesIds);
    }


    public boolean isGammaMatrix(Tensor tensor) {
        return (tensor instanceof SimpleTensor) && ((SimpleTensor) tensor).getName() == gammaMatrixName.getId();
    }

    public IndexType getDiracMatrixIndexType() {
        return DiracMatrixIndexType;
    }

    public IndexType getLorentzIndexType() {
        return LorentzIndexType;
    }

    public static FeynCalcContext current() {
        return FeynCalcContextManager.getCurrentContext();
    }

    public IntArray getGammaMatrixIndicesIds() {
        return gammaMatrixIndicesIds;
    }

    public String getGammaMatrixStringName() {
        return gammaMatrixStringName;
    }
}
