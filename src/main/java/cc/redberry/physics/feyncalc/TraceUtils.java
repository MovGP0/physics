package cc.redberry.physics.feyncalc;

import cc.redberry.core.context.CC;
import cc.redberry.core.context.NameDescriptor;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.indices.IndicesTypeStructure;
import cc.redberry.core.tensor.SimpleTensor;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
final class TraceUtils {
    static final IndexType[] extractTypesFromMatrix(SimpleTensor matrix) {
        if (matrix.getIndices().size() != 3)
            throw new IllegalArgumentException("Not a SU(N) matrix: " + matrix + ".");
        NameDescriptor descriptor = CC.getNameDescriptor(matrix.getName());
        IndicesTypeStructure typeStructure = descriptor.getIndicesTypeStructure();
        byte metricType = -1, matrixType = -1;
        int typeCount;
        for (byte type = 0; type < IndexType.TYPES_COUNT; ++type) {
            typeCount = typeStructure.typeCount(type);
            if (typeCount == 0)
                continue;
            else if (typeCount == 2) {
                if (matrixType != -1)
                    throw new IllegalArgumentException("Not a SU(N) matrix: " + matrix + ".");
                matrixType = type;
                if (CC.isMetric(matrixType))
                    throw new IllegalArgumentException("Not a SU(N) matrix: " + matrix + ".");
            } else if (typeCount == 1) {
                if (metricType != -1)
                    throw new IllegalArgumentException("Not a SU(N) matrix: " + matrix + ".");
                metricType = type;
                if (!CC.isMetric(metricType))
                    throw new IllegalArgumentException("Not a SU(N) matrix: " + matrix + ".");
            } else
                throw new IllegalArgumentException("Not a SU(N) matrix: " + matrix + ".");
        }
        return new IndexType[]{IndexType.getType(metricType), IndexType.getType(matrixType)};
    }
}
