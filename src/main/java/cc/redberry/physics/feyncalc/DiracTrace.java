package cc.redberry.physics.feyncalc;

import cc.redberry.core.context.CC;
import cc.redberry.core.context.NameDescriptor;
import cc.redberry.core.indexgenerator.IndexGenerator;
import cc.redberry.core.indices.*;
import cc.redberry.core.number.Complex;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterator.TensorLastIterator;
import cc.redberry.core.transformations.ContractIndices;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.transformations.expand.Expand;
import cc.redberry.core.transformations.expand.ExpandAll;
import cc.redberry.core.transformations.substitutions.Substitution;

import static cc.redberry.core.indices.IndicesUtils.setType;
import static cc.redberry.core.tensor.Tensors.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class DiracTrace implements Transformation {
    private final SimpleTensor gammaMatrix;

    public DiracTrace(SimpleTensor gammaMatrix) {
        this.gammaMatrix = gammaMatrix;
    }

    @Override
    public Tensor transform(Tensor t) {
        return trace(t, gammaMatrix);
    }

    public static Tensor trace(Tensor tensor) {
        return trace(tensor, Tensors.parseSimple("G^a'_b'a"));
    }

    public static Tensor trace(Tensor tensor, SimpleTensor gammaMatrix) {
        IndexType[] types = TraceUtils.extractTypesFromMatrix(gammaMatrix);
        //todo check for contains gammas
        tensor = ExpandAll.expandAll(tensor, ContractIndices.ContractIndices);
        tensor = ContractIndices.contract(tensor);
        TensorLastIterator iterator = new TensorLastIterator(tensor);
        Tensor current;
        while ((current = iterator.next()) != null) {
            if (current instanceof SimpleTensor && ((SimpleTensor) current).getName() == gammaMatrix.getName() && current.getIndices().getFree().size(types[1]) == 0) {
                iterator.set(Complex.ZERO);
            } else if (current instanceof Product) {
                if (current.getIndices().getFree().size(types[1]) != 0)
                    continue;

                int gammasCount = 0;
                for (Tensor t : current) {
                    if (t instanceof SimpleTensor && ((SimpleTensor) t).getName() == gammaMatrix.getName())
                        ++gammasCount;
                }
                if (gammasCount == 0)
                    continue;

                if (gammasCount % 2 == 1)
                    iterator.set(Complex.ZERO);

                current = createSubstitution(gammaMatrix.getName(),
                        gammasCount, types[0], types[1]).transform(current);
                current = Expand.expand(current, ContractIndices.ContractIndices);
                current = ContractIndices.contract(current);
                iterator.set(current);
            }
        }
        return iterator.result();
    }

    private static Substitution createSubstitution(int gammaName, int gammasCount, IndexType metricType, IndexType matrixType) {
        Tensor[] gammas = new Tensor[gammasCount];
        IndexGenerator generator = new IndexGenerator();
        int firstUpper, u = firstUpper = generator.generate(matrixType), i;
        for (i = 0; i < gammasCount; ++i) {
            gammas[i] = Tensors.simpleTensor(gammaName,
                    IndicesFactory.createSimple(null,
                            u | 0x80000000,
                            i == gammasCount - 1 ? firstUpper : (u = generator.generate(matrixType)),
                            generator.generate(metricType)));

        }
        return new Substitution(Tensors.multiply(gammas), traceOfArray(gammas, metricType));
    }

    private static Substitution gamma5(int gammaName, int gamma5Name, byte metricType, byte matrixType) {
        Expression[] gamma5 = new Expression[3];

        //G5*G_m = -G_m*G5
        SimpleTensor a, b;
        Tensor[] lrhs = new Tensor[2];
        int[] gammaNames = {gamma5Name, gammaName};
        SimpleIndices aIndices = IndicesFactory.createSimple(null,
                0x80000000 | setType(matrixType, 0),
                setType(matrixType, 1),
                setType(metricType, 0)),
                bIndices = IndicesFactory.createSimple(null,
                        0x80000000 | setType(matrixType, 1),
                        setType(matrixType, 2),
                        setType(metricType, 1));
        for (int i = 0; i < 2; ++i) {
            a = simpleTensor(gammaNames[i], aIndices);
            b = simpleTensor(gammaNames[1 - i], bIndices);
            lrhs[i] = multiply(a, b);
        }
        gamma5[0] = expression(lrhs[0], lrhs[1]);

        //G5*G5 = 1
        lrhs[0] = multiply(simpleTensor(gamma5Name, aIndices), simpleTensor(gamma5Name, bIndices));
        lrhs[1] = createKronecker(0x80000000 | setType(matrixType, 0), 0x80000000 | setType(matrixType, 2));
        gamma5[1] = expression(lrhs[0], lrhs[1]);

        //Tr[G5] = 0
        gamma5[2] = expression(simpleTensor(gamma5Name,
                IndicesFactory.createSimple(null, 0x80000000 | setType(matrixType, 0), 0)), Complex.ZERO);

        return new Substitution(gamma5);
    }

    static Tensor traceOfArray(Tensor[] product, IndexType metricType) {
        if (product.length == 1)
            return Complex.ZERO;
        if (product.length == 2)
            return multiply(Complex.FOUR,
                    createMetricOrKronecker(product[0].getIndices().get(metricType, 0),
                            product[1].getIndices().get(metricType, 0)));
        if (product.length % 2 != 0)
            return Complex.ZERO;
        SumBuilder sb = new SumBuilder();
        Tensor temp;
        for (int i = 0; i < product.length - 1; ++i) {
            temp = multiply(Complex.TWO,
                    createMetricOrKronecker(product[i].getIndices().get(metricType, 0),
                            product[i + 1].getIndices().get(metricType, 0)),
                    traceOfArray(subArray(product, i, i + 1), metricType));
            if (i % 2 != 0)
                temp = negate(temp);
            sb.put(temp);
            swap(product, i, i + 1);
        }
        return multiply(Complex.ONE_HALF, sb.build());
    }

    static Tensor traceOfArray5(Tensor[] product, IndexType metricType, int epsilonName) {
        if (product.length < 5 || product.length % 2 == 0)
            return Complex.ZERO;
        if (product.length == 5) {
            int[] indices = new int[4];
            for (int i = 0; i < 4; ++i)
                indices[i] = product[i].getIndices().get(metricType, 0);
            return multiply(Complex.FOUR.negate(), Complex.IMAGE_ONE, simpleTensor(epsilonName, IndicesFactory.createSimple(null, indices)));
        }
        return null;
    }

    private static Tensor[] subArray(Tensor[] array, int a, int b) {
        Tensor[] result = new Tensor[array.length - 2];
        int k = 0;
        for (int i = 0; i < array.length; ++i) {
            if (i == a || i == b)
                continue;
            result[k++] = array[i];
        }
        return result;
    }

    private static void swap(Tensor[] array, int a, int b) {
        Tensor temp = array[a];
        array[a] = array[b];
        array[b] = temp;
    }
}
