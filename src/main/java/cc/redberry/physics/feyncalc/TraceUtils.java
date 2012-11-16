package cc.redberry.physics.feyncalc;

import cc.redberry.core.number.Complex;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.tensor.SumBuilder;
import cc.redberry.core.tensor.Tensor;

import static cc.redberry.core.tensor.Tensors.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class TraceUtils {
    private TraceUtils() {
    }

    public static final class Commutator {
        final boolean isAntiCommutator;
        final SimpleTensor A, B;
        final Tensor traceOfA, traceOfB;
        Tensor C;

        /**
         * Commutato container [A, B] = C
         *
         * @param A              matrix
         * @param B              matrix
         * @param antiCommutator is specified then the relation treats as AB+BA = C
         * @param C              matrix
         */
        public Commutator(SimpleTensor A, SimpleTensor B, Tensor C,
                          boolean antiCommutator,
                          Tensor traceOfA, Tensor traceOfB) {
            this.isAntiCommutator = antiCommutator;
            this.A = A;
            this.B = B;
            this.C = C;
            this.traceOfA = traceOfA;
            this.traceOfB = traceOfB;
        }


        final boolean applicableTo(SimpleTensor A, SimpleTensor B) {
            return (A.getName() == this.A.getName() && B.getName() == this.B.getName())
                    || (A.getName() == this.B.getName() && B.getName() == this.A.getName());
        }

        final Tensor traceOf(SimpleTensor tensor) {
            SimpleTensor match = null;
            if (tensor.getName() == A.getName())
                match = A;
            else if (tensor.getName() == B.getName())
                match = B;
            if (match == null)
                return null;

        }
    }

    private static void checkCommutatorConsistency(Expression commutator) {
    }

    public Tensor traceOfArray(Tensor[] array, Commutator[] commutators) {
        if (array.length == 1)
            return Complex.ZERO;
        if (array.length == 2)
            return multiply(Complex.FOUR,
                    createMetricOrKronecker(array[0].getIndices().get(metricType, 0),
                            array[1].getIndices().get(metricType, 0)));
        if (array.length % 2 != 0)
            return Complex.ZERO;
        SumBuilder sb = new SumBuilder();
        Tensor temp;
        for (int i = 0; i < array.length - 1; ++i) {
            temp = multiply(Complex.TWO,
                    createMetricOrKronecker(array[i].getIndices().get(metricType, 0),
                            array[i + 1].getIndices().get(metricType, 0)),
                    traceOfArray(subArray(array, i, i + 1), metricType));
            if (i % 2 != 0)
                temp = negate(temp);
            sb.put(temp);
            swap(array, i, i + 1);
        }
        return multiply(Complex.ONE_HALF, sb.build());
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
